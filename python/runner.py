"""
runner.py
=========
Command-line interface for the inelastic-rebound-md project.

Subcommands:
    equilibrate   Run only the surface equilibration phase.
    impact        Run a single impact simulation (requires restart).
    sweep         Run a parametric sweep over v0 and/or angle.
    analyze       Analyse an existing rebound_data.dat without re-running.
    plot          Generate figures from existing output data.

Usage examples
--------------
    python -m python.runner equilibrate --T 0.1 --Nx 8 --Ny 8 --Nz 6
    python -m python.runner impact --v0 2.5 --angle 0.0
    python -m python.runner sweep --config sweep.json
    python -m python.runner analyze --data data/outputs/v0_2.50_ang_0.0/rebound_data.dat
    python -m python.runner plot trajectory --data data/outputs/.../rebound_data.dat
"""

from __future__ import annotations

import argparse
import concurrent.futures
import logging
import sys
from pathlib import Path

import numpy as np

from .analysis import ParametricAnalyzer, ReboundAnalyzer
from .config import (
    ConfigRegistry,
    EquilibrationConfig,
    ImpactConfig,
    ParametricSweepConfig,
)
from .exceptions import ReboundMDError
from .simulation import EquilibrationSimulation, ImpactSimulation, LAMMPSRunner, SimulationPipeline
from .visualization import plot_cor_vs_v0, plot_eac_vs_Estar, plot_trajectory

# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

_BASE_OUTPUT = Path("data/outputs")


# ==============================================================================
# Sub-command handlers
# ==============================================================================

def cmd_equilibrate(args: argparse.Namespace) -> int:
    """Run the equilibration phase."""
    cfg = EquilibrationConfig(
        T_target=args.T,
        Nx=args.Nx,
        Ny=args.Ny,
        Nz=args.Nz,
        steps_equil=args.steps_equil,
        steps_prod=args.steps_prod,
        dt=args.dt,
        seed=args.seed,
    )
    runner = LAMMPSRunner(mpi_procs=args.mpi)
    workdir = _BASE_OUTPUT / f"equilibration_T{args.T:.3f}"

    with EquilibrationSimulation(cfg, workdir, runner) as sim:
        sim.run()
        logger.info("Restart written to: %s", sim.restart_path)

    _print_config_summary("Equilibration", cfg)
    return 0


def cmd_impact(args: argparse.Namespace) -> int:
    """Run a single impact simulation."""
    equil_cfg = EquilibrationConfig(T_target=args.T)
    impact_cfg = ImpactConfig(
        v0=args.v0,
        angle=args.angle,
        T_target=args.T,
        dt=args.dt,
        steps_approach=args.steps_approach,
        steps_impact=args.steps_impact,
    )

    runner = LAMMPSRunner(mpi_procs=args.mpi)
    equil_dir = _BASE_OUTPUT / f"equilibration_T{args.T:.3f}"
    pipeline = SimulationPipeline.from_configs(
        equil_cfg, impact_cfg,
        base_workdir=_BASE_OUTPUT,
        runner=runner,
        equil_dir=equil_dir,
    )

    with pipeline:
        pipeline.run(skip_equilibration_if_exists=not args.force_equil)

    _print_config_summary("Impact", impact_cfg)

    # Immediate analysis
    analyzer = ReboundAnalyzer(v0=args.v0, angle=args.angle, T_surface_init=args.T)
    analyzer.load(pipeline.rebound_data_path)
    result = analyzer.compute_all()
    _print_result(result)

    return 0


def _run_one_job(job: tuple) -> "ReboundAnalyzer | None":
    """Module-level worker for ProcessPoolExecutor (must be picklable).

    Parameters
    ----------
    job : (ImpactConfig, replica_index, restart_path, LAMMPSRunner)
    """
    cfg, replica, restart_path, runner = job
    tag = f"T{cfg.T_target:.3f}_v0_{cfg.v0:.3f}_ang_{cfg.angle:.1f}_rep{replica}"
    workdir = _BASE_OUTPUT / tag
    sim = ImpactSimulation(cfg, workdir, restart_path, runner)
    try:
        with sim:
            sim.run()
        analyzer = ReboundAnalyzer(v0=cfg.v0, angle=cfg.angle,
                                   T_surface_init=cfg.T_target)
        analyzer.load(sim.rebound_data_path)
        return analyzer
    except ReboundMDError as exc:
        logger.error("Job %s failed: %s", tag, exc)
        return None


def cmd_sweep(args: argparse.Namespace) -> int:
    """Run a parametric sweep over v0, angle, and/or surface temperature."""
    if args.config:
        sweep_cfg = ConfigRegistry.from_json_auto(args.config)
        if not isinstance(sweep_cfg, ParametricSweepConfig):
            logger.error("Config file must be a ParametricSweepConfig, got %s", type(sweep_cfg).__name__)
            return 1
    else:
        sweep_cfg = ParametricSweepConfig(
            v0_min=args.v0_min,
            v0_max=args.v0_max,
            v0_steps=args.v0_steps,
            angle_list=",".join(str(a) for a in args.angles),
            T_surface=args.T,
            T_surface_list=args.T_list,
            n_replicas=args.replicas,
            parallel_jobs=args.jobs,
        )

    temperatures = sweep_cfg.temperatures
    runner = LAMMPSRunner(mpi_procs=args.mpi)

    # --- Phase 1: one equilibration per temperature (sequential) --------------
    # Each T* requires its own equilibrated surface.
    restart_paths: dict[float, Path] = {}
    for T in temperatures:
        equil_cfg = EquilibrationConfig(T_target=T)
        equil_workdir = _BASE_OUTPUT / f"equilibration_T{T:.3f}"
        logger.info("Equilibrating at T* = %.3f ...", T)
        with EquilibrationSimulation(equil_cfg, equil_workdir, runner) as sim:
            sim.run()
            restart_paths[T] = sim.restart_path

    # --- Phase 2: build all impact jobs (all T × v0 × angle × replica) -------
    Job = tuple[ImpactConfig, int, Path]  # (config, replica_index, restart_path)
    jobs: list[Job] = []
    for T in temperatures:
        for v0 in sweep_cfg.v0_values:
            for angle in sweep_cfg.angles:
                for replica in range(sweep_cfg.n_replicas):
                    cfg = ImpactConfig(
                        v0=v0,
                        angle=angle,
                        T_target=T,
                    )
                    jobs.append((cfg, replica, restart_paths[T]))

    n_T = len(temperatures)
    logger.info(
        "Parametric sweep: %d T* × %d v0 × %d angles × %d replicas = %d simulations",
        n_T, sweep_cfg.v0_steps, len(sweep_cfg.angles), sweep_cfg.n_replicas, len(jobs),
    )

    # --- Phase 3: run impacts in parallel -------------------------------------
    # runner is included in each job tuple so _run_one_job (module-level,
    # picklable) can access it from a worker process.
    parametric = ParametricAnalyzer()
    jobs_with_runner = [(cfg, replica, restart_path, runner) for cfg, replica, restart_path in jobs]

    with concurrent.futures.ProcessPoolExecutor(
        max_workers=sweep_cfg.parallel_jobs
    ) as pool:
        for analyzer in pool.map(_run_one_job, jobs_with_runner):
            if analyzer is not None:
                parametric.add(analyzer.compute_all())

    logger.info("Sweep complete: %d successful runs", len(parametric))
    print(parametric.summary())

    if len(parametric) < 2:
        return 0

    # --- COR vs v0 plot -------------------------------------------------------
    v0_u, cor_m, cor_s = parametric.cor_vs_v0()
    try:
        A_cor, beta_cor = parametric.fit_cor_power_law()
        logger.info("COR power-law fit: e = %.4f * v0^(-%.4f)", A_cor, beta_cor)
        cor_fit = (A_cor, beta_cor)
    except Exception:
        cor_fit = None
    plot_cor_vs_v0(
        v0_u, cor_m, cor_s,
        fit_params=cor_fit,
        output_path=_BASE_OUTPUT / "cor_vs_v0.pdf",
    )

    # --- α_E vs E* plot (Gonzalo's general observable) -----------------------
    if len(parametric) >= 4:
        E_u, eac_m, eac_s = parametric.eac_vs_Estar()
        try:
            A_eac, beta_eac = parametric.fit_eac_power_law()
            logger.info("EAC power-law fit: α_E = %.4f * E*^(%.4f)", A_eac, beta_eac)
            eac_fit = (A_eac, beta_eac)
        except Exception:
            eac_fit = None

        # One curve per T* if multiple temperatures were swept
        unique_T = np.unique(parametric.T_surface_array)
        if len(unique_T) > 1:
            # Plot α_E vs E* with one series per T*
            plot_eac_vs_Estar(
                E_u, eac_m, eac_s,
                fit_params=eac_fit,
                output_path=_BASE_OUTPUT / "eac_vs_Estar.pdf",
                T_label="all T*",
            )
            # Also plot α_E vs T* to show thermal dependence
            T_u, eac_T_m, eac_T_s = parametric.eac_vs_T()
            plot_eac_vs_Estar(
                T_u, eac_T_m, eac_T_s,
                output_path=_BASE_OUTPUT / "eac_vs_T.pdf",
                xlabel=r"Surface temperature  T*  [ε/k$_B$]",
                title=r"Energy Accommodation Coefficient vs. T*",
            )
        else:
            plot_eac_vs_Estar(
                E_u, eac_m, eac_s,
                fit_params=eac_fit,
                output_path=_BASE_OUTPUT / "eac_vs_Estar.pdf",
            )

    return 0


def cmd_analyze(args: argparse.Namespace) -> int:
    """Analyse an existing rebound_data.dat."""
    data_path = Path(args.data)
    analyzer = ReboundAnalyzer(
        v0=args.v0,
        angle=args.angle,
        T_surface_init=args.T,
        surface_top_z=args.surface_top_z,
    )
    analyzer.load(data_path)
    result = analyzer.compute_all()
    _print_result(result)
    return 0


def cmd_plot(args: argparse.Namespace) -> int:
    """Generate publication figures from existing output."""
    from .analysis import ReboundDataParser

    data_path = Path(args.data)
    data = ReboundDataParser.parse(data_path)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.type in ("trajectory", "all"):
        plot_trajectory(data, output_path=out_dir / "trajectory.pdf")
        logger.info("Saved trajectory.pdf")

    return 0


# ==============================================================================
# Helpers
# ==============================================================================

def _print_config_summary(phase: str, cfg: object) -> None:
    print(f"\n{'='*50}")
    print(f"  {phase} Configuration")
    print(f"{'='*50}")
    for k, v in cfg.to_dict().items():  # type: ignore[attr-defined]
        print(f"  {k:25s} = {v}")
    print()


def _print_result(result: object) -> None:
    from .analysis import ReboundResult
    assert isinstance(result, ReboundResult)
    print(f"\n{'='*54}")
    print(f"  Rebound Analysis Results")
    print(f"{'='*54}")
    print(f"  Incident velocity   v0_z    : {result.v0_z:.4f}  σ/τ")
    print(f"  Rebound velocity    vf_z    : {result.vf_z:.4f}  σ/τ")
    print(f"  COR                  e      : {result.COR:.5f}")
    print(f"  Reduced energy       E*     : {result.E_star:.5f}  ε")
    print(f"{'─'*54}")
    print(f"  EAC (total)          α_E    : {result.EAC:.5f}")
    print(f"  EAC (normal)         α_n    : {result.EAC_normal:.5f}")
    print(f"  EAC (tangential)     α_t    : {result.EAC_tangential:.5f}")
    print(f"{'─'*54}")
    print(f"  Energy dissipated    ΔE     : {result.delta_E:.5f}  ε")
    print(f"  Contact time         τ_c    : {result.contact_time:.4f}  τ")
    print(f"  Max penetration      δ      : {result.max_penetration:.4f}  σ")
    print(f"  Surface ΔT                  : {result.dT_surface:.5f}  ε/k_B")
    print()


# ==============================================================================
# Argument parser
# ==============================================================================

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m python.runner",
        description="Inelastic Rebound MD – simulation runner and analyser",
    )
    parser.add_argument("--mpi", type=int, default=1,
                        help="Number of MPI ranks for LAMMPS (default: 1)")

    sub = parser.add_subparsers(dest="command", required=True)

    # ---------- equilibrate ----------
    p_equil = sub.add_parser("equilibrate", help="Run surface equilibration")
    p_equil.add_argument("--T",           type=float, default=0.1)
    p_equil.add_argument("--Nx",          type=int,   default=8)
    p_equil.add_argument("--Ny",          type=int,   default=8)
    p_equil.add_argument("--Nz",          type=int,   default=6)
    p_equil.add_argument("--steps-equil", type=int,   default=50000, dest="steps_equil")
    p_equil.add_argument("--steps-prod",  type=int,   default=30000, dest="steps_prod")
    p_equil.add_argument("--dt",          type=float, default=0.005)
    p_equil.add_argument("--seed",        type=int,   default=12345)

    # ---------- impact ----------
    p_impact = sub.add_parser("impact", help="Run single impact simulation")
    p_impact.add_argument("--v0",             type=float, default=2.0)
    p_impact.add_argument("--angle",          type=float, default=0.0)
    p_impact.add_argument("--T",              type=float, default=0.1)
    p_impact.add_argument("--dt",             type=float, default=0.002)
    p_impact.add_argument("--steps-approach", type=int,   default=10000, dest="steps_approach")
    p_impact.add_argument("--steps-impact",   type=int,   default=80000, dest="steps_impact")
    p_impact.add_argument("--force-equil",    action="store_true",  dest="force_equil",
                          help="Re-run equilibration even if restart exists")

    # ---------- sweep ----------
    p_sweep = sub.add_parser("sweep", help="Parametric sweep over v0/angle/T")
    p_sweep.add_argument("--config",   type=str, default=None,
                         help="JSON config file (ParametricSweepConfig)")
    p_sweep.add_argument("--v0-min",   type=float, default=0.5,  dest="v0_min")
    p_sweep.add_argument("--v0-max",   type=float, default=5.0,  dest="v0_max")
    p_sweep.add_argument("--v0-steps", type=int,   default=10,   dest="v0_steps")
    p_sweep.add_argument("--angles",   type=float, nargs="+",    default=[0.0])
    p_sweep.add_argument("--T",        type=float, default=0.1,
                         help="Single surface temperature [eps/kB] (default: 0.1)")
    p_sweep.add_argument("--T-list",   type=str,   default="",  dest="T_list",
                         help="Comma-separated T* values, e.g. '0.05,0.1,0.2,0.3'. "
                              "Overrides --T; runs one equilibration per temperature.")
    p_sweep.add_argument("--replicas", type=int,   default=3)
    p_sweep.add_argument("--jobs",     type=int,   default=4)

    # ---------- analyze ----------
    p_an = sub.add_parser("analyze", help="Analyse existing rebound_data.dat")
    p_an.add_argument("--data",            required=True)
    p_an.add_argument("--v0",             type=float, required=True)
    p_an.add_argument("--angle",          type=float, default=0.0)
    p_an.add_argument("--T",              type=float, default=0.1)
    p_an.add_argument("--surface-top-z",  type=float, default=9.3, dest="surface_top_z")

    # ---------- plot ----------
    p_plot = sub.add_parser("plot", help="Generate figures from output data")
    p_plot.add_argument("type", choices=["trajectory", "cor", "energy", "all"])
    p_plot.add_argument("--data",       required=True)
    p_plot.add_argument("--output-dir", default="data/outputs/figures", dest="output_dir")

    return parser


# ==============================================================================
# Entry point
# ==============================================================================

def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    dispatch = {
        "equilibrate": cmd_equilibrate,
        "impact":      cmd_impact,
        "sweep":       cmd_sweep,
        "analyze":     cmd_analyze,
        "plot":        cmd_plot,
    }

    try:
        return dispatch[args.command](args)
    except ReboundMDError as exc:
        logger.error("Simulation error: %s", exc)
        return 1
    except KeyboardInterrupt:
        logger.info("Interrupted by user.")
        return 130
    except Exception as exc:
        logger.exception("Unexpected error: %s", exc)
        return 2


if __name__ == "__main__":
    sys.exit(main())
