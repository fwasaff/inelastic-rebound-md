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

from .analysis import ParametricAnalyzer, ReboundAnalyzer
from .config import (
    ConfigRegistry,
    EquilibrationConfig,
    ImpactConfig,
    ParametricSweepConfig,
)
from .exceptions import ReboundMDError
from .simulation import EquilibrationSimulation, ImpactSimulation, LAMMPSRunner, SimulationPipeline
from .visualization import plot_cor_vs_v0, plot_trajectory

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
    workdir = _BASE_OUTPUT / "equilibration"

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
    pipeline = SimulationPipeline.from_configs(
        equil_cfg, impact_cfg,
        base_workdir=_BASE_OUTPUT,
        runner=runner,
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


def cmd_sweep(args: argparse.Namespace) -> int:
    """Run a parametric sweep."""
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
            n_replicas=args.replicas,
            parallel_jobs=args.jobs,
        )

    equil_cfg = EquilibrationConfig(T_target=sweep_cfg.T_surface)
    runner = LAMMPSRunner(mpi_procs=args.mpi)

    # Equilibrate once; all sweeps reuse the same surface
    equil_workdir = _BASE_OUTPUT / "equilibration"
    with EquilibrationSimulation(equil_cfg, equil_workdir, runner) as sim:
        sim.run()
        restart_path = sim.restart_path

    # Build all impact jobs
    jobs: list[tuple[ImpactConfig, int]] = []
    for v0 in sweep_cfg.v0_values:
        for angle in sweep_cfg.angles:
            for replica in range(sweep_cfg.n_replicas):
                cfg = ImpactConfig(
                    v0=v0,
                    angle=angle,
                    T_target=sweep_cfg.T_surface,
                    seed=replica * 1000,  # different initial conditions per replica
                )
                jobs.append((cfg, replica))

    logger.info("Parametric sweep: %d total simulations", len(jobs))

    # Run in parallel using a process pool
    parametric = ParametricAnalyzer()

    def run_one(job: tuple[ImpactConfig, int]) -> ReboundAnalyzer | None:
        cfg, replica = job
        tag = f"v0_{cfg.v0:.3f}_ang_{cfg.angle:.1f}_rep{replica}"
        workdir = _BASE_OUTPUT / tag
        sim = ImpactSimulation(cfg, workdir, restart_path, runner)
        try:
            with sim:
                sim.run()
            analyzer = ReboundAnalyzer(v0=cfg.v0, angle=cfg.angle,
                                        T_surface_init=sweep_cfg.T_surface)
            analyzer.load(sim.rebound_data_path)
            return analyzer
        except ReboundMDError as exc:
            logger.error("Job %s failed: %s", tag, exc)
            return None

    with concurrent.futures.ProcessPoolExecutor(
        max_workers=sweep_cfg.parallel_jobs
    ) as pool:
        for analyzer in pool.map(run_one, jobs):
            if analyzer is not None:
                parametric.add(analyzer.compute_all())

    logger.info("Sweep complete: %d successful runs", len(parametric))
    print(parametric.summary())

    if len(parametric) >= 2:
        v0_u, cor_m, cor_s = parametric.cor_vs_v0()
        try:
            A, beta = parametric.fit_cor_power_law()
            logger.info("Power-law fit: COR = %.4f * v0^(-%.4f)", A, beta)
            fit_params = (A, beta)
        except Exception:
            fit_params = None

        fig = plot_cor_vs_v0(
            v0_u, cor_m, cor_s,
            fit_params=fit_params,
            output_path=_BASE_OUTPUT / "cor_vs_v0.pdf",
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
    print(f"\n{'='*50}")
    print(f"  Rebound Analysis Results")
    print(f"{'='*50}")
    print(f"  Incident velocity v0_z : {result.v0_z:.4f}  σ/τ")
    print(f"  Rebound velocity  vf_z : {result.vf_z:.4f}  σ/τ")
    print(f"  COR                    : {result.COR:.5f}")
    print(f"  Energy dissipated  ΔE  : {result.delta_E:.5f}  ε")
    print(f"  Contact time       τ_c : {result.contact_time:.4f}  τ")
    print(f"  Max penetration    δ   : {result.max_penetration:.4f}  σ")
    print(f"  Surface ΔT             : {result.dT_surface:.5f}  ε/k_B")
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
    p_sweep.add_argument("--T",        type=float, default=0.1)
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
