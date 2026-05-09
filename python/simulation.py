"""
simulation.py
=============
OOP wrapper around LAMMPS execution for the inelastic-rebound-md project.

Class hierarchy:
    LAMMPSRunner           – low-level subprocess launcher
    EquilibrationSimulation – Phase 1 wrapper
    ImpactSimulation        – Phase 2 wrapper
    SimulationPipeline      – orchestrates both phases end-to-end

Context manager support:
    Each simulation class is a context manager that manages a working directory
    and cleans up temporary files on exit.

Usage
-----
>>> from python.config import ImpactConfig
>>> from python.simulation import SimulationPipeline
>>> with SimulationPipeline.from_configs(equil_cfg, impact_cfg) as pipeline:
...     pipeline.run()
...     cor = pipeline.coefficient_of_restitution
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Optional

from .config import ClusterImpactConfig, EquilibrationConfig, ImpactConfig
from .exceptions import (
    LAMMPSNotFoundError,
    RestartNotFoundError,
    SimulationRuntimeError,
    TimeoutError as SimTimeoutError,
)

logger = logging.getLogger(__name__)

# Default LAMMPS executable: override with env var LAMMPS_EXE
_DEFAULT_LAMMPS = os.environ.get("LAMMPS_EXE", "lmp")

# Location of LAMMPS input scripts relative to project root
_LAMMPS_DIR = Path(__file__).parent.parent / "lammps"


# ==============================================================================
# Low-level subprocess launcher
# ==============================================================================

class LAMMPSRunner:
    """Manages a single LAMMPS subprocess invocation.

    Parameters
    ----------
    executable : str
        Path or name of the LAMMPS binary (default: ``lmp`` or ``$LAMMPS_EXE``).
    mpi_procs : int
        Number of MPI ranks (1 = serial).
    timeout : float | None
        Wall-time limit in seconds.  ``None`` = no limit.
    """

    def __init__(
        self,
        executable: str = _DEFAULT_LAMMPS,
        mpi_procs: int = 1,
        timeout: Optional[float] = None,
    ) -> None:
        self.executable = executable
        self.mpi_procs = mpi_procs
        self.timeout = timeout
        self._verify_executable()

    def _verify_executable(self) -> None:
        """Raise LAMMPSNotFoundError if the executable is not on PATH."""
        if shutil.which(self.executable) is None:
            raise LAMMPSNotFoundError(self.executable)

    def _build_command(self, script: Path, extra_vars: dict[str, str]) -> list[str]:
        """Build the subprocess command list."""
        cmd: list[str] = []
        if self.mpi_procs > 1:
            cmd += ["mpirun", "-n", str(self.mpi_procs)]
        # Resolve to absolute path: LAMMPS runs with cwd=workdir, so relative
        # paths would be mis-interpreted if script lives outside the workdir.
        cmd += [self.executable, "-in", str(script.resolve())]
        for key, val in extra_vars.items():
            cmd += ["-var", key, val]
        return cmd

    def run(
        self,
        script: Path,
        workdir: Path,
        extra_vars: Optional[dict[str, str]] = None,
    ) -> subprocess.CompletedProcess:
        """Run LAMMPS synchronously and return the completed process.

        Parameters
        ----------
        script : Path
            Absolute path to the LAMMPS input script.
        workdir : Path
            Working directory; output files are written here.
        extra_vars : dict, optional
            Additional ``-var key value`` command-line overrides.

        Raises
        ------
        SimulationRuntimeError
            If LAMMPS exits with non-zero return code.
        SimTimeoutError
            If wall-time limit is exceeded.
        """
        vars_ = extra_vars or {}
        cmd = self._build_command(script, vars_)

        logger.info("Launching LAMMPS: %s", " ".join(cmd))
        logger.info("Working directory: %s", workdir)

        # When lmp is MPI-linked, concurrent subprocess invocations conflict on
        # shared PMIx/ORTE resources.  Point each job at its own tmpdir so the
        # MPI daemon doesn't collide with sibling jobs.
        env = os.environ.copy()
        env["OMPI_MCA_orte_tmpdir_base"] = str(workdir)
        env["PMIX_MCA_psec"] = "native"

        t0 = time.perf_counter()
        try:
            result = subprocess.run(
                cmd,
                cwd=workdir,
                capture_output=True,
                text=True,
                timeout=self.timeout,
                env=env,
            )
        except subprocess.TimeoutExpired as exc:
            elapsed = time.perf_counter() - t0
            raise SimTimeoutError(elapsed, self.timeout, str(script)) from exc

        elapsed = time.perf_counter() - t0
        logger.info("LAMMPS finished in %.2f s (returncode=%d)", elapsed, result.returncode)

        if result.returncode != 0:
            # With Open MPI, MPI_Abort flushes log.lammps to disk but may not
            # flush stdout to the pipe.  Read the log file directly for the
            # real LAMMPS error message.
            log_file = workdir / "log.lammps"
            if log_file.exists():
                combined = log_file.read_text()
            else:
                combined = result.stdout + "\n" + result.stderr
            raise SimulationRuntimeError(result.returncode, combined, str(script))

        return result


# ==============================================================================
# Phase 1: Equilibration
# ==============================================================================

class EquilibrationSimulation:
    """Wraps the LAMMPS equilibration phase (in.equilibrate).

    Parameters
    ----------
    config : EquilibrationConfig
        Physical/numerical parameters.
    workdir : Path | str
        Directory where output files will be written.
    runner : LAMMPSRunner | None
        Pre-configured runner; created with defaults if ``None``.
    """

    SCRIPT_NAME = "in.equilibrate"
    RESTART_NAME = "restart.equilibrated"

    def __init__(
        self,
        config: EquilibrationConfig,
        workdir: Path | str,
        runner: Optional[LAMMPSRunner] = None,
    ) -> None:
        self.config = config
        self.workdir = Path(workdir)
        self.runner = runner or LAMMPSRunner()
        self._completed = False

    # --- Context manager -------------------------------------------------------

    def __enter__(self) -> "EquilibrationSimulation":
        self.workdir.mkdir(parents=True, exist_ok=True)
        logger.debug("EquilibrationSimulation workdir: %s", self.workdir)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        if exc_type is not None:
            logger.error("EquilibrationSimulation exited with error: %s", exc_val)
        return False   # do not suppress exceptions

    # --- Public API ------------------------------------------------------------

    @property
    def restart_path(self) -> Path:
        return self.workdir / self.RESTART_NAME

    @property
    def is_complete(self) -> bool:
        return self._completed and self.restart_path.exists()

    def run(self) -> None:
        """Execute the equilibration simulation.

        The method copies the LAMMPS script into ``workdir``, then runs LAMMPS,
        overriding script variables with values from ``self.config``.
        """
        script_src = _LAMMPS_DIR / self.SCRIPT_NAME
        script_dst = self.workdir / self.SCRIPT_NAME
        shutil.copy2(script_src, script_dst)

        extra_vars = self._config_to_lammps_vars()
        self.runner.run(script_dst, self.workdir, extra_vars)

        if not self.restart_path.exists():
            raise FileNotFoundError(
                f"LAMMPS did not produce expected restart file: {self.restart_path}"
            )
        self._completed = True
        logger.info("Equilibration complete. Restart: %s", self.restart_path)

    def _config_to_lammps_vars(self) -> dict[str, str]:
        """Convert config fields to LAMMPS -var overrides."""
        cfg = self.config
        return {
            "T_target":    str(cfg.T_target),
            "T_damp":      str(cfg.T_damp),
            "a0":          str(cfg.a0),
            "Nx":          str(cfg.Nx),
            "Ny":          str(cfg.Ny),
            "Nz":          str(cfg.Nz),
            "steps_equil": str(cfg.steps_equil),
            "steps_prod":  str(cfg.steps_prod),
            "dt":          str(cfg.dt),
            "seed":        str(cfg.seed),
        }

    def __repr__(self) -> str:
        return (
            f"EquilibrationSimulation(config={self.config!r}, "
            f"workdir={self.workdir!r}, completed={self._completed})"
        )


# ==============================================================================
# Phase 2: Impact
# ==============================================================================

class ImpactSimulation:
    """Wraps the LAMMPS impact simulation phase (in.impact).

    Requires that equilibration has been run first (reads restart file).

    Parameters
    ----------
    config : ImpactConfig
        Physical/numerical parameters for the impact.
    workdir : Path | str
        Directory where output files will be written.
    restart_path : Path | str
        Path to the restart file from equilibration.
    runner : LAMMPSRunner | None
        Pre-configured runner; created with defaults if ``None``.
    """

    SCRIPT_NAME = "in.impact"
    REBOUND_DATA = "rebound_data.dat"
    DUMP_FILE = "dump.impact.lammpstrj"

    def __init__(
        self,
        config: ImpactConfig,
        workdir: Path | str,
        restart_path: Path | str,
        runner: Optional[LAMMPSRunner] = None,
    ) -> None:
        self.config = config
        self.workdir = Path(workdir)
        self.restart_path = Path(restart_path)
        self.runner = runner or LAMMPSRunner()
        self._completed = False
        self._raw_result: Optional[subprocess.CompletedProcess] = None

    # --- Context manager -------------------------------------------------------

    def __enter__(self) -> "ImpactSimulation":
        self.workdir.mkdir(parents=True, exist_ok=True)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        return False

    # --- Public API ------------------------------------------------------------

    @property
    def rebound_data_path(self) -> Path:
        return self.workdir / self.REBOUND_DATA

    @property
    def dump_path(self) -> Path:
        return self.workdir / self.DUMP_FILE

    @property
    def is_complete(self) -> bool:
        return self._completed and self.rebound_data_path.exists()

    def run(self, skip_if_exists: bool = False) -> None:
        """Execute the impact simulation.

        Parameters
        ----------
        skip_if_exists : bool
            If True and ``rebound_data.dat`` already exists, skip the LAMMPS
            run.  Useful for re-running analysis without re-simulating.
        """
        if skip_if_exists and self.rebound_data_path.exists():
            logger.info("rebound_data.dat found; skipping impact simulation.")
            self._completed = True
            return

        if not self.restart_path.exists():
            raise RestartNotFoundError(self.restart_path)

        script_src = _LAMMPS_DIR / self.SCRIPT_NAME
        script_dst = self.workdir / self.SCRIPT_NAME
        shutil.copy2(script_src, script_dst)

        # Symlink restart into workdir so LAMMPS can find it by relative path
        restart_link = self.workdir / "restart.equilibrated"
        if not restart_link.exists():
            restart_link.symlink_to(self.restart_path.resolve())

        extra_vars = self._config_to_lammps_vars()
        self._raw_result = self.runner.run(script_dst, self.workdir, extra_vars)
        self._completed = True
        logger.info("Impact simulation complete. Data: %s", self.rebound_data_path)

    def _config_to_lammps_vars(self) -> dict[str, str]:
        cfg = self.config
        return {
            "v0":              str(cfg.v0),
            "angle":           str(cfg.angle),
            "T_target":        str(cfg.T_target),
            "T_damp":          str(cfg.T_damp),
            "dt":              str(cfg.dt),
            "z_gap":           str(cfg.z_gap),
            "steps_approach":  str(cfg.steps_approach),
            "steps_impact":    str(cfg.steps_impact),
            "thermo_freq":     str(cfg.thermo_freq),
            "dump_freq":       str(cfg.dump_freq),
            "eps_ps":          str(cfg.eps_ps),
        }

    def __repr__(self) -> str:
        return (
            f"ImpactSimulation(v0={self.config.v0}, angle={self.config.angle}, "
            f"T={self.config.T_target}, completed={self._completed})"
        )


# ==============================================================================
# Pipeline: orchestrates both phases
# ==============================================================================

class SimulationPipeline:
    """End-to-end pipeline: equilibration -> impact -> analysis-ready output.

    Parameters
    ----------
    equil_sim : EquilibrationSimulation
    impact_sim : ImpactSimulation
    """

    def __init__(
        self,
        equil_sim: EquilibrationSimulation,
        impact_sim: ImpactSimulation,
    ) -> None:
        self.equil_sim = equil_sim
        self.impact_sim = impact_sim

    @classmethod
    def from_configs(
        cls,
        equil_config: EquilibrationConfig,
        impact_config: ImpactConfig,
        base_workdir: Path | str = Path("data/outputs"),
        runner: Optional[LAMMPSRunner] = None,
        run_id: Optional[str] = None,
        equil_dir: Optional[Path | str] = None,
    ) -> "SimulationPipeline":
        """Construct a pipeline from config objects.

        Automatically creates unique subdirectories under ``base_workdir``.
        ``equil_dir`` overrides the default equilibration directory, which is
        useful when running at multiple temperatures.
        """
        base = Path(base_workdir)
        tag = run_id or f"v0_{impact_config.v0:.2f}_ang_{impact_config.angle:.1f}"

        if equil_dir is not None:
            equil_dir = Path(equil_dir)
        else:
            equil_dir = base / f"equilibration_T{equil_config.T_target:.3f}"
        impact_dir = base / tag

        runner_ = runner or LAMMPSRunner()

        equil_sim = EquilibrationSimulation(equil_config, equil_dir, runner_)
        impact_sim = ImpactSimulation(
            impact_config,
            impact_dir,
            restart_path=equil_dir / EquilibrationSimulation.RESTART_NAME,
            runner=runner_,
        )
        return cls(equil_sim, impact_sim)

    # --- Context manager -------------------------------------------------------

    def __enter__(self) -> "SimulationPipeline":
        self.equil_sim.__enter__()
        self.impact_sim.__enter__()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        self.equil_sim.__exit__(exc_type, exc_val, exc_tb)
        self.impact_sim.__exit__(exc_type, exc_val, exc_tb)
        return False

    # --- Public API ------------------------------------------------------------

    def run(self, skip_equilibration_if_exists: bool = True) -> None:
        """Run the full pipeline.

        Parameters
        ----------
        skip_equilibration_if_exists : bool
            If ``True`` and a restart file already exists, skip equilibration.
            Useful for parametric sweeps that reuse the same surface.
        """
        if skip_equilibration_if_exists and self.equil_sim.restart_path.exists():
            logger.info("Restart file found; skipping equilibration.")
        else:
            logger.info("Phase 1: Running equilibration...")
            self.equil_sim.run()

        logger.info("Phase 2: Running impact simulation...")
        self.impact_sim.run()
        logger.info("Pipeline complete.")

    @property
    def rebound_data_path(self) -> Path:
        """Path to the per-timestep projectile tracking file."""
        return self.impact_sim.rebound_data_path

    def __repr__(self) -> str:
        return (
            f"SimulationPipeline(\n"
            f"  equil={self.equil_sim!r},\n"
            f"  impact={self.impact_sim!r}\n)"
        )


# ==============================================================================
# Cluster minimisation
# ==============================================================================

class ClusterMinimiser:
    """Minimise an N-atom LJ cluster using LAMMPS in vacuum.

    Reads ``clusters/cluster_N{N}_init.dat`` (written by ``cluster_gen.py``)
    and produces ``clusters/cluster_N{N}_min.dat`` (minimised coordinates).

    Parameters
    ----------
    N : int
        Cluster size.
    workdir : Path | str
        Directory where the minimisation is run (LAMMPS writes log there).
    runner : LAMMPSRunner | None
        Pre-configured runner.
    """

    SCRIPT_NAME = "in.minimize_cluster"

    def __init__(
        self,
        N: int,
        workdir: Path | str,
        runner: Optional[LAMMPSRunner] = None,
    ) -> None:
        self.N = N
        self.workdir = Path(workdir)
        self.runner = runner or LAMMPSRunner()

    @property
    def init_dat(self) -> Path:
        return _LAMMPS_DIR.parent / "lammps" / "clusters" / f"cluster_N{self.N}_init.dat"

    @property
    def min_dat(self) -> Path:
        return _LAMMPS_DIR.parent / "lammps" / "clusters" / f"cluster_N{self.N}_min.dat"

    def run(self, skip_if_exists: bool = True) -> Path:
        """Run the minimisation and return the path to the minimised data file.

        Parameters
        ----------
        skip_if_exists : bool
            If True and the minimised file already exists, skip the LAMMPS run.
        """
        if skip_if_exists and self.min_dat.exists():
            logger.info("Minimised cluster N=%d found; skipping.", self.N)
            return self.min_dat

        if not self.init_dat.exists():
            raise FileNotFoundError(
                f"Initial cluster data not found: {self.init_dat}\n"
                "Run: python -m python.runner cluster-setup"
            )

        self.workdir.mkdir(parents=True, exist_ok=True)

        script = _LAMMPS_DIR / self.SCRIPT_NAME
        extra_vars = {
            "N_cluster":  str(self.N),
            "init_file":  str(self.init_dat.resolve()),
            "out_file":   str(self.min_dat.resolve()),
        }
        self.runner.run(script, self.workdir, extra_vars)

        if not self.min_dat.exists():
            raise FileNotFoundError(
                f"LAMMPS minimisation did not produce: {self.min_dat}"
            )
        logger.info("Cluster N=%d minimised → %s", self.N, self.min_dat)
        return self.min_dat


# ==============================================================================
# Phase 2 (cluster variant): Cluster Impact Simulation
# ==============================================================================

class ClusterImpactSimulation:
    """Wraps the LAMMPS cluster impact simulation (in.impact_cluster).

    Tracks both the cluster COM motion and its internal kinetic energy,
    enabling computation of α_E_int (energy to internal modes) and
    α_E_surf (energy to surface phonons) — the central observables for
    the cluster size dependence study.

    Parameters
    ----------
    config : ClusterImpactConfig
        Physical/numerical parameters (includes N_cluster, cluster_file).
    workdir : Path | str
        Directory where output files are written.
    restart_path : Path | str
        Path to the restart file from equilibration.
    runner : LAMMPSRunner | None
        Pre-configured runner.
    """

    SCRIPT_NAME = "in.impact_cluster"
    REBOUND_DATA = "cluster_rebound_data.dat"
    DUMP_FILE = "dump.cluster_impact.lammpstrj"

    def __init__(
        self,
        config: ClusterImpactConfig,
        workdir: Path | str,
        restart_path: Path | str,
        runner: Optional[LAMMPSRunner] = None,
    ) -> None:
        self.config = config
        self.workdir = Path(workdir)
        self.restart_path = Path(restart_path)
        self.runner = runner or LAMMPSRunner()
        self._completed = False

    # --- Context manager -------------------------------------------------------

    def __enter__(self) -> "ClusterImpactSimulation":
        self.workdir.mkdir(parents=True, exist_ok=True)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        return False

    # --- Public API ------------------------------------------------------------

    @property
    def rebound_data_path(self) -> Path:
        return self.workdir / self.REBOUND_DATA

    @property
    def dump_path(self) -> Path:
        return self.workdir / self.DUMP_FILE

    @property
    def is_complete(self) -> bool:
        return self._completed and self.rebound_data_path.exists()

    def run(self, skip_if_exists: bool = False) -> None:
        """Execute the cluster impact simulation.

        Parameters
        ----------
        skip_if_exists : bool
            If True and ``cluster_rebound_data.dat`` already exists, skip.
        """
        if skip_if_exists and self.rebound_data_path.exists():
            logger.info("cluster_rebound_data.dat found; skipping.")
            self._completed = True
            return

        if not self.restart_path.exists():
            raise RestartNotFoundError(self.restart_path)

        script_src = _LAMMPS_DIR / self.SCRIPT_NAME
        script_dst = self.workdir / self.SCRIPT_NAME
        shutil.copy2(script_src, script_dst)

        # Symlink restart into workdir
        restart_link = self.workdir / "restart.equilibrated"
        if not restart_link.exists():
            restart_link.symlink_to(self.restart_path.resolve())

        extra_vars = self._config_to_lammps_vars()
        self.runner.run(script_dst, self.workdir, extra_vars)
        self._completed = True
        logger.info("Cluster impact complete. Data: %s", self.rebound_data_path)

    def _config_to_lammps_vars(self) -> dict[str, str]:
        cfg = self.config
        return {
            "v0":              str(cfg.v0),
            "angle":           str(cfg.angle),
            "T_target":        str(cfg.T_target),
            "T_damp":          str(cfg.T_damp),
            "dt":              str(cfg.dt),
            "z_gap":           str(cfg.z_gap),
            "steps_approach":  str(cfg.steps_approach),
            "steps_impact":    str(cfg.steps_impact),
            "thermo_freq":     str(cfg.thermo_freq),
            "dump_freq":       str(cfg.dump_freq),
            "eps_ps":          str(cfg.eps_ps),
            "N_cluster":       str(cfg.N_cluster),
            "cluster_file":    str(Path(cfg.cluster_file).resolve()),
            "seed":            str(cfg.seed),
        }

    def __repr__(self) -> str:
        return (
            f"ClusterImpactSimulation(N={self.config.N_cluster}, "
            f"v0={self.config.v0}, T={self.config.T_target}, "
            f"completed={self._completed})"
        )
