"""
sweep.py
========
Parametric sweep engine for the inelastic-rebound-md project.

Provides utilities for designing, executing and aggregating batches of
``ImpactSimulation`` runs that explore the (v0, angle, T) parameter space.

Usage
-----
>>> from python.config import EquilibrationConfig, ParametricSweepConfig
>>> from python.sweep import ParametricSweep
>>> sweep = ParametricSweep.from_config(sweep_cfg, equil_cfg)
>>> sweep.run(jobs=4)
>>> sweep.analyzer.summary()
"""

from __future__ import annotations

import concurrent.futures
import itertools
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from .analysis import ParametricAnalyzer, ReboundAnalyzer, ReboundResult
from .config import EquilibrationConfig, ImpactConfig, ParametricSweepConfig
from .exceptions import ReboundMDError
from .simulation import (
    EquilibrationSimulation,
    ImpactSimulation,
    LAMMPSRunner,
)

logger = logging.getLogger(__name__)

_BASE_OUTPUT = Path("data/outputs")


# ==============================================================================
# Job descriptor
# ==============================================================================

@dataclass(frozen=True)
class SweepJob:
    """Immutable description of a single simulation job in the sweep.

    Attributes
    ----------
    v0     : float   – incident speed [sigma/tau]
    angle  : float   – impact angle [degrees]
    T      : float   – surface temperature [eps/kB]
    replica: int     – replica index (different seeds / initial velocities)
    run_id : str     – unique identifier used as output subdirectory name
    """

    v0:      float
    angle:   float
    T:       float
    replica: int
    run_id:  str

    def to_impact_config(self) -> ImpactConfig:
        return ImpactConfig(v0=self.v0, angle=self.angle, T_target=self.T)


# ==============================================================================
# Main sweep class
# ==============================================================================

class ParametricSweep:
    """Plans and executes a parametric sweep, collects results.

    Parameters
    ----------
    jobs : list[SweepJob]
        All jobs to run.
    equil_config : EquilibrationConfig
        Configuration for the shared surface equilibration.
    base_workdir : Path
        Root output directory.
    runner : LAMMPSRunner | None
        Shared LAMMPS runner.
    """

    def __init__(
        self,
        jobs: list[SweepJob],
        equil_config: EquilibrationConfig,
        base_workdir: Path = _BASE_OUTPUT,
        runner: Optional[LAMMPSRunner] = None,
    ) -> None:
        self.jobs = jobs
        self.equil_config = equil_config
        self.base_workdir = Path(base_workdir)
        self.runner = runner or LAMMPSRunner()
        self.analyzer = ParametricAnalyzer()
        self._restart_path: Optional[Path] = None

    @classmethod
    def from_config(
        cls,
        sweep_cfg: ParametricSweepConfig,
        equil_cfg: EquilibrationConfig,
        base_workdir: Path = _BASE_OUTPUT,
    ) -> "ParametricSweep":
        """Build a ``ParametricSweep`` from ``ParametricSweepConfig``."""
        jobs = cls._build_jobs(sweep_cfg)
        return cls(jobs, equil_cfg, base_workdir)

    @staticmethod
    def _build_jobs(cfg: ParametricSweepConfig) -> list[SweepJob]:
        """Enumerate all (v0, angle, replica) combinations."""
        jobs = []
        for v0, angle, replica in itertools.product(
            cfg.v0_values, cfg.angles, range(cfg.n_replicas)
        ):
            run_id = f"v0_{v0:.3f}_ang_{angle:.1f}_T_{cfg.T_surface:.3f}_rep{replica}"
            jobs.append(SweepJob(v0=v0, angle=angle, T=cfg.T_surface,
                                  replica=replica, run_id=run_id))
        return jobs

    # --- Public API -----------------------------------------------------------

    def run(
        self,
        jobs: Optional[int] = 4,
        skip_equilibration_if_exists: bool = True,
    ) -> ParametricAnalyzer:
        """Execute the full sweep and return a populated ``ParametricAnalyzer``.

        Steps:
        1. Run equilibration once (or reuse existing restart).
        2. Dispatch all impact jobs in parallel via ``ProcessPoolExecutor``.
        3. Collect results into ``self.analyzer``.

        Parameters
        ----------
        jobs : int
            Maximum concurrent LAMMPS processes.
        skip_equilibration_if_exists : bool
            If True and restart file already exists, skip equilibration.
        """
        self._ensure_equilibrated(skip_equilibration_if_exists)

        logger.info("Starting parametric sweep: %d jobs, %d workers", len(self.jobs), jobs)

        with concurrent.futures.ProcessPoolExecutor(max_workers=jobs) as pool:
            futures = {
                pool.submit(self._run_one_job, job): job
                for job in self.jobs
            }
            for future in concurrent.futures.as_completed(futures):
                job = futures[future]
                try:
                    result = future.result()
                    if result is not None:
                        self.analyzer.add(result)
                        logger.info(
                            "Job complete: %s  COR=%.4f",
                            job.run_id, result.COR,
                        )
                except Exception as exc:
                    logger.error("Job %s raised: %s", job.run_id, exc)

        logger.info(
            "Sweep finished. %d / %d jobs succeeded.",
            len(self.analyzer), len(self.jobs),
        )
        return self.analyzer

    def _ensure_equilibrated(self, skip_if_exists: bool) -> None:
        equil_dir = self.base_workdir / "equilibration"
        restart = equil_dir / EquilibrationSimulation.RESTART_NAME
        if skip_if_exists and restart.exists():
            logger.info("Reusing existing restart: %s", restart)
            self._restart_path = restart
            return

        with EquilibrationSimulation(self.equil_config, equil_dir, self.runner) as sim:
            sim.run()
            self._restart_path = sim.restart_path

    def _run_one_job(self, job: SweepJob) -> Optional[ReboundResult]:
        """Execute a single impact job. Called in a subprocess."""
        workdir = self.base_workdir / job.run_id
        config  = job.to_impact_config()

        sim = ImpactSimulation(
            config,
            workdir,
            restart_path=self._restart_path,
            runner=self.runner,
        )
        try:
            with sim:
                sim.run()
            analyzer = ReboundAnalyzer(
                v0=job.v0,
                angle=job.angle,
                T_surface_init=job.T,
            )
            analyzer.load(sim.rebound_data_path)
            return analyzer.compute_all()
        except ReboundMDError as exc:
            logger.error("SweepJob %s failed: %s", job.run_id, exc)
            return None

    # --- Summary --------------------------------------------------------------

    def design_summary(self) -> str:
        """Human-readable description of what the sweep will run."""
        v0_vals    = sorted({j.v0 for j in self.jobs})
        ang_vals   = sorted({j.angle for j in self.jobs})
        T_vals     = sorted({j.T for j in self.jobs})
        n_replicas = max({j.replica for j in self.jobs}) + 1
        lines = [
            f"ParametricSweep design",
            f"  v0     : {v0_vals}  ({len(v0_vals)} values)",
            f"  angles : {ang_vals}  ({len(ang_vals)} values)",
            f"  T*     : {T_vals}",
            f"  replicas: {n_replicas}",
            f"  Total jobs: {len(self.jobs)}",
        ]
        return "\n".join(lines)

    def __repr__(self) -> str:
        return (
            f"ParametricSweep(n_jobs={len(self.jobs)}, "
            f"n_results={len(self.analyzer)})"
        )
