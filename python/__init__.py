"""
inelastic-rebound-md Python package
=====================================
Simulation management and analysis tools for LJ inelastic rebound MD.

Submodules
----------
config       – metaclass-driven configuration dataclasses
exceptions   – custom exception hierarchy
simulation   – LAMMPS subprocess wrappers and pipeline
analysis     – physical observables, decorators, parametric tools
visualization – matplotlib plotting utilities
runner       – CLI entry point
"""

from .config import (
    EquilibrationConfig,
    ImpactConfig,
    ParametricSweepConfig,
    ConfigRegistry,
)
from .exceptions import ReboundMDError
from .simulation import SimulationPipeline, LAMMPSRunner
from .analysis import ReboundAnalyzer, ParametricAnalyzer, ReboundResult

__all__ = [
    "EquilibrationConfig",
    "ImpactConfig",
    "ParametricSweepConfig",
    "ConfigRegistry",
    "ReboundMDError",
    "SimulationPipeline",
    "LAMMPSRunner",
    "ReboundAnalyzer",
    "ParametricAnalyzer",
    "ReboundResult",
]
