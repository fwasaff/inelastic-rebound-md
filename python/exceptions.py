"""
exceptions.py
=============
Custom exception hierarchy for the inelastic-rebound-md project.

Design principles:
- All project exceptions inherit from ReboundMDError for easy catching.
- Each exception carries structured context (not just a message string).
- __str__ is overridden to produce human-readable, actionable messages.
- Exceptions are categorised by subsystem: IO, simulation, physics, analysis.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional


# ==============================================================================
# Base
# ==============================================================================

class ReboundMDError(Exception):
    """Root exception for all inelastic-rebound-md errors.

    Attributes
    ----------
    message : str
        Human-readable description of what went wrong.
    context : dict
        Arbitrary key-value pairs for structured diagnostic information.
    """

    def __init__(self, message: str, **context: Any) -> None:
        super().__init__(message)
        self.message = message
        self.context = context

    def __str__(self) -> str:
        base = self.message
        if self.context:
            details = "; ".join(f"{k}={v!r}" for k, v in self.context.items())
            return f"{base} [{details}]"
        return base

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.message!r}, **{self.context!r})"


# ==============================================================================
# I/O and file handling
# ==============================================================================

class IOError(ReboundMDError):
    """Raised when a file or directory operation fails."""


class FileNotFoundError(IOError):
    """Raised when a required input file does not exist.

    Parameters
    ----------
    path : Path | str
        The missing file path.
    role : str, optional
        Human-readable description of the file's role (e.g. "LAMMPS restart").
    """

    def __init__(self, path: Path | str, role: Optional[str] = None) -> None:
        self.path = Path(path)
        self.role = role
        label = f" ({role})" if role else ""
        super().__init__(
            f"Required file not found{label}: {self.path}",
            path=str(self.path),
            role=role,
        )


class ParseError(IOError):
    """Raised when parsing a LAMMPS dump, thermo log, or config file fails.

    Parameters
    ----------
    filename : str
        Name of the file being parsed.
    line_number : int | None
        Line number where parsing failed, if known.
    raw_line : str | None
        The offending raw text.
    """

    def __init__(
        self,
        filename: str,
        *,
        line_number: Optional[int] = None,
        raw_line: Optional[str] = None,
        detail: str = "",
    ) -> None:
        self.filename = filename
        self.line_number = line_number
        self.raw_line = raw_line
        loc = f" at line {line_number}" if line_number is not None else ""
        suffix = f": {detail}" if detail else ""
        super().__init__(
            f"Parse error in '{filename}'{loc}{suffix}",
            filename=filename,
            line_number=line_number,
            raw_line=raw_line,
        )


class RestartNotFoundError(FileNotFoundError):
    """Raised when the LAMMPS restart file expected from equilibration is absent."""

    def __init__(self, path: Path | str) -> None:
        super().__init__(path, role="LAMMPS restart (run equilibration first)")


# ==============================================================================
# Simulation / LAMMPS execution
# ==============================================================================

class SimulationError(ReboundMDError):
    """Raised when a LAMMPS simulation fails to run or produces bad output."""


class LAMMPSNotFoundError(SimulationError):
    """Raised when the 'lmp' (or 'lmp_mpi') executable cannot be located."""

    def __init__(self, executable: str = "lmp") -> None:
        super().__init__(
            f"LAMMPS executable '{executable}' not found in PATH. "
            "Install LAMMPS or set the LAMMPS_EXE environment variable.",
            executable=executable,
        )


class SimulationRuntimeError(SimulationError):
    """Raised when LAMMPS exits with a non-zero return code.

    Parameters
    ----------
    returncode : int
        Process return code.
    stderr : str
        Captured stderr output from LAMMPS.
    script : str
        Path to the input script that was run.
    """

    def __init__(self, returncode: int, stderr: str, script: str) -> None:
        self.returncode = returncode
        self.stderr = stderr
        self.script = script
        snippet = stderr.strip().splitlines()[-5:]  # last 5 lines
        super().__init__(
            f"LAMMPS terminated with code {returncode} while running '{script}'.\n"
            f"Last stderr lines:\n" + "\n".join(f"  {l}" for l in snippet),
            returncode=returncode,
            script=script,
        )


class TimeoutError(SimulationError):
    """Raised when a LAMMPS run exceeds the configured wall-time limit."""

    def __init__(self, elapsed: float, limit: float, script: str) -> None:
        super().__init__(
            f"Simulation '{script}' exceeded time limit "
            f"({elapsed:.1f}s > {limit:.1f}s).",
            elapsed=elapsed,
            limit=limit,
            script=script,
        )


# ==============================================================================
# Configuration / validation
# ==============================================================================

class ConfigError(ReboundMDError):
    """Raised when simulation configuration is invalid or self-inconsistent."""


class ValidationError(ConfigError):
    """Raised when a configuration parameter fails a physical sanity check.

    Parameters
    ----------
    param : str
        Name of the failing parameter.
    value : Any
        The invalid value.
    constraint : str
        Human-readable description of the constraint violated.
    """

    def __init__(self, param: str, value: Any, constraint: str) -> None:
        self.param = param
        self.value = value
        self.constraint = constraint
        super().__init__(
            f"Invalid value for '{param}': got {value!r}, expected {constraint}.",
            param=param,
            value=value,
            constraint=constraint,
        )


class MissingParameterError(ConfigError):
    """Raised when a required configuration key is absent."""

    def __init__(self, param: str, section: Optional[str] = None) -> None:
        loc = f" in section '{section}'" if section else ""
        super().__init__(
            f"Required parameter '{param}' is missing{loc}.",
            param=param,
            section=section,
        )


# ==============================================================================
# Physics / analysis
# ==============================================================================

class PhysicsError(ReboundMDError):
    """Raised when a physically impossible or undefined result is computed."""


class InsufficientDataError(PhysicsError):
    """Raised when there are too few data points to perform a calculation.

    Parameters
    ----------
    quantity : str
        What was being calculated (e.g. "coefficient of restitution").
    n_found : int
        Number of data points available.
    n_required : int
        Minimum required.
    """

    def __init__(self, quantity: str, n_found: int, n_required: int) -> None:
        super().__init__(
            f"Cannot compute '{quantity}': need at least {n_required} samples, "
            f"found {n_found}.",
            quantity=quantity,
            n_found=n_found,
            n_required=n_required,
        )


class NonPhysicalResultError(PhysicsError):
    """Raised when a computed quantity violates physical bounds.

    Examples: COR > 1 (super-elastic), negative kinetic energy.
    """

    def __init__(self, quantity: str, value: float, bounds: tuple[float, float]) -> None:
        lo, hi = bounds
        super().__init__(
            f"Non-physical value for '{quantity}': {value:.6g} "
            f"is outside allowed range [{lo}, {hi}].",
            quantity=quantity,
            value=value,
            bounds=bounds,
        )


class ConvergenceError(PhysicsError):
    """Raised when an iterative analysis algorithm fails to converge."""

    def __init__(self, algorithm: str, iterations: int, residual: float) -> None:
        super().__init__(
            f"'{algorithm}' did not converge after {iterations} iterations "
            f"(residual={residual:.3e}).",
            algorithm=algorithm,
            iterations=iterations,
            residual=residual,
        )
