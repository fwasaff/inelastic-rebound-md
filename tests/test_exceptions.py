"""
tests/test_exceptions.py
========================
Tests for the custom exception hierarchy.

Verifies:
- All exceptions are subclasses of ReboundMDError.
- __str__ produces informative messages with context.
- Structured attributes are accessible on exception instances.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from python.exceptions import (
    ConvergenceError,
    FileNotFoundError,
    InsufficientDataError,
    LAMMPSNotFoundError,
    MissingParameterError,
    NonPhysicalResultError,
    ParseError,
    ReboundMDError,
    RestartNotFoundError,
    SimulationRuntimeError,
    TimeoutError,
    ValidationError,
)


# ==============================================================================
# Hierarchy
# ==============================================================================

class TestHierarchy:

    @pytest.mark.parametrize("cls", [
        FileNotFoundError,
        ParseError,
        RestartNotFoundError,
        LAMMPSNotFoundError,
        SimulationRuntimeError,
        TimeoutError,
        ValidationError,
        MissingParameterError,
        InsufficientDataError,
        NonPhysicalResultError,
        ConvergenceError,
    ])
    def test_is_subclass_of_root(self, cls):
        assert issubclass(cls, ReboundMDError)

    def test_restart_is_file_not_found(self):
        from python.exceptions import FileNotFoundError as RFnF
        assert issubclass(RestartNotFoundError, RFnF)


# ==============================================================================
# FileNotFoundError
# ==============================================================================

class TestFileNotFoundError:

    def test_str_contains_path(self):
        exc = FileNotFoundError("/some/path/data.dat", role="trajectory")
        s = str(exc)
        assert "/some/path/data.dat" in s
        assert "trajectory" in s

    def test_path_attribute(self):
        exc = FileNotFoundError(Path("/tmp/x.dat"))
        assert exc.path == Path("/tmp/x.dat")

    def test_no_role(self):
        exc = FileNotFoundError("/file")
        assert "None" not in str(exc)


# ==============================================================================
# ParseError
# ==============================================================================

class TestParseError:

    def test_line_number_in_str(self):
        exc = ParseError("dump.lammpstrj", line_number=42, detail="bad token")
        assert "42" in str(exc)
        assert "bad token" in str(exc)

    def test_attributes(self):
        exc = ParseError("f.dat", line_number=10, raw_line="abc", detail="x")
        assert exc.filename == "f.dat"
        assert exc.line_number == 10
        assert exc.raw_line == "abc"


# ==============================================================================
# SimulationRuntimeError
# ==============================================================================

class TestSimulationRuntimeError:

    def test_returncode_in_str(self):
        exc = SimulationRuntimeError(returncode=1, stderr="ERROR: bad syntax\n", script="in.run")
        assert "1" in str(exc)
        assert "in.run" in str(exc)

    def test_last_lines_included(self):
        stderr = "\n".join([f"line {i}" for i in range(20)])
        exc = SimulationRuntimeError(1, stderr, "in.run")
        # Should include last few lines only
        assert "line 19" in str(exc)
        assert "line 0" not in str(exc)


# ==============================================================================
# ValidationError
# ==============================================================================

class TestValidationError:

    def test_message_format(self):
        exc = ValidationError("v0", -1.0, "positive (> 0)")
        s = str(exc)
        assert "v0" in s
        assert "-1.0" in s
        assert "positive" in s

    def test_attributes(self):
        exc = ValidationError("angle", 95.0, "in [0, 89.9]")
        assert exc.param == "angle"
        assert exc.value == 95.0
        assert "89.9" in exc.constraint


# ==============================================================================
# Physics errors
# ==============================================================================

class TestPhysicsErrors:

    def test_non_physical_result(self):
        exc = NonPhysicalResultError("COR", 1.5, (0.0, 1.0))
        s = str(exc)
        assert "COR" in s
        assert "1.5" in s
        assert "[0.0, 1.0]" in s

    def test_insufficient_data(self):
        exc = InsufficientDataError("COR calculation", n_found=2, n_required=10)
        assert "10" in str(exc)
        assert "2" in str(exc)

    def test_convergence_error(self):
        exc = ConvergenceError("Newton-Raphson", iterations=1000, residual=1e-3)
        assert "1000" in str(exc)
        assert "1e-3" in str(exc) or "0.001" in str(exc)


# ==============================================================================
# Context attributes
# ==============================================================================

class TestContextAttributes:

    def test_context_in_repr(self):
        exc = ValidationError("dt", -0.001, "positive")
        r = repr(exc)
        assert "ValidationError" in r

    def test_catch_as_parent(self):
        """All project exceptions can be caught as ReboundMDError."""
        with pytest.raises(ReboundMDError):
            raise ValidationError("x", 0, "positive")
