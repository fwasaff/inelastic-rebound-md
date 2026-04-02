"""
tests/test_config.py
====================
Dedicated tests for the metaclass-based configuration system.

Covers:
- FieldDescriptor: coercion, validation, default, REQUIRED sentinel
- ConfigMeta:      __init__, __repr__, __eq__, to_dict/from_dict, JSON round-trip
- ConfigRegistry:  __init_subclass__ registration, from_json_auto
- All concrete config classes
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from python.config import (
    ConfigRegistry,
    EquilibrationConfig,
    FieldDescriptor,
    ImpactConfig,
    ParametricSweepConfig,
)
from python.exceptions import MissingParameterError, ValidationError


# ==============================================================================
# FieldDescriptor
# ==============================================================================

class TestFieldDescriptor:

    def test_coercion_int_to_float(self):
        cfg = ImpactConfig(v0=3)    # int 3 should coerce to float 3.0
        assert isinstance(cfg.v0, float)
        assert cfg.v0 == 3.0

    def test_invalid_type_raises_validation(self):
        with pytest.raises(ValidationError):
            ImpactConfig(v0="fast")

    def test_default_returned_when_not_set(self):
        cfg = ImpactConfig()
        assert cfg.v0 == 2.0      # default
        assert cfg.angle == 0.0   # default

    def test_required_field_raises_without_default(self):
        """FieldDescriptor.REQUIRED sentinel should trigger MissingParameterError."""
        from python.config import ConfigMeta, FieldDescriptor, BaseConfig

        class _Strict(BaseConfig):
            mandatory = FieldDescriptor(float, FieldDescriptor.REQUIRED,
                                         doc="no default")

        with pytest.raises(MissingParameterError):
            _Strict()


# ==============================================================================
# ConfigMeta – generated methods
# ==============================================================================

class TestConfigMeta:

    def test_repr_contains_field_values(self):
        cfg = ImpactConfig(v0=1.5, angle=30.0)
        r = repr(cfg)
        assert "1.5" in r
        assert "30.0" in r
        assert "ImpactConfig" in r

    def test_eq_same_values(self):
        assert ImpactConfig(v0=2.0) == ImpactConfig(v0=2.0)

    def test_eq_different_values(self):
        assert ImpactConfig(v0=2.0) != ImpactConfig(v0=3.0)

    def test_eq_different_types(self):
        cfg = ImpactConfig()
        result = cfg.__eq__("not a config")
        assert result is NotImplemented

    def test_to_dict_keys(self):
        cfg = ImpactConfig()
        d = cfg.to_dict()
        assert "v0" in d
        assert "angle" in d
        assert "T_target" in d

    def test_from_dict_roundtrip(self):
        cfg = ImpactConfig(v0=3.1, angle=45.0)
        cfg2 = ImpactConfig.from_dict(cfg.to_dict())
        assert cfg == cfg2

    def test_unknown_kwarg_raises_type_error(self):
        with pytest.raises(TypeError, match="unexpected keyword"):
            ImpactConfig(nonexistent=99)

    def test_inheritance_collects_parent_fields(self):
        """EquilibrationConfig inherits nothing from ImpactConfig but both
        inherit from BaseConfig; their _fields should be independent."""
        equil_fields = set(EquilibrationConfig._fields)
        impact_fields = set(ImpactConfig._fields)
        assert "T_target" in equil_fields
        assert "v0" not in equil_fields
        assert "v0" in impact_fields


# ==============================================================================
# JSON serialisation
# ==============================================================================

class TestJSONSerialisation:

    def test_to_json_contains_type_key(self, tmp_path):
        cfg = ImpactConfig(v0=2.5)
        path = tmp_path / "cfg.json"
        cfg.to_json(path)
        raw = json.loads(path.read_text())
        assert raw["_type"] == "ImpactConfig"

    def test_from_json_auto_dispatch(self, tmp_path):
        cfg = EquilibrationConfig(T_target=0.2, Nx=6)
        path = tmp_path / "equil.json"
        cfg.to_json(path)
        loaded = ConfigRegistry.from_json_auto(path)
        assert isinstance(loaded, EquilibrationConfig)
        assert loaded == cfg

    def test_from_json_auto_missing_type_raises(self, tmp_path):
        path = tmp_path / "notype.json"
        path.write_text('{"v0": 2.0}')
        with pytest.raises(ValueError, match="_type"):
            ConfigRegistry.from_json_auto(path)

    def test_from_json_auto_unknown_type_raises(self, tmp_path):
        path = tmp_path / "unknown.json"
        path.write_text('{"_type": "NoSuchConfig", "v0": 2.0}')
        with pytest.raises(KeyError):
            ConfigRegistry.from_json_auto(path)


# ==============================================================================
# EquilibrationConfig
# ==============================================================================

class TestEquilibrationConfig:

    def test_total_atoms_formula(self):
        cfg = EquilibrationConfig(Nx=4, Ny=4, Nz=4)
        assert cfg.total_atoms == 4 * 4 * 4 * 4  # 4 atoms / unit cell

    def test_slab_height(self):
        cfg = EquilibrationConfig(Nz=6, a0=1.5496)
        assert abs(cfg.slab_height - 6 * 1.5496) < 1e-10

    def test_negative_temperature_raises(self):
        with pytest.raises(ValidationError):
            EquilibrationConfig(T_target=-0.1)

    def test_zero_Nz_raises(self):
        with pytest.raises(ValidationError):
            EquilibrationConfig(Nz=0)


# ==============================================================================
# ImpactConfig
# ==============================================================================

class TestImpactConfig:

    @pytest.mark.parametrize("angle,expected_vz,expected_vx", [
        (0.0,  -2.0,  0.0),
        (90.0,  0.0,  2.0),
        (60.0, -1.0,  math.sqrt(3)),
        (30.0, -math.sqrt(3), 1.0),
    ])
    def test_velocity_components(self, angle, expected_vz, expected_vx):
        cfg = ImpactConfig(v0=2.0, angle=angle)
        assert abs(cfg.vz0 - expected_vz) < 1e-9, f"vz0 wrong for angle={angle}"
        assert abs(cfg.vx0 - expected_vx) < 1e-9, f"vx0 wrong for angle={angle}"

    def test_total_time(self):
        cfg = ImpactConfig(dt=0.002, steps_approach=1000, steps_impact=2000)
        assert abs(cfg.total_time - 6.0) < 1e-12

    def test_angle_90_raises(self):
        """angle=90 means grazing incidence; vz0=0 which is degenerate."""
        with pytest.raises(ValidationError):
            ImpactConfig(angle=90.0)  # constraint is < 89.9


# ==============================================================================
# ParametricSweepConfig
# ==============================================================================

class TestParametricSweepConfig:

    def test_v0_values_linspace(self):
        cfg = ParametricSweepConfig(v0_min=1.0, v0_max=5.0, v0_steps=5)
        vals = cfg.v0_values
        assert len(vals) == 5
        assert abs(vals[0] - 1.0) < 1e-10
        assert abs(vals[-1] - 5.0) < 1e-10
        for i in range(1, len(vals)):
            assert vals[i] > vals[i - 1]

    def test_v0_values_single_step(self):
        cfg = ParametricSweepConfig(v0_min=2.5, v0_max=4.0, v0_steps=1)
        assert cfg.v0_values == [2.5]

    def test_angles_parse(self):
        cfg = ParametricSweepConfig(angle_list="0.0, 15.0, 30.0, 45.0")
        assert cfg.angles == [0.0, 15.0, 30.0, 45.0]
