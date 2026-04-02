"""
tests/test_analysis.py
======================
Unit tests for the analysis module using synthetic trajectory data.

Tests are designed to run WITHOUT a LAMMPS installation. All
simulated data is generated in-memory using numpy.

Run with:   pytest tests/ -v
"""

from __future__ import annotations

import math
import tempfile
from pathlib import Path

import numpy as np
import pytest

# We test only pure-Python / numpy modules — no LAMMPS required.
from python.analysis import (
    ParametricAnalyzer,
    ReboundAnalyzer,
    ReboundDataParser,
    ReboundResult,
    TrajectoryData,
)
from python.config import (
    ConfigRegistry,
    EquilibrationConfig,
    ImpactConfig,
    ParametricSweepConfig,
)
from python.exceptions import (
    InsufficientDataError,
    MissingParameterError,
    NonPhysicalResultError,
    ParseError,
    ValidationError,
)


# ==============================================================================
# Fixtures – synthetic trajectory generation
# ==============================================================================

def _make_synthetic_trajectory(
    v0: float = 2.0,
    COR_true: float = 0.75,
    dt: float = 0.002,
    n_approach: int = 5000,
    n_contact: int = 2000,
    n_rebound: int = 3000,
    surface_top_z: float = 9.3,
    T_init: float = 0.1,
    dT: float = 0.005,
    noise: float = 0.001,
) -> TrajectoryData:
    """Generate a physically plausible synthetic trajectory.

    The projectile:
      - Approaches at constant velocity -v0.
      - Contacts (parabolic z, velocity reversal with COR).
      - Rebounds at constant velocity +COR*v0.
    """
    rng = np.random.default_rng(42)
    n_total = n_approach + n_contact + n_rebound

    # Approach phase: z decreases linearly from z_start to surface_top
    z_start = surface_top_z + v0 * n_approach * dt + 5.0
    z_approach = np.linspace(z_start, surface_top_z + 0.01, n_approach)
    vz_approach = np.full(n_approach, -v0)
    ke_approach = np.full(n_approach, 0.5 * v0 ** 2)

    # Contact phase: parabolic z (simple penetration model), velocity transitions
    t_contact = np.linspace(0, n_contact * dt, n_contact)
    delta_max = 0.3  # penetration depth
    z_contact = surface_top_z - delta_max * np.sin(math.pi * t_contact / (n_contact * dt))
    # Velocity linearly reverses from -v0 to +COR*v0
    vz_contact = np.linspace(-v0, COR_true * v0, n_contact)
    ke_contact = 0.5 * vz_contact ** 2

    # Rebound phase: z increases linearly
    vf = COR_true * v0
    z_rebound = surface_top_z + vf * np.arange(n_rebound) * dt
    vz_rebound = np.full(n_rebound, vf)
    ke_rebound = np.full(n_rebound, 0.5 * vf ** 2)

    # Concatenate
    z  = np.concatenate([z_approach, z_contact, z_rebound])
    vz = np.concatenate([vz_approach, vz_contact, vz_rebound])
    ke = np.concatenate([ke_approach, ke_contact, ke_rebound])

    # Add Gaussian noise
    z  += rng.normal(0, noise, n_total)
    vz += rng.normal(0, noise * 0.1, n_total)

    # Time
    time = np.arange(n_total) * dt

    # Surface temperature: rises during contact, levels off
    T_surf = np.full(n_total, T_init)
    T_surf[n_approach:n_approach + n_contact] += np.linspace(0, dT, n_contact)
    T_surf[n_approach + n_contact:] = T_init + dT

    # Distance from surface top
    d_contact = z - surface_top_z

    return TrajectoryData(
        time=time,
        z=z,
        vz=vz,
        ke_proj=ke,
        T_surf=T_surf,
        d_contact=d_contact,
    )


def _write_trajectory_file(data: TrajectoryData, path: Path) -> None:
    """Write a TrajectoryData to a rebound_data.dat text file."""
    with path.open("w") as fh:
        fh.write("# time z vz ke_proj T_surf d_contact\n")
        for i in range(len(data)):
            fh.write(
                f"{data.time[i]:.6f} {data.z[i]:.6f} {data.vz[i]:.6f} "
                f"{data.ke_proj[i]:.6f} {data.T_surf[i]:.6f} {data.d_contact[i]:.6f}\n"
            )


@pytest.fixture
def synthetic_data() -> TrajectoryData:
    return _make_synthetic_trajectory()


@pytest.fixture
def rebound_data_file(tmp_path, synthetic_data) -> Path:
    p = tmp_path / "rebound_data.dat"
    _write_trajectory_file(synthetic_data, p)
    return p


# ==============================================================================
# Tests: ReboundDataParser
# ==============================================================================

class TestReboundDataParser:

    def test_parse_valid_file(self, rebound_data_file):
        data = ReboundDataParser.parse(rebound_data_file)
        assert isinstance(data, TrajectoryData)
        assert len(data) == 10000   # 5000 + 2000 + 3000

    def test_parse_column_counts(self, rebound_data_file):
        data = ReboundDataParser.parse(rebound_data_file)
        assert data.time.shape == data.z.shape == data.vz.shape

    def test_missing_file_raises(self, tmp_path):
        from python.exceptions import FileNotFoundError as ReboundFNF
        with pytest.raises(ReboundFNF):
            ReboundDataParser.parse(tmp_path / "nonexistent.dat")

    def test_bad_column_count_raises(self, tmp_path):
        p = tmp_path / "bad.dat"
        p.write_text("# time z vz ke_proj T_surf d_contact\n1.0 2.0 3.0\n")
        with pytest.raises(ParseError):
            ReboundDataParser.parse(p)

    def test_too_few_rows_raises(self, tmp_path):
        p = tmp_path / "short.dat"
        lines = ["# time z vz ke_proj T_surf d_contact\n"]
        for i in range(5):
            lines.append(f"{i} {i} {i} {i} 0.1 {i}\n")
        p.write_text("".join(lines))
        with pytest.raises(InsufficientDataError):
            ReboundDataParser.parse(p)

    def test_non_numeric_data_raises(self, tmp_path):
        p = tmp_path / "bad_numeric.dat"
        p.write_text("# time z vz ke_proj T_surf d_contact\nfoo bar baz qux quux quuz\n")
        with pytest.raises(ParseError):
            ReboundDataParser.parse(p)


# ==============================================================================
# Tests: ReboundAnalyzer
# ==============================================================================

class TestReboundAnalyzer:

    def test_load_and_analyze(self, rebound_data_file):
        analyzer = ReboundAnalyzer(v0=2.0)
        analyzer.load(rebound_data_file)
        result = analyzer.compute_all()
        assert isinstance(result, ReboundResult)

    def test_cor_in_physical_range(self, rebound_data_file):
        analyzer = ReboundAnalyzer(v0=2.0)
        analyzer.load(rebound_data_file)
        cor = analyzer.coefficient_of_restitution
        assert 0.0 <= cor <= 1.0, f"COR={cor} out of range"

    def test_cor_approximately_correct(self, rebound_data_file):
        """COR should be close to the synthetic ground truth of 0.75."""
        analyzer = ReboundAnalyzer(v0=2.0)
        analyzer.load(rebound_data_file)
        cor = analyzer.coefficient_of_restitution
        assert abs(cor - 0.75) < 0.05, f"COR={cor:.4f}, expected ~0.75"

    def test_energy_dissipated_positive(self, rebound_data_file):
        analyzer = ReboundAnalyzer(v0=2.0)
        analyzer.load(rebound_data_file)
        assert analyzer.energy_dissipated > 0

    def test_contact_time_positive(self, rebound_data_file):
        analyzer = ReboundAnalyzer(v0=2.0)
        analyzer.load(rebound_data_file)
        assert analyzer.contact_time > 0

    def test_surface_temperature_rise(self, rebound_data_file):
        analyzer = ReboundAnalyzer(v0=2.0, T_surface_init=0.1)
        analyzer.load(rebound_data_file)
        dT = analyzer.surface_temperature_rise
        assert dT >= 0, f"Surface temperature decreased: dT={dT}"

    def test_require_data_raises_before_load(self):
        analyzer = ReboundAnalyzer(v0=2.0)
        with pytest.raises(RuntimeError, match="load()"):
            _ = analyzer.coefficient_of_restitution

    def test_non_physical_cor_raises(self, tmp_path):
        """Force a super-elastic result to trigger NonPhysicalResultError."""
        data = _make_synthetic_trajectory(v0=2.0, COR_true=1.1)  # COR > 1
        p = tmp_path / "super_elastic.dat"
        _write_trajectory_file(data, p)

        analyzer = ReboundAnalyzer(v0=2.0)
        analyzer.load(p)
        with pytest.raises(NonPhysicalResultError):
            _ = analyzer.coefficient_of_restitution

    def test_caching(self, rebound_data_file):
        """Repeated calls to cached properties should return identical objects."""
        analyzer = ReboundAnalyzer(v0=2.0)
        analyzer.load(rebound_data_file)
        cor1 = analyzer.coefficient_of_restitution
        cor2 = analyzer.coefficient_of_restitution
        assert cor1 == cor2


# ==============================================================================
# Tests: Config metaclass
# ==============================================================================

class TestConfig:

    def test_impact_config_defaults(self):
        cfg = ImpactConfig()
        assert cfg.v0 == 2.0
        assert cfg.angle == 0.0

    def test_impact_config_custom(self):
        cfg = ImpactConfig(v0=3.5, angle=15.0)
        assert cfg.v0 == 3.5
        assert cfg.angle == 15.0

    def test_validation_negative_v0(self):
        with pytest.raises(ValidationError):
            ImpactConfig(v0=-1.0)

    def test_validation_angle_out_of_range(self):
        with pytest.raises(ValidationError):
            ImpactConfig(angle=95.0)

    def test_vz0_normal_incidence(self):
        cfg = ImpactConfig(v0=2.0, angle=0.0)
        assert abs(cfg.vz0 - (-2.0)) < 1e-12

    def test_vz0_oblique_incidence(self):
        cfg = ImpactConfig(v0=2.0, angle=60.0)
        assert abs(cfg.vz0 - (-1.0)) < 1e-6   # cos(60°) = 0.5

    def test_vx0_oblique_incidence(self):
        cfg = ImpactConfig(v0=2.0, angle=30.0)
        assert abs(cfg.vx0 - 1.0) < 1e-6      # sin(30°) = 0.5 -> vx = 1.0

    def test_to_dict_round_trip(self):
        cfg = ImpactConfig(v0=1.5, angle=10.0)
        d = cfg.to_dict()
        cfg2 = ImpactConfig.from_dict(d)
        assert cfg == cfg2

    def test_to_json_round_trip(self, tmp_path):
        cfg = ImpactConfig(v0=2.5, angle=30.0)
        path = tmp_path / "impact.json"
        cfg.to_json(path)
        cfg2 = ConfigRegistry.from_json_auto(path)
        assert cfg == cfg2

    def test_unknown_kwarg_raises(self):
        with pytest.raises(TypeError, match="unexpected keyword"):
            ImpactConfig(nonexistent_param=99)

    def test_equil_config_total_atoms(self):
        cfg = EquilibrationConfig(Nx=4, Ny=4, Nz=4)
        assert cfg.total_atoms == 4 * 4 * 4 * 4  # 4 atoms/unit cell

    def test_parametric_sweep_v0_values(self):
        cfg = ParametricSweepConfig(v0_min=1.0, v0_max=3.0, v0_steps=3)
        vals = cfg.v0_values
        assert len(vals) == 3
        assert abs(vals[0] - 1.0) < 1e-10
        assert abs(vals[2] - 3.0) < 1e-10

    def test_parametric_sweep_angles(self):
        cfg = ParametricSweepConfig(angle_list="0.0,30.0,60.0")
        assert cfg.angles == [0.0, 30.0, 60.0]


# ==============================================================================
# Tests: ParametricAnalyzer
# ==============================================================================

class TestParametricAnalyzer:

    def _make_results(self, v0_list, cor_list):
        return [
            ReboundResult(
                v0_z=v0, vf_z=v0 * cor, COR=cor,
                delta_E=0.5 * v0 ** 2 * (1 - cor ** 2),
                contact_time=0.1, max_penetration=0.2,
                dT_surface=0.001
            )
            for v0, cor in zip(v0_list, cor_list)
        ]

    def test_cor_vs_v0_basic(self):
        v0s = [1.0, 2.0, 3.0, 1.0, 2.0, 3.0]
        cors = [0.9, 0.8, 0.7, 0.88, 0.82, 0.72]
        pa = ParametricAnalyzer(self._make_results(v0s, cors))
        v0_u, cor_m, cor_s = pa.cor_vs_v0()
        assert len(v0_u) == 3
        assert cor_m[0] > cor_m[2]   # COR decreases with v0

    def test_fit_power_law(self):
        # Generate COR = 0.9 * v0^(-0.1)  (weak velocity dependence)
        v0s = np.linspace(1.0, 5.0, 10).tolist() * 2
        cors = [0.9 * v ** (-0.1) for v in v0s]
        pa = ParametricAnalyzer(self._make_results(v0s, cors))
        A, beta = pa.fit_cor_power_law()
        assert abs(A - 0.9) < 0.05
        assert abs(beta - 0.1) < 0.02

    def test_insufficient_data_raises(self):
        pa = ParametricAnalyzer()
        with pytest.raises(InsufficientDataError):
            pa.cor_vs_v0()
