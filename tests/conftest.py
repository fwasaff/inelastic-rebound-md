"""
tests/conftest.py
=================
Shared pytest fixtures, helpers, and configuration.

All fixtures here are available to every test module without import.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from python.analysis import TrajectoryData
from python.config import EquilibrationConfig, ImpactConfig


# ==============================================================================
# Synthetic trajectory factory (shared across test modules)
# ==============================================================================

def make_trajectory(
    v0: float = 2.0,
    angle: float = 0.0,
    COR_true: float = 0.75,
    COR_t_true: float = 1.0,
    dt: float = 0.002,
    n_approach: int = 5000,
    n_contact: int = 2000,
    n_rebound: int = 3000,
    surface_top_z: float = 9.3,
    lj_cutoff: float = 2.5,
    T_init: float = 0.1,
    dT: float = 0.005,
    noise: float = 0.001,
    seed: int = 42,
) -> TrajectoryData:
    """Generate a physically plausible synthetic trajectory.

    The projectile:
    - Approaches at velocity (-v0_x, -v0_z) (approach phase).
    - Penetrates and rebounds: vz reverses with COR, vx changes with COR_t.
    - Departs at constant velocity (vf_x, +vf_z) (rebound phase).

    Parameters
    ----------
    v0 : float
        Incident speed [sigma/tau].
    angle : float
        Impact angle from surface normal [degrees].
    COR_true : float
        Ground-truth normal coefficient of restitution.
    COR_t_true : float
        Ground-truth tangential restitution coefficient vf_x/v0_x (typically ~1).
    dt : float
        Time step [tau].
    n_approach, n_contact, n_rebound : int
        Number of frames in each phase.
    surface_top_z : float
        Z-coordinate of unperturbed surface top [sigma].
    lj_cutoff : float
        LJ cutoff used for contact detection [sigma].
    T_init : float
        Initial surface temperature [eps/kB].
    dT : float
        Surface temperature rise during impact [eps/kB].
    noise : float
        Gaussian noise amplitude on z.
    seed : int
        NumPy random seed for reproducibility.
    """
    rng = np.random.default_rng(seed)
    n_total = n_approach + n_contact + n_rebound

    v0_z = v0 * math.cos(math.radians(angle))
    v0_x = v0 * math.sin(math.radians(angle))
    vf_z = COR_true * v0_z
    vf_x = COR_t_true * v0_x

    # --- Approach ---
    z_start = surface_top_z + v0_z * n_approach * dt + 5.0
    z_approach  = np.linspace(z_start, surface_top_z + 0.01, n_approach)
    vz_approach = np.full(n_approach, -v0_z)
    vx_approach = np.full(n_approach, v0_x)
    ke_approach = np.full(n_approach, 0.5 * (v0_z ** 2 + v0_x ** 2))

    # --- Contact ---
    t_c = np.linspace(0.0, n_contact * dt, n_contact)
    delta_max = min(0.3, lj_cutoff * 0.9)
    z_contact  = surface_top_z - delta_max * np.sin(math.pi * t_c / (n_contact * dt))
    vz_contact = np.linspace(-v0_z, vf_z, n_contact)
    vx_contact = np.linspace(v0_x, vf_x, n_contact)
    ke_contact = 0.5 * (vz_contact ** 2 + vx_contact ** 2)

    # --- Rebound ---
    z_rebound  = surface_top_z + vf_z * np.arange(n_rebound) * dt
    vz_rebound = np.full(n_rebound, vf_z)
    vx_rebound = np.full(n_rebound, vf_x)
    ke_rebound = np.full(n_rebound, 0.5 * (vf_z ** 2 + vf_x ** 2))

    # Assemble
    z  = np.concatenate([z_approach,  z_contact,  z_rebound])
    vz = np.concatenate([vz_approach, vz_contact, vz_rebound])
    vx = np.concatenate([vx_approach, vx_contact, vx_rebound])
    ke = np.concatenate([ke_approach, ke_contact, ke_rebound])

    z  += rng.normal(0.0, noise, n_total)
    vz += rng.normal(0.0, noise * 0.1, n_total)
    vx += rng.normal(0.0, noise * 0.1, n_total)

    time = np.arange(n_total, dtype=float) * dt

    T_surf = np.full(n_total, T_init)
    T_surf[n_approach : n_approach + n_contact] += np.linspace(0.0, dT, n_contact)
    T_surf[n_approach + n_contact :] = T_init + dT

    d_contact = z - surface_top_z

    return TrajectoryData(
        time=time,
        z=z,
        vx=vx,
        vz=vz,
        ke_proj=ke,
        T_surf=T_surf,
        d_contact=d_contact,
    )


def write_trajectory_file(data: TrajectoryData, path: Path) -> None:
    """Write ``TrajectoryData`` to a ``rebound_data.dat`` text file."""
    with path.open("w") as fh:
        fh.write("# time z vx vz ke_proj T_surf d_contact\n")
        for i in range(len(data)):
            fh.write(
                f"{data.time[i]:.8f} {data.z[i]:.8f} {data.vx[i]:.8f} "
                f"{data.vz[i]:.8f} {data.ke_proj[i]:.8f} "
                f"{data.T_surf[i]:.8f} {data.d_contact[i]:.8f}\n"
            )


# ==============================================================================
# Pytest fixtures
# ==============================================================================

@pytest.fixture(scope="session")
def default_trajectory() -> TrajectoryData:
    """Synthetic trajectory with v0=2.0, COR=0.75 (session-scoped, shared)."""
    return make_trajectory()


@pytest.fixture
def trajectory_file(tmp_path, default_trajectory) -> Path:
    """Write the default trajectory to a temp file and return its path."""
    p = tmp_path / "rebound_data.dat"
    write_trajectory_file(default_trajectory, p)
    return p


@pytest.fixture
def default_equil_config() -> EquilibrationConfig:
    """Small equilibration config for fast tests."""
    return EquilibrationConfig(Nx=4, Ny=4, Nz=4, steps_equil=100, steps_prod=100)


@pytest.fixture
def default_impact_config() -> ImpactConfig:
    return ImpactConfig(v0=2.0, angle=0.0, steps_approach=100, steps_impact=500)


# ==============================================================================
# Pytest configuration
# ==============================================================================

def pytest_configure(config):
    """Register custom markers."""
    config.addinivalue_line(
        "markers", "lammps: tests that require a LAMMPS installation (skip if absent)"
    )
    config.addinivalue_line(
        "markers", "slow: tests that take more than a few seconds"
    )


def pytest_collection_modifyitems(config, items):
    """Auto-skip lammps-marked tests if 'lmp' is not on PATH."""
    import shutil
    lmp_available = shutil.which("lmp") is not None or shutil.which("lmp_mpi") is not None
    skip_lammps = pytest.mark.skip(reason="LAMMPS not installed (lmp not on PATH)")
    for item in items:
        if "lammps" in item.keywords and not lmp_available:
            item.add_marker(skip_lammps)
