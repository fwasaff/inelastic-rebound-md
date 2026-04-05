"""
analysis.py
===========
Physical analysis of inelastic rebound simulations.

Advanced Python features used:
- Decorators: ``@timed``, ``@cached_property_invalidated``, ``@require_data``
- Properties: lazy-loaded, auto-invalidated on data reload
- Dataclasses: structured result containers
- __slots__: memory-efficient storage for large trajectory arrays
- Type annotations throughout

Key physics computed:
- Coefficient of Restitution (COR):  e = |v_f_z| / |v_0_z|
- Energy Accommodation Coefficient (EAC):  α_E = 1 - E_f/E_inc  (general)
- Normal EAC:       α_n = 1 - e²  = ΔE_n / E_inc_n
- Tangential EAC:   α_t = 1 - (v_f_x/v_0_x)²  (oblique impacts)
- Reduced energy:   E* = ½ m v₀²  (= ½ v₀² in LJ units where m=ε=1)
- Energy dissipation:  ΔE = E_inc - E_fin  [eps]
- Contact time:  τ_c (time projectile is within LJ cutoff of surface)
- Maximum penetration depth:  δ_max
- Surface temperature rise:  ΔT = T_final - T_initial
"""

from __future__ import annotations

import functools
import logging
import time
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional, TypeVar

import numpy as np

from .exceptions import (
    FileNotFoundError,
    InsufficientDataError,
    NonPhysicalResultError,
    ParseError,
)

logger = logging.getLogger(__name__)

F = TypeVar("F", bound=Callable)


# ==============================================================================
# Decorators
# ==============================================================================

def timed(func: F) -> F:
    """Decorator: log the wall-clock time taken by a method or function."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        t0 = time.perf_counter()
        result = func(*args, **kwargs)
        elapsed = time.perf_counter() - t0
        logger.debug("%s completed in %.4f s", func.__qualname__, elapsed)
        return result
    return wrapper  # type: ignore[return-value]


def require_data(func: F) -> F:
    """Decorator: raise RuntimeError if the analyzer has no data loaded yet."""
    @functools.wraps(func)
    def wrapper(self: "ReboundAnalyzer", *args, **kwargs):
        if self._data is None:
            raise RuntimeError(
                f"'{func.__name__}' called before data was loaded. "
                "Call ReboundAnalyzer.load() first."
            )
        return func(self, *args, **kwargs)
    return wrapper  # type: ignore[return-value]


def validate_physical(quantity: str, lo: float, hi: float) -> Callable[[F], F]:
    """Parametric decorator factory: check that a numeric result is in [lo, hi].

    Usage
    -----
    @validate_physical("COR", 0.0, 1.0)
    def coefficient_of_restitution(self) -> float:
        ...
    """
    def decorator(func: F) -> F:
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            result = func(*args, **kwargs)
            if not (lo <= result <= hi):
                raise NonPhysicalResultError(quantity, result, (lo, hi))
            return result
        return wrapper  # type: ignore[return-value]
    return decorator


def cached(func: F) -> F:
    """Decorator: cache the result on the instance as ``_cache_<name>``.

    Invalidated automatically when ``ReboundAnalyzer._data`` is reassigned.
    """
    cache_attr = f"_cache_{func.__name__}"

    @functools.wraps(func)
    def wrapper(self: "ReboundAnalyzer", *args, **kwargs):
        if not hasattr(self, cache_attr):
            setattr(self, cache_attr, func(self, *args, **kwargs))
        return getattr(self, cache_attr)
    return wrapper  # type: ignore[return-value]


# ==============================================================================
# Data containers (dataclasses)
# ==============================================================================

@dataclass(frozen=True)
class ReboundResult:
    """Immutable result container for a single impact event.

    Attributes
    ----------
    v0_z : float
        Incident normal velocity [sigma/tau].
    vf_z : float
        Final normal velocity after rebound [sigma/tau] (positive = upward).
    COR : float
        Coefficient of restitution = |vf_z| / |v0_z|.
    delta_E : float
        Energy dissipated = E_inc - E_fin  [eps].
    contact_time : float
        Duration of contact phase [tau].
    max_penetration : float
        Maximum penetration below surface equilibrium [sigma].
    dT_surface : float
        Surface temperature rise during impact [eps/kB].
    angle : float
        Impact angle from normal [degrees].
    T_surface : float
        Surface temperature [eps/kB].
    EAC : float
        Energy Accommodation Coefficient α_E = ΔE/E_inc  (total, all components).
    EAC_normal : float
        Normal EAC α_n = 1 - COR²  (fraction of normal KE transferred).
    EAC_tangential : float
        Tangential EAC α_t = 1 - (vf_x/v0_x)²  (0 for normal incidence).
    E_star : float
        Reduced incident energy E* = ½v₀²  [eps]  (in LJ units m=ε=1).
    """

    v0_z: float
    vf_z: float
    COR: float
    delta_E: float
    contact_time: float
    max_penetration: float
    dT_surface: float
    angle: float = 0.0
    T_surface: float = 0.1
    EAC: float = 0.0
    EAC_normal: float = 0.0
    EAC_tangential: float = 0.0
    E_star: float = 0.0


@dataclass
class TrajectoryData:
    """Parsed rebound_data.dat contents as numpy arrays.

    All arrays have the same length (number of recorded timesteps).

    Attributes
    ----------
    time     : 1-D array, simulation time [tau]
    z        : 1-D array, projectile z-position [sigma]
    vx       : 1-D array, projectile x-velocity (tangential) [sigma/tau]
    vz       : 1-D array, projectile z-velocity (normal) [sigma/tau]
    ke_proj  : 1-D array, projectile total kinetic energy [eps]
    T_surf   : 1-D array, surface temperature [eps/kB]
    d_contact: 1-D array, signed distance from surface top [sigma]
    """

    __slots__ = ("time", "z", "vx", "vz", "ke_proj", "T_surf", "d_contact")

    time:      np.ndarray
    z:         np.ndarray
    vx:        np.ndarray
    vz:        np.ndarray
    ke_proj:   np.ndarray
    T_surf:    np.ndarray
    d_contact: np.ndarray

    def __len__(self) -> int:
        return len(self.time)


# ==============================================================================
# Parser
# ==============================================================================

class ReboundDataParser:
    """Parses the ``rebound_data.dat`` file written by LAMMPS ``fix print``.

    The file has one header line (``# time z vx vz ke_proj T_surf d_contact``)
    followed by whitespace-delimited numeric data.
    """

    EXPECTED_COLUMNS = ("time", "z", "vx", "vz", "ke_proj", "T_surf", "d_contact")
    N_COLS = len(EXPECTED_COLUMNS)

    @classmethod
    @timed
    def parse(cls, path: Path | str) -> TrajectoryData:
        """Parse ``rebound_data.dat`` and return a ``TrajectoryData`` instance.

        Parameters
        ----------
        path : Path | str
            Path to the data file.

        Raises
        ------
        FileNotFoundError
            If the file does not exist.
        ParseError
            If the file format is unexpected.
        InsufficientDataError
            If fewer than 10 data rows are found.
        """
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(path, role="rebound tracking data")

        rows: list[list[float]] = []
        with path.open() as fh:
            for lineno, raw in enumerate(fh, start=1):
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                tokens = line.split()
                if len(tokens) != cls.N_COLS:
                    raise ParseError(
                        path.name,
                        line_number=lineno,
                        raw_line=raw.rstrip(),
                        detail=f"expected {cls.N_COLS} columns, got {len(tokens)}",
                    )
                try:
                    rows.append([float(t) for t in tokens])
                except ValueError as exc:
                    raise ParseError(
                        path.name,
                        line_number=lineno,
                        raw_line=raw.rstrip(),
                        detail=str(exc),
                    ) from exc

        if len(rows) < 10:
            raise InsufficientDataError("trajectory parsing", len(rows), 10)

        arr = np.array(rows, dtype=float).T
        return TrajectoryData(
            time=arr[0],
            z=arr[1],
            vx=arr[2],
            vz=arr[3],
            ke_proj=arr[4],
            T_surf=arr[5],
            d_contact=arr[6],
        )


# ==============================================================================
# Main analyzer
# ==============================================================================

class ReboundAnalyzer:
    """Computes physical observables from a parsed rebound trajectory.

    Parameters
    ----------
    v0 : float
        Incident speed [sigma/tau] (used to verify COR bounds).
    angle : float
        Impact angle from normal [degrees].
    surface_top_z : float
        Z-coordinate of the surface top [sigma] (used for contact detection).
    lj_cutoff : float
        LJ cutoff distance [sigma] (contact criterion: d < lj_cutoff).
    T_surface_init : float
        Initial surface temperature before impact [eps/kB].
    """

    def __init__(
        self,
        v0: float,
        angle: float = 0.0,
        surface_top_z: float = 9.3,
        lj_cutoff: float = 2.5,
        T_surface_init: float = 0.1,
    ) -> None:
        self.v0 = v0
        self.angle = angle
        self.surface_top_z = surface_top_z
        self.lj_cutoff = lj_cutoff
        self.T_surface_init = T_surface_init

        self._data: Optional[TrajectoryData] = None
        self._invalidate_cache()

    def _invalidate_cache(self) -> None:
        """Remove all cached computation results."""
        for attr in list(vars(self)):
            if attr.startswith("_cache_"):
                delattr(self, attr)

    @timed
    def load(self, path: Path | str) -> "ReboundAnalyzer":
        """Load trajectory data from ``rebound_data.dat``.

        Returns ``self`` for method chaining.
        """
        self._data = ReboundDataParser.parse(path)
        self._invalidate_cache()
        logger.info("Loaded %d trajectory frames from %s", len(self._data), path)
        return self

    # --- Core physics ---------------------------------------------------------

    @property
    @require_data
    def v0_z(self) -> float:
        """Incident normal velocity (magnitude, positive) from trajectory."""
        import math
        return self.v0 * math.cos(math.radians(self.angle))

    @require_data
    @cached
    def _contact_mask(self) -> np.ndarray:
        """Boolean mask: True when projectile is within LJ cutoff of surface."""
        return self._data.d_contact < self.lj_cutoff

    @require_data
    @cached
    def _approach_end_idx(self) -> int:
        """Index of last frame before contact begins."""
        mask = self._contact_mask()
        contact_indices = np.where(mask)[0]
        if len(contact_indices) == 0:
            raise InsufficientDataError("contact detection", 0, 1)
        return int(contact_indices[0]) - 1

    @require_data
    @cached
    def _rebound_start_idx(self) -> int:
        """Index of first frame after contact ends (projectile moving away).

        If the projectile never clearly leaves the contact zone (sticking
        regime), fall back to the last 50 frames so callers get valid data
        rather than an out-of-bounds empty window.
        """
        mask = self._contact_mask()
        contact_indices = np.where(mask)[0]
        if len(contact_indices) == 0:
            raise InsufficientDataError("rebound detection", 0, 1)
        last_contact = int(contact_indices[-1])
        n_frames = len(self._data.time)
        if last_contact >= n_frames - 10:
            # Projectile still in contact at end of simulation → sticking regime
            warnings.warn(
                "Projectile still in contact at end of simulation "
                "(sticking / COR ≈ 0). Using final trajectory frames.",
                stacklevel=3,
            )
            return max(0, n_frames - 50)
        return last_contact + 1

    @property
    @require_data
    @cached
    @validate_physical("contact_time", 0.0, float("inf"))
    def contact_time(self) -> float:
        """Duration of contact phase [tau]."""
        t = self._data.time
        i_start = self._approach_end_idx()
        i_end = self._rebound_start_idx()
        return float(t[i_end] - t[i_start])

    @property
    @require_data
    @cached
    def final_velocity_z(self) -> float:
        """Z-velocity of projectile after clear rebound [sigma/tau] (positive = up)."""
        i_end = self._rebound_start_idx()
        # Average over a window of 50 frames after rebound to reduce noise
        window = self._data.vz[i_end : i_end + 50]
        if len(window) < 10:
            warnings.warn(
                "Less than 10 post-rebound frames available; COR may be noisy.",
                stacklevel=2,
            )
        vf = float(np.mean(window))
        logger.debug("v_final_z = %.4f (mean over %d frames)", vf, len(window))
        return vf

    @property
    @require_data
    @cached
    @validate_physical("COR", 0.0, 1.05)   # allow tiny super-elasticity from noise
    def coefficient_of_restitution(self) -> float:
        """Normal coefficient of restitution: e = |vf_z| / v0_z."""
        vf = self.final_velocity_z
        v0z = self.v0_z
        if v0z == 0.0:
            raise ZeroDivisionError("Incident normal velocity is zero; COR undefined.")
        return abs(vf) / v0z

    @property
    @require_data
    @cached
    def energy_dissipated(self) -> float:
        """Energy dissipated in the collision [eps]: ΔE = ½m(v0²-vf²)."""
        v0z = self.v0_z
        vf = self.final_velocity_z
        return 0.5 * (v0z ** 2 - vf ** 2)  # mass = 1 in LJ units

    @property
    @require_data
    @cached
    def max_penetration(self) -> float:
        """Maximum penetration below the unperturbed surface top [sigma]."""
        z_min = float(np.min(self._data.z))
        return max(0.0, self.surface_top_z - z_min)

    @property
    @require_data
    @cached
    def surface_temperature_rise(self) -> float:
        """Rise in surface temperature during impact [eps/kB]."""
        # Average T_surf over the last 200 frames (post-rebound)
        T_final = float(np.mean(self._data.T_surf[-200:]))
        return T_final - self.T_surface_init

    # --- Energy accommodation (generalised) -----------------------------------

    @property
    @require_data
    @cached
    def incident_energy(self) -> float:
        """Total incident kinetic energy E* = ½v₀² [eps].

        Measured as mean ke_proj over the pre-contact approach frames,
        providing E* directly from the trajectory rather than assuming v0.
        """
        i_end = self._approach_end_idx()
        # Use last 200 approach frames to avoid startup transients
        start = max(0, i_end - 200)
        return float(np.mean(self._data.ke_proj[start:i_end]))

    @property
    @require_data
    @cached
    def final_energy(self) -> float:
        """Total final kinetic energy E_f [eps] after rebound.

        Measured as mean ke_proj over the first 50 post-rebound frames.
        """
        i_start = self._rebound_start_idx()
        window = self._data.ke_proj[i_start: i_start + 50]
        return float(np.mean(window))

    @property
    @require_data
    @cached
    @validate_physical("EAC", -1.05, 1.05)   # negative EAC is physical in sub-thermal regime
    def energy_accommodation(self) -> float:
        """Energy Accommodation Coefficient α_E = (E_inc - E_fin) / E_inc.

        General observable for Gonzalo's study: measures the fraction of the
        total incident kinetic energy transferred to the surface (phonons + heat).
        Independent of impact angle — valid for all (v0, T*, θ) combinations.

        Range: [0, 1]
          0 → elastic (no energy transfer)
          1 → full accommodation (projectile sticks)
        """
        E_inc = self.incident_energy
        if E_inc < 1e-12:
            raise ZeroDivisionError("Incident energy is essentially zero; EAC undefined.")
        E_fin = self.final_energy
        return float((E_inc - E_fin) / E_inc)

    @property
    @require_data
    @cached
    def eac_normal(self) -> float:
        """Normal EAC α_n = 1 - COR² = ΔE_n / E_inc_n.

        Fraction of normal kinetic energy transferred during the collision.
        For normal incidence (angle=0), α_n = α_E.
        """
        return float(1.0 - self.coefficient_of_restitution ** 2)

    @property
    @require_data
    @cached
    def final_velocity_x(self) -> float:
        """Tangential (x) velocity of projectile after rebound [sigma/tau]."""
        i_start = self._rebound_start_idx()
        window = self._data.vx[i_start: i_start + 50]
        return float(np.mean(window))

    @property
    @require_data
    @cached
    def eac_tangential(self) -> float:
        """Tangential EAC α_t = 1 - (vf_x / v0_x)².

        Fraction of tangential kinetic energy transferred.
        Returns 0.0 for normal incidence (angle ≈ 0, v0_x ≈ 0) where α_t
        is undefined (no tangential energy to transfer).
        """
        import math
        v0_x = self.v0 * math.sin(math.radians(self.angle))
        if abs(v0_x) < 1e-10:
            return 0.0
        vf_x = self.final_velocity_x
        alpha_t = 1.0 - (vf_x / v0_x) ** 2
        # Clip to physical range; noise can push it slightly outside [0, 1]
        return float(np.clip(alpha_t, 0.0, 1.0))

    # --- Summary --------------------------------------------------------------

    @require_data
    @timed
    def compute_all(self) -> ReboundResult:
        """Compute and return all physical observables as a ``ReboundResult``."""
        return ReboundResult(
            v0_z=self.v0_z,
            vf_z=self.final_velocity_z,
            COR=self.coefficient_of_restitution,
            delta_E=self.energy_dissipated,
            contact_time=self.contact_time,
            max_penetration=self.max_penetration,
            dT_surface=self.surface_temperature_rise,
            angle=self.angle,
            T_surface=self.T_surface_init,
            EAC=self.energy_accommodation,
            EAC_normal=self.eac_normal,
            EAC_tangential=self.eac_tangential,
            E_star=self.incident_energy,
        )

    def __repr__(self) -> str:
        loaded = self._data is not None
        return (
            f"ReboundAnalyzer(v0={self.v0}, angle={self.angle}, "
            f"data_loaded={loaded})"
        )


# ==============================================================================
# Parametric analysis
# ==============================================================================

class ParametricAnalyzer:
    """Aggregates results from multiple ImpactSimulation runs.

    Collects ``ReboundResult`` objects and provides numpy-backed aggregate
    statistics for COR vs. v0 curves, temperature dependencies, etc.

    Parameters
    ----------
    results : list[ReboundResult]
        Individual impact results (can be empty initially).
    """

    def __init__(self, results: Optional[list[ReboundResult]] = None) -> None:
        self.results: list[ReboundResult] = results or []

    def add(self, result: ReboundResult) -> None:
        self.results.append(result)

    def __len__(self) -> int:
        return len(self.results)

    def _require_results(self, n: int = 2) -> None:
        if len(self.results) < n:
            raise InsufficientDataError("parametric analysis", len(self.results), n)

    @property
    def v0_array(self) -> np.ndarray:
        self._require_results()
        return np.array([r.v0_z for r in self.results])

    @property
    def COR_array(self) -> np.ndarray:
        self._require_results()
        return np.array([r.COR for r in self.results])

    @property
    def delta_E_array(self) -> np.ndarray:
        self._require_results()
        return np.array([r.delta_E for r in self.results])

    @property
    def EAC_array(self) -> np.ndarray:
        """Energy Accommodation Coefficients α_E for all results."""
        self._require_results()
        return np.array([r.EAC for r in self.results])

    @property
    def E_star_array(self) -> np.ndarray:
        """Reduced incident energies E* = ½v₀² for all results [eps]."""
        self._require_results()
        return np.array([r.E_star for r in self.results])

    @property
    def T_surface_array(self) -> np.ndarray:
        """Surface temperatures T* for all results [eps/kB]."""
        self._require_results()
        return np.array([r.T_surface for r in self.results])

    @timed
    def cor_vs_v0(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Group results by v0, return (v0_unique, COR_mean, COR_std).

        Returns
        -------
        v0_unique : 1-D ndarray
            Sorted unique incident velocities.
        COR_mean : 1-D ndarray
            Mean COR at each velocity (across replicas).
        COR_std : 1-D ndarray
            Standard deviation of COR at each velocity.
        """
        self._require_results()
        v0s = self.v0_array
        cors = self.COR_array
        unique_v0 = np.unique(v0s)
        means = np.array([cors[v0s == v].mean() for v in unique_v0])
        stds  = np.array([cors[v0s == v].std()  for v in unique_v0])
        return unique_v0, means, stds

    @timed
    def eac_vs_Estar(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Group results by E*, return (E*_unique, EAC_mean, EAC_std).

        Primary curve for Gonzalo's generalisation study: α_E(E*) collapses
        data from different v0 onto a single reduced-variable axis.

        Returns
        -------
        E_star_unique : 1-D ndarray
            Sorted unique reduced incident energies E* = ½v₀² [eps].
        EAC_mean : 1-D ndarray
            Mean α_E at each E* (across replicas and angles if present).
        EAC_std : 1-D ndarray
            Standard deviation of α_E at each E*.
        """
        self._require_results()
        E_stars = self.E_star_array
        eacs = self.EAC_array
        unique_E = np.unique(np.round(E_stars, decimals=6))
        means = np.array([eacs[np.abs(E_stars - E) < 1e-5].mean() for E in unique_E])
        stds  = np.array([eacs[np.abs(E_stars - E) < 1e-5].std()  for E in unique_E])
        return unique_E, means, stds

    @timed
    def eac_vs_T(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Group results by surface temperature T*, return (T*_unique, EAC_mean, EAC_std).

        Secondary curve for generalisation: α_E(T*) at fixed E* shows how
        thermal surface agitation affects energy accommodation.

        Returns
        -------
        T_unique : 1-D ndarray
            Sorted unique surface temperatures [eps/kB].
        EAC_mean : 1-D ndarray
            Mean α_E at each T*.
        EAC_std : 1-D ndarray
            Standard deviation of α_E at each T*.
        """
        self._require_results()
        Ts = self.T_surface_array
        eacs = self.EAC_array
        unique_T = np.unique(np.round(Ts, decimals=6))
        means = np.array([eacs[np.abs(Ts - T) < 1e-5].mean() for T in unique_T])
        stds  = np.array([eacs[np.abs(Ts - T) < 1e-5].std()  for T in unique_T])
        return unique_T, means, stds

    @timed
    def fit_eac_power_law(self) -> tuple[float, float]:
        """Fit α_E ~ A * E*^β to the data via log-log regression.

        Returns
        -------
        A : float
            Prefactor.
        beta : float
            Power-law exponent (positive → α_E increases with E*).

        Notes
        -----
        Linearises via ln(α_E) = ln(A) + β·ln(E*).
        Requires >= 4 distinct E* values.
        """
        self._require_results(4)
        E_stars, eacs, _ = self.eac_vs_Estar()
        mask = (E_stars > 0) & (eacs > 0) & (eacs < 1)
        if mask.sum() < 4:
            raise InsufficientDataError("EAC power-law fit", int(mask.sum()), 4)
        ln_E   = np.log(E_stars[mask])
        ln_eac = np.log(eacs[mask])
        A_mat = np.column_stack([np.ones_like(ln_E), ln_E])
        coeffs, *_ = np.linalg.lstsq(A_mat, ln_eac, rcond=None)
        ln_A, beta = coeffs
        return float(np.exp(ln_A)), float(beta)

    @timed
    def fit_cor_power_law(self) -> tuple[float, float]:
        """Fit COR ~ A * v0^(-beta) to the data.

        Returns
        -------
        A : float
            Prefactor.
        beta : float
            Power-law exponent (positive for decreasing COR with v0).

        Notes
        -----
        Linearises via log-log regression:  ln(COR) = ln(A) - beta*ln(v0).
        Requires >= 4 data points.
        """
        self._require_results(4)
        v0s, cors, _ = self.cor_vs_v0()
        # Remove zeros/negatives to avoid log issues
        mask = (v0s > 0) & (cors > 0)
        if mask.sum() < 4:
            raise InsufficientDataError("power-law fit", int(mask.sum()), 4)
        ln_v0  = np.log(v0s[mask])
        ln_cor = np.log(cors[mask])
        # Least-squares: ln_cor = ln_A - beta * ln_v0
        A_mat = np.column_stack([np.ones_like(ln_v0), ln_v0])
        coeffs, *_ = np.linalg.lstsq(A_mat, ln_cor, rcond=None)
        ln_A, neg_beta = coeffs
        return float(np.exp(ln_A)), float(-neg_beta)

    def summary(self) -> str:
        """Return a human-readable summary string."""
        lines = [
            f"ParametricAnalyzer: {len(self.results)} results",
            f"  v0 range:    {self.v0_array.min():.3f} – {self.v0_array.max():.3f}  σ/τ",
            f"  E* range:    {self.E_star_array.min():.4f} – {self.E_star_array.max():.4f}  ε",
            f"  COR range:   {self.COR_array.min():.4f} – {self.COR_array.max():.4f}",
            f"  α_E range:   {self.EAC_array.min():.4f} – {self.EAC_array.max():.4f}",
            f"  mean α_E:    {self.EAC_array.mean():.4f} ± {self.EAC_array.std():.4f}",
        ]
        unique_T = np.unique(self.T_surface_array)
        if len(unique_T) > 1:
            lines.append(f"  T* values:   {', '.join(f'{t:.3f}' for t in unique_T)}  ε/k_B")
        return "\n".join(lines)

    def __repr__(self) -> str:
        return f"ParametricAnalyzer(n_results={len(self.results)})"
