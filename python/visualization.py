"""
visualization.py
================
Matplotlib-based visualisation for inelastic rebound simulation data.

All plot functions accept a ``TrajectoryData`` or aggregated data and return
``matplotlib.figure.Figure`` objects. Figures are never shown automatically
(no ``plt.show()`` calls); the caller decides whether to display or save.

Advanced patterns:
- ``PlotStyle`` context manager: temporarily switches rcParams.
- ``FigureBuilder``: fluent builder pattern for complex multi-panel figures.
- ``@save_figure`` decorator: auto-saves if an output path is given.
"""

from __future__ import annotations

import functools
import logging
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Callable, Generator, Optional, Sequence

import numpy as np

logger = logging.getLogger(__name__)

# Lazy import matplotlib so the module can be imported without a display
def _mpl():
    try:
        import matplotlib.pyplot as plt
        import matplotlib
        return plt, matplotlib
    except ImportError as exc:
        raise ImportError(
            "matplotlib is required for visualisation. "
            "Install it with: pip install matplotlib"
        ) from exc


# ==============================================================================
# Decorator: auto-save
# ==============================================================================

def save_figure(func: Callable) -> Callable:
    """Decorator: if ``output_path`` kwarg is given, save the figure there.

    The wrapped function must return a ``matplotlib.figure.Figure``.
    """
    @functools.wraps(func)
    def wrapper(*args, output_path: Optional[Path | str] = None, **kwargs):
        fig = func(*args, **kwargs)
        if output_path is not None:
            path = Path(output_path)
            path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(path, dpi=150, bbox_inches="tight")
            logger.info("Figure saved to %s", path)
        return fig
    return wrapper


# ==============================================================================
# Plot style context manager
# ==============================================================================

@contextmanager
def publication_style() -> Generator[None, None, None]:
    """Context manager: apply publication-quality rcParams within the block.

    Usage
    -----
    >>> with publication_style():
    ...     fig = plot_trajectory(data)
    """
    plt, mpl = _mpl()
    params = {
        "font.family":       "serif",
        "font.size":         11,
        "axes.labelsize":    12,
        "axes.titlesize":    13,
        "xtick.labelsize":   10,
        "ytick.labelsize":   10,
        "legend.fontsize":   10,
        "figure.dpi":        150,
        "lines.linewidth":   1.5,
        "axes.grid":         True,
        "grid.alpha":        0.35,
        "savefig.bbox":      "tight",
    }
    with mpl.rc_context(rc=params):
        yield


# ==============================================================================
# Individual plot functions
# ==============================================================================

@save_figure
def plot_trajectory(
    data: Any,   # TrajectoryData (avoid circular import)
    surface_top_z: float = 9.3,
    lj_cutoff: float = 2.5,
) -> Any:  # matplotlib.figure.Figure
    """Four-panel trajectory overview.

    Panels:
    1. z(t) – projectile height vs. time
    2. vz(t) – normal velocity vs. time
    3. KE(t) – projectile kinetic energy vs. time
    4. T_surf(t) – surface temperature vs. time

    Contact region is shaded on all panels.
    """
    plt, _ = _mpl()

    with publication_style():
        fig, axes = plt.subplots(4, 1, figsize=(8, 10), sharex=True)

        contact_mask = data.d_contact < lj_cutoff
        t_contact_start = data.time[contact_mask][0] if contact_mask.any() else None
        t_contact_end   = data.time[contact_mask][-1] if contact_mask.any() else None

        def shade_contact(ax):
            if t_contact_start is not None:
                ax.axvspan(t_contact_start, t_contact_end,
                           color="orange", alpha=0.15, label="contact")

        ax1, ax2, ax3, ax4 = axes

        # Panel 1: z(t)
        ax1.plot(data.time, data.z, color="steelblue")
        ax1.axhline(surface_top_z, color="brown", linestyle="--", lw=1.0,
                    label=f"surface top (z={surface_top_z})")
        shade_contact(ax1)
        ax1.set_ylabel("z  [σ]")
        ax1.set_title("Projectile Trajectory – Inelastic Rebound MD")
        ax1.legend(loc="upper right", fontsize=9)

        # Panel 2: vz(t)
        ax2.plot(data.time, data.vz, color="darkorange")
        ax2.axhline(0, color="gray", linestyle=":", lw=0.8)
        shade_contact(ax2)
        ax2.set_ylabel("v_z  [σ/τ]")

        # Panel 3: KE(t)
        ax3.plot(data.time, data.ke_proj, color="seagreen")
        shade_contact(ax3)
        ax3.set_ylabel("KE_proj  [ε]")

        # Panel 4: T_surf(t)
        ax4.plot(data.time, data.T_surf, color="firebrick")
        shade_contact(ax4)
        ax4.set_ylabel("T_surf  [ε/k_B]")
        ax4.set_xlabel("time  [τ]")

        for ax in axes:
            ax.legend(loc="best", fontsize=8)

        fig.tight_layout()

    return fig


@save_figure
def plot_cor_vs_v0(
    v0_values: np.ndarray,
    cor_mean: np.ndarray,
    cor_std: Optional[np.ndarray] = None,
    fit_params: Optional[tuple[float, float]] = None,
) -> Any:
    """COR vs. incident velocity plot with optional errorbars and power-law fit.

    Parameters
    ----------
    v0_values : 1-D array
        Unique incident velocities [sigma/tau].
    cor_mean : 1-D array
        Mean COR at each v0.
    cor_std : 1-D array | None
        Standard deviation of COR (error bars).
    fit_params : (A, beta) | None
        Power-law fit parameters; if given, overlay the fit curve.
    """
    plt, _ = _mpl()

    with publication_style():
        fig, ax = plt.subplots(figsize=(6, 5))

        if cor_std is not None:
            ax.errorbar(v0_values, cor_mean, yerr=cor_std,
                        fmt="o", color="steelblue", capsize=4,
                        label="COR (mean ± std)")
        else:
            ax.plot(v0_values, cor_mean, "o-", color="steelblue", label="COR")

        if fit_params is not None:
            A, beta = fit_params
            v0_fit = np.linspace(v0_values.min(), v0_values.max(), 200)
            cor_fit = A * v0_fit ** (-beta)
            ax.plot(v0_fit, cor_fit, "--", color="darkorange",
                    label=f"fit: {A:.3f}·v₀^(-{beta:.3f})")

        ax.axhline(1.0, color="gray", linestyle=":", lw=0.8, label="elastic (e=1)")
        ax.set_xlabel("Incident velocity  v₀  [σ/τ]")
        ax.set_ylabel("Coefficient of Restitution  e")
        ax.set_title("COR vs. Incident Velocity")
        ax.set_ylim(0, 1.1)
        ax.legend()

        fig.tight_layout()

    return fig


@save_figure
def plot_eac_vs_Estar(
    E_star_values: np.ndarray,
    eac_mean: np.ndarray,
    eac_std: Optional[np.ndarray] = None,
    fit_params: Optional[tuple[float, float]] = None,
    T_label: str = "",
    xlabel: str = r"Reduced incident energy  E* = ½v₀²  [ε]",
    title: str = r"Energy Accommodation Coefficient vs. E*",
) -> Any:
    """α_E vs. reduced incident energy E* — Gonzalo's general observable.

    Parameters
    ----------
    E_star_values : 1-D array
        Unique reduced incident energies E* = ½v₀² [eps].
    eac_mean : 1-D array
        Mean α_E at each E*.
    eac_std : 1-D array | None
        Standard deviation of α_E (error bars).
    fit_params : (A, beta) | None
        Power-law fit α_E ~ A * E*^beta; if given, overlay the fit curve.
    T_label : str
        Legend label suffix, e.g. "T*=0.1" to identify a temperature series.
    xlabel : str
        X-axis label (override for α_E vs T* use case).
    title : str
        Figure title.
    """
    plt, _ = _mpl()

    with publication_style():
        fig, ax = plt.subplots(figsize=(6, 5))

        label = f"α_E" + (f"  ({T_label})" if T_label else "")
        if eac_std is not None:
            ax.errorbar(E_star_values, eac_mean, yerr=eac_std,
                        fmt="o", color="darkorange", capsize=4, label=label)
        else:
            ax.plot(E_star_values, eac_mean, "o-", color="darkorange", label=label)

        if fit_params is not None:
            A, beta = fit_params
            E_fit = np.linspace(E_star_values.min(), E_star_values.max(), 300)
            eac_fit = A * E_fit ** beta
            ax.plot(E_fit, eac_fit, "--", color="steelblue",
                    label=rf"fit: {A:.3f}·E*$^{{{beta:+.3f}}}$")

        ax.axhline(1.0, color="gray", linestyle=":", lw=0.8, label="full accommodation")
        ax.axhline(0.0, color="gray", linestyle=":", lw=0.8)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r"Energy Accommodation Coefficient  α$_E$")
        ax.set_title(title)
        ax.set_ylim(-0.05, 1.1)
        ax.legend()

        fig.tight_layout()

    return fig


@save_figure
def plot_energy_partition(
    data: Any,   # TrajectoryData
    v0: float,
) -> Any:
    """Energy partition: projectile KE + transferred (surface) energy vs. time.

    Shows conservation of energy breakdown: KE_proj + ΔE_surf ≈ KE_initial.
    """
    plt, _ = _mpl()

    with publication_style():
        fig, ax = plt.subplots(figsize=(7, 4))

        KE_0 = float(data.ke_proj[0])
        # Approximate transferred energy as rise in surface thermal energy
        E_surf_rel = (data.T_surf - data.T_surf[0]) * 3.0  # 3/2 * NkB * dT (approx)

        ax.plot(data.time, data.ke_proj / KE_0, color="steelblue",
                label="KE_proj / KE₀")
        ax.plot(data.time, E_surf_rel / KE_0, color="firebrick", linestyle="--",
                label="ΔE_surf / KE₀  (approx)")
        ax.axhline(1.0, color="gray", linestyle=":", lw=0.8, label="initial KE")

        ax.set_xlabel("time  [τ]")
        ax.set_ylabel("Normalised energy")
        ax.set_title("Energy Partition During Impact")
        ax.legend()

        fig.tight_layout()

    return fig


@save_figure
def plot_parametric_heatmap(
    v0_values: Sequence[float],
    T_values: Sequence[float],
    COR_matrix: np.ndarray,
) -> Any:
    """2-D heatmap of COR as a function of v0 and surface temperature.

    Parameters
    ----------
    v0_values : sequence of float
        Incident velocities (x-axis).
    T_values : sequence of float
        Surface temperatures (y-axis).
    COR_matrix : 2-D ndarray, shape (len(T_values), len(v0_values))
        COR values.
    """
    plt, _ = _mpl()

    with publication_style():
        fig, ax = plt.subplots(figsize=(7, 5))

        im = ax.imshow(
            COR_matrix,
            origin="lower",
            aspect="auto",
            extent=[
                min(v0_values), max(v0_values),
                min(T_values),  max(T_values),
            ],
            cmap="plasma",
            vmin=0.0,
            vmax=1.0,
        )
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label("COR  e")

        ax.set_xlabel("Incident velocity  v₀  [σ/τ]")
        ax.set_ylabel("Surface temperature  T*  [ε/k_B]")
        ax.set_title("COR Heatmap: v₀ vs T_surface")

        fig.tight_layout()

    return fig


# ==============================================================================
# Fluent FigureBuilder
# ==============================================================================

class FigureBuilder:
    """Fluent builder for constructing multi-panel figures incrementally.

    Usage
    -----
    >>> fig = (
    ...     FigureBuilder(nrows=2, ncols=2, figsize=(10, 8))
    ...     .title("Impact Summary")
    ...     .panel(0, 0, data.time, data.z, ylabel="z [σ]", color="blue")
    ...     .panel(0, 1, data.time, data.vz, ylabel="vz [σ/τ]")
    ...     .panel(1, 0, v0_arr, cor_arr, ylabel="COR", xlabel="v0 [σ/τ]")
    ...     .build()
    ... )
    """

    def __init__(self, nrows: int = 2, ncols: int = 2, figsize: tuple = (10, 8)) -> None:
        plt, _ = _mpl()
        self._plt = plt
        self._fig, self._axes = plt.subplots(nrows, ncols, figsize=figsize)
        if nrows == 1 and ncols == 1:
            self._axes = np.array([[self._axes]])
        elif nrows == 1 or ncols == 1:
            self._axes = self._axes.reshape(nrows, ncols)
        self._nrows = nrows
        self._ncols = ncols

    def title(self, text: str) -> "FigureBuilder":
        self._fig.suptitle(text, fontsize=14)
        return self

    def panel(
        self,
        row: int,
        col: int,
        x: np.ndarray,
        y: np.ndarray,
        *,
        xlabel: str = "",
        ylabel: str = "",
        label: str = "",
        color: str = "steelblue",
        linestyle: str = "-",
    ) -> "FigureBuilder":
        ax = self._axes[row, col]
        ax.plot(x, y, color=color, linestyle=linestyle, label=label or None)
        if xlabel:
            ax.set_xlabel(xlabel)
        if ylabel:
            ax.set_ylabel(ylabel)
        if label:
            ax.legend()
        return self

    def build(self) -> Any:  # matplotlib.figure.Figure
        self._fig.tight_layout()
        return self._fig
