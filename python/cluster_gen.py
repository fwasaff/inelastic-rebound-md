"""
cluster_gen.py
==============
Generate approximate initial geometries for Lennard-Jones clusters.

These are starting points for LAMMPS energy minimisation (in.minimize_cluster).
Supported sizes: 1, 2, 7, 13, 55, 147  (Mackay icosahedral series + small).

Physics
-------
The global minima for LJ clusters follow the Mackay icosahedral series:
  N = 1  : monomer
  N = 2  : dimer
  N = 7  : pentagonal bipyramid  (D_5h)
  N = 13 : Mackay icosahedron, 1 shell  (I_h)
  N = 55 : Mackay icosahedron, 2 shells (I_h)
  N = 147: Mackay icosahedron, 3 shells (I_h)

Shell k of a Mackay icosahedron has:
  12 vertex atoms  +  30(k-1) edge atoms  +  20·T(k-2) face atoms
  where T(n) = n(n-1)/2  →  total shell k = 10k² + 2

Reference: Mackay, Acta Cryst. 15, 916 (1962)

Usage
-----
    from python.cluster_gen import generate_all_clusters
    generate_all_clusters("lammps/clusters", sizes=[1, 2, 7, 13, 55])
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)

# LJ pair equilibrium distance (2^(1/6) sigma)
_R_EQ: float = 2.0 ** (1.0 / 6.0)

# Golden ratio
_PHI: float = (1.0 + np.sqrt(5.0)) / 2.0

# Standard icosahedron has vertices at cyclic permutations of (0, ±1, ±φ)
# with edge length = 2.  We scale so that edge length = _R_EQ.
_ICO_EDGE_SCALE: float = _R_EQ / 2.0

# Cached icosahedron data (computed once)
_ICO_VERTS: np.ndarray | None = None
_ICO_EDGES: list[tuple[int, int]] | None = None
_ICO_FACES: list[tuple[int, int, int]] | None = None


# ==============================================================================
# Icosahedron geometry helpers
# ==============================================================================

def _get_ico_verts() -> np.ndarray:
    """12 vertices of a regular icosahedron (edge length = _R_EQ)."""
    global _ICO_VERTS
    if _ICO_VERTS is None:
        raw = []
        for s1 in (1, -1):
            for s2 in (1, -1):
                raw += [
                    [0.0,        float(s1),        float(s2) * _PHI],
                    [float(s1),  float(s2) * _PHI, 0.0             ],
                    [float(s2) * _PHI, 0.0,        float(s1)       ],
                ]
        _ICO_VERTS = np.array(raw) * _ICO_EDGE_SCALE
    return _ICO_VERTS


def _get_ico_edges() -> list[tuple[int, int]]:
    """30 nearest-neighbor edges of the icosahedron (by index into _get_ico_verts)."""
    global _ICO_EDGES
    if _ICO_EDGES is None:
        verts = _get_ico_verts()
        nn_dist = 2.0 * _ICO_EDGE_SCALE          # nearest-neighbor distance = _R_EQ
        threshold = nn_dist * 1.05               # 5 % tolerance
        edges: list[tuple[int, int]] = []
        for i in range(12):
            for j in range(i + 1, 12):
                if np.linalg.norm(verts[i] - verts[j]) < threshold:
                    edges.append((i, j))
        assert len(edges) == 30, f"Expected 30 icosahedron edges, got {len(edges)}"
        _ICO_EDGES = edges
    return _ICO_EDGES


def _get_ico_faces() -> list[tuple[int, int, int]]:
    """20 triangular faces of the icosahedron (vertex index triples)."""
    global _ICO_FACES
    if _ICO_FACES is None:
        edge_set = {frozenset(e) for e in _get_ico_edges()}
        faces: list[tuple[int, int, int]] = []
        for i in range(12):
            for j in range(i + 1, 12):
                for k in range(j + 1, 12):
                    if (frozenset({i, j}) in edge_set and
                            frozenset({j, k}) in edge_set and
                            frozenset({i, k}) in edge_set):
                        faces.append((i, j, k))
        assert len(faces) == 20, f"Expected 20 icosahedron faces, got {len(faces)}"
        _ICO_FACES = faces
    return _ICO_FACES


def _icosahedron_shell(k: int) -> np.ndarray:
    """Generate the k-th shell of a Mackay icosahedron.

    Shell k contains:
      - 12 vertex atoms  (at scale k)
      - 30(k-1) edge atoms  ((k-1) atoms per each of 30 edges)
      - 20 · T(k-2) face atoms  (T(n) = n(n-1)/2 interior atoms per face)

    Total: 10k² + 2 atoms.

    All nearest-neighbor distances are approximately _R_EQ (LJ sigma units).
    """
    scale = _ICO_EDGE_SCALE * k
    verts = _get_ico_verts() / _ICO_EDGE_SCALE * scale   # scale shell-k vertices
    atoms: list[np.ndarray] = list(verts)

    # Edge atoms: (k-1) equally-spaced atoms along each edge
    for (i, j) in _get_ico_edges():
        v1, v2 = verts[i], verts[j]
        for m in range(1, k):
            t = float(m) / float(k)
            atoms.append((1.0 - t) * v1 + t * v2)

    # Face atoms: strictly-interior points of each triangular face.
    # For side-k subdivisions: interior count = (k-2)(k-1)/2 per face.
    # Barycentric coords (p, q, r) with p+q+r=k, p,q,r >= 1.
    for (fi, fj, fk) in _get_ico_faces():
        v0, v1, v2 = verts[fi], verts[fj], verts[fk]
        for p in range(1, k - 1):       # p: 1 … k-2
            for q in range(1, k - p):   # q: 1 … k-p-1  (ensures r=k-p-q >= 1)
                r = k - p - q
                atoms.append((p * v0 + q * v1 + r * v2) / float(k))

    result = np.array(atoms)
    expected = 10 * k * k + 2
    assert len(result) == expected, (
        f"Shell k={k}: expected {expected} atoms, got {len(result)}"
    )
    return result


# ==============================================================================
# Cluster geometry generators
# ==============================================================================

def generate_cluster(N: int) -> np.ndarray:
    """Return approximate initial positions for an N-atom LJ cluster.

    Positions are centred at the origin.  Pass through LAMMPS minimisation
    (in.minimize_cluster) to obtain the true local minimum.

    Parameters
    ----------
    N : int
        Number of atoms.  Supported values: 1, 2, 7, 13, 55, 147.

    Returns
    -------
    positions : ndarray, shape (N, 3)
        Atom positions in LJ sigma units.

    Raises
    ------
    ValueError
        If N is not in the supported set.
    """
    supported = {1, 2, 7, 13, 55, 147}
    if N not in supported:
        raise ValueError(
            f"Unsupported cluster size N={N}.  Supported: {sorted(supported)}"
        )

    if N == 1:
        return np.zeros((1, 3))

    if N == 2:
        return np.array([[0.0, 0.0, -_R_EQ / 2.0],
                         [0.0, 0.0,  _R_EQ / 2.0]])

    if N == 7:
        # Pentagonal bipyramid: 5 equatorial + 2 apical
        # Ring radius chosen so equatorial–apical distance ≈ _R_EQ
        r_ring = 0.951 * _R_EQ
        h = float(np.sqrt(max(0.0, _R_EQ ** 2 - r_ring ** 2)))
        equatorial = np.array([
            [r_ring * np.cos(2 * np.pi * m / 5),
             r_ring * np.sin(2 * np.pi * m / 5),
             0.0]
            for m in range(5)
        ])
        apical = np.array([[0.0, 0.0, h], [0.0, 0.0, -h]])
        return np.vstack([equatorial, apical])

    if N == 13:
        center = np.zeros((1, 3))
        return np.vstack([center, _icosahedron_shell(1)])

    if N == 55:
        core = generate_cluster(13)
        return np.vstack([core, _icosahedron_shell(2)])

    if N == 147:
        core = generate_cluster(55)
        return np.vstack([core, _icosahedron_shell(3)])

    raise ValueError(f"N={N} not handled")   # unreachable


def cluster_radius(N: int) -> float:
    """Approximate outer radius of the N-atom LJ cluster [sigma].

    Used to set z_gap so that the cluster starts above the surface.
    """
    positions = generate_cluster(N)
    return float(np.max(np.linalg.norm(positions, axis=1))) + 0.5


# ==============================================================================
# LAMMPS data file writer
# ==============================================================================

def write_lammps_cluster_data(
    positions: np.ndarray,
    output_path: Path | str,
    atom_type: int = 2,
    padding: float = 6.0,
) -> None:
    """Write cluster atom positions to a LAMMPS data file.

    Atoms are written as type ``atom_type`` (default 2 = projectile type).
    The file declares 2 atom types so it is compatible with ``read_data add merge``
    into a surface restart that already has types 1 and 2.

    Parameters
    ----------
    positions : ndarray, shape (N, 3)
        Cluster atom positions.  COM is recentred at origin before writing.
    output_path : Path | str
        Where to write the ``.dat`` file.
    atom_type : int
        Atom type index written for every atom (default: 2).
    padding : float
        Vacuum padding added around the cluster bounding box [sigma].
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Re-centre at COM
    com = positions.mean(axis=0)
    pos = positions - com

    N = len(pos)
    xlo = pos[:, 0].min() - padding
    xhi = pos[:, 0].max() + padding
    ylo = pos[:, 1].min() - padding
    yhi = pos[:, 1].max() + padding
    zlo = pos[:, 2].min() - padding
    zhi = pos[:, 2].max() + padding

    lines = [
        f"# LJ cluster N={N}, centered at origin, atom type={atom_type}",
        "",
        f"{N} atoms",
        "2 atom types",
        "",
        f"{xlo:.8f} {xhi:.8f} xlo xhi",
        f"{ylo:.8f} {yhi:.8f} ylo yhi",
        f"{zlo:.8f} {zhi:.8f} zlo zhi",
        "",
        "Masses",
        "",
        "1 1.0",
        f"2 1.0",
        "",
        "Atoms  # atomic",
        "",
    ]
    for atom_id, (x, y, z) in enumerate(pos, start=1):
        lines.append(f"{atom_id} {atom_type} {x:.10f} {y:.10f} {z:.10f}")

    output_path.write_text("\n".join(lines) + "\n")
    logger.debug("Wrote %d-atom cluster to %s", N, output_path)


# ==============================================================================
# High-level: generate all cluster data files
# ==============================================================================

def generate_all_clusters(
    output_dir: Path | str = "lammps/clusters",
    sizes: list[int] | None = None,
) -> dict[int, Path]:
    """Generate initial LAMMPS data files for all requested cluster sizes.

    Parameters
    ----------
    output_dir : Path | str
        Directory where files are written (created if absent).
    sizes : list[int] | None
        Cluster sizes to generate.  Defaults to [1, 2, 7, 13, 55].

    Returns
    -------
    paths : dict[int, Path]
        Mapping from N to the written file path.
    """
    if sizes is None:
        sizes = [1, 2, 7, 13, 55]

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    paths: dict[int, Path] = {}
    for N in sizes:
        positions = generate_cluster(N)
        out = output_dir / f"cluster_N{N}_init.dat"
        write_lammps_cluster_data(positions, out)
        logger.info("Generated N=%d cluster → %s", N, out)
        paths[N] = out

    return paths
