"""Tetrahedral meshes (stage S4b): mesh helpers and the facade to the C++ transport.

The transport itself is C++ (`cpp/include/scattuccino/photon_tet.hpp`); it is checked against the
voxel walk on the same geometry (a voxel grid split into tetrahedra) and against geometric optics
for a curved surface (glass sphere).

- `mesh_from_voxels(labels, voxel)`: every voxel split into 6 tetrahedra (Kuhn triangulation,
  conforming across voxels), same coordinates as `scattuccino.voxels`;
- `icosphere(radius, level)`: a sphere as a fan of tetrahedra from the centre;
- `simulate_tets(...)`: runs the C++ transport and returns a `Protocol`.
"""

from __future__ import annotations

import itertools
import os
from dataclasses import dataclass

import numpy as np

from scattuccino._types import FloatArray, IntArray
from scattuccino.layers import Layer
from scattuccino.transport import Protocol


@dataclass
class TetMesh:
    nodes: FloatArray  # (N, 3) mm
    tets: IntArray  # (M, 4)
    labels: IntArray  # (M,)


def mesh_from_voxels(labels: IntArray, voxel: tuple[float, float, float]) -> TetMesh:
    lab = np.asarray(labels, dtype=np.int64)
    nx, ny, nz = lab.shape
    dx, dy, dz = voxel
    xs = (np.arange(nx + 1) - nx / 2) * dx
    ys = (np.arange(ny + 1) - ny / 2) * dy
    zs = np.arange(nz + 1) * dz
    gx, gy, gz = np.meshgrid(xs, ys, zs, indexing="ij")
    nodes = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])

    def nid(i: IntArray, j: IntArray, k: IntArray) -> IntArray:
        return np.asarray((i * (ny + 1) + j) * (nz + 1) + k)

    ii, jj, kk = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")
    ii, jj, kk = ii.ravel(), jj.ravel(), kk.ravel()
    tets, tl = [], []
    for perm in itertools.permutations(range(3)):
        corner = [np.zeros_like(ii), np.zeros_like(ii), np.zeros_like(ii)]
        verts = [nid(ii, jj, kk)]
        for axis in perm:
            corner[axis] = corner[axis] + 1
            verts.append(nid(ii + corner[0], jj + corner[1], kk + corner[2]))
        tets.append(np.column_stack(verts))
        tl.append(lab[ii, jj, kk])
    return TetMesh(nodes, np.vstack(tets).astype(np.int64), np.concatenate(tl).astype(np.int64))


def icosphere(radius: float, level: int) -> TetMesh:
    """Sphere (centre at the origin): tetrahedra (centre, triangle) of a subdivided icosahedron."""
    t = (1 + 5**0.5) / 2
    v = [
        (-1, t, 0),
        (1, t, 0),
        (-1, -t, 0),
        (1, -t, 0),
        (0, -1, t),
        (0, 1, t),
        (0, -1, -t),
        (0, 1, -t),
        (t, 0, -1),
        (t, 0, 1),
        (-t, 0, -1),
        (-t, 0, 1),
    ]
    verts = [np.array(p, dtype=np.float64) / np.linalg.norm(p) for p in v]
    faces = [
        (0, 11, 5),
        (0, 5, 1),
        (0, 1, 7),
        (0, 7, 10),
        (0, 10, 11),
        (1, 5, 9),
        (5, 11, 4),
        (11, 10, 2),
        (10, 7, 6),
        (7, 1, 8),
        (3, 9, 4),
        (3, 4, 2),
        (3, 2, 6),
        (3, 6, 8),
        (3, 8, 9),
        (4, 9, 5),
        (2, 4, 11),
        (6, 2, 10),
        (8, 6, 7),
        (9, 8, 1),
    ]
    for _ in range(level):
        cache: dict[tuple[int, int], int] = {}

        def mid(a: int, b: int, cache: dict[tuple[int, int], int] = cache) -> int:
            key = (min(a, b), max(a, b))
            if key not in cache:
                m = verts[a] + verts[b]
                verts.append(m / np.linalg.norm(m))
                cache[key] = len(verts) - 1
            return cache[key]

        new = []
        for a, b, c in faces:
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
            new += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        faces = new
    nodes = np.vstack([np.zeros(3), radius * np.array(verts)])
    tets = np.array([(0, a + 1, b + 1, c + 1) for a, b, c in faces], dtype=np.int64)
    return TetMesh(nodes, tets, np.zeros(len(tets), dtype=np.int64))


def simulate_tets(
    mesh: TetMesh,
    materials: tuple[Layer, ...],
    packets: int,
    seed: int = 1,
    entry_point: tuple[float, float, float] = (0.0, 0.0, 0.0),
    direction: tuple[float, float, float] = (0.0, 0.0, 1.0),
    n_top: float = 1.0,
    n_bottom: float = 1.0,
    n_side: float | None = None,
    n_all: float = 0.0,
    threads: int | None = None,
) -> tuple[Protocol, dict[str, object]]:
    """C++ transport in the mesh; boundary faces get n_top (normal z < −0.5), n_bottom (z > 0.5),
    n_side (else, default: index-matched to the material) or n_all for all."""
    from scattuccino import fast

    if fast._core is None:
        raise RuntimeError("tetrahedral transport needs scattuccino._core")
    core = fast._core
    inp = core.TetInput()
    inp.nodes = np.asarray(mesh.nodes, dtype=np.float64).ravel().tolist()
    inp.tets = np.asarray(mesh.tets, dtype=np.int32).ravel().tolist()
    inp.labels = np.asarray(mesh.labels, dtype=np.int32).tolist()
    mats = []
    for layer in materials:
        m = core.Material()
        m.mua, m.mus, m.n = layer.mua, layer.mus, layer.n
        m.phase = fast._phase(layer.phase_function)
        mats.append(m)
    inp.materials = mats
    inp.n_top, inp.n_bottom = n_top, n_bottom
    inp.n_side = materials[0].n if n_side is None else n_side
    inp.n_all = n_all
    inp.entry_point, inp.direction = list(entry_point), list(direction)
    inp.seed, inp.packets = seed, packets
    inp.threads = threads if threads and threads > 0 else (os.cpu_count() or 1)
    out = core.simulate_tets(inp)
    order = np.argsort(out["id"], kind="stable")
    nm = len(materials)
    proto = Protocol(
        np.asarray(out["side"], dtype=np.int64)[order],
        out["x"][order],
        out["y"][order],
        out["ux"][order],
        out["uy"][order],
        out["uz"][order],
        out["w"][order] / packets,
        np.asarray(out["path"]).reshape(-1, nm)[order],
        np.asarray(out["scatterings"], dtype=np.int64)[order],
        out["z"][order],
        np.array([m.n for m in materials]),
        packets,
        np.array([m.mua for m in materials]),
        np.asarray(out["id"], dtype=np.int64)[order],
    )
    info = {
        "r_specular": float(out["r_specular"]),
        "entry_tet": int(out["entry_tet"]),
        "absorbed": out["absorbed"] / packets,
        "absorbed_tets": out["absorbed_tets"] / packets,
    }
    return proto, info
