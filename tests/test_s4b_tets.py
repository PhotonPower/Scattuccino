"""S4b: tetrahedral meshes (T6b). Needs the C++ core (skipped otherwise, required in CI).

- the same geometry as voxels split into tetrahedra gives the same result per packet as the voxel
  walk (same random numbers, same Fresnel decisions);
- curved boundary: a glass sphere (faceted, 20 480 tetrahedra), pencil beam at 30° to the local
  normal: specular reflection = Fresnel, the directly transmitted light leaves at the deviation
  2(θᵢ − θₜ) with weight (1 − R₁)(1 − R₂) (geometric optics), energy is conserved.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

from scattuccino import fast
from scattuccino.layers import Layer, fresnel_unpolarised
from scattuccino.phase import ModifiedHG
from scattuccino.voxels import VoxelMedium

if not fast.CORE_AVAILABLE:
    if os.environ.get("SCATTUCCINO_REQUIRE_CORE") == "1":
        raise RuntimeError("SCATTUCCINO_REQUIRE_CORE=1 but scattuccino._core is not built")
    pytest.skip("scattuccino._core not built", allow_module_level=True)

from scattuccino.tets import icosphere, mesh_from_voxels, simulate_tets


def test_tetrahedra_equal_voxels_packet_by_packet() -> None:
    tissue = Layer(0.0, 0.1, 10.0, 0.9, 1.4)
    vessel = Layer(0.0, 2.0, 10.0, 0.0, 1.5, phase=ModifiedHG(0.95, 0.9))
    n = 11
    lab = np.zeros((n, n, n), dtype=np.int64)
    c = np.arange(n) - n // 2
    xx, yy, zz = np.meshgrid(c, c, np.arange(n), indexing="ij")
    lab[(xx**2 + yy**2 + (zz - 4) ** 2) <= 4] = 1
    vox = fast.simulate_voxels(
        VoxelMedium(lab, (0.1, 0.1, 0.1), (tissue, vessel), 1.0, 1.33, lateral="open"),
        20_000,
        seed=3,
        backend="cpu",
    )
    mesh = mesh_from_voxels(lab, (0.1, 0.1, 0.1))
    p, info = simulate_tets(mesh, (tissue, vessel), 20_000, seed=3, n_top=1.0, n_bottom=1.33)
    pv = vox.protocol
    assert pv.packet is not None and p.packet is not None
    o = np.argsort(pv.packet)
    assert np.array_equal(pv.packet[o], p.packet)
    assert np.array_equal(pv.weight[o], p.weight)
    assert np.array_equal(pv.scatterings[o], p.scatterings)
    for a, b in (
        (pv.x[o], p.x),
        (pv.y[o], p.y),
        (pv.ux[o], p.ux),
        (pv.uz[o], p.uz),
        (pv.path[o], p.path),
    ):
        assert np.max(np.abs(a - b)) < 1e-12
    assert info["r_specular"] == vox.r_specular
    assert np.allclose(info["absorbed"], vox.absorbed, rtol=1e-12)  # type: ignore[arg-type]


def test_glass_sphere_against_geometric_optics() -> None:
    mesh = icosphere(1.0, 5)
    tri = mesh.nodes[mesh.tets[:, 1:]]
    centre = tri.mean(axis=1)
    normal = centre / np.linalg.norm(centre, axis=1)[:, None]
    k = int(np.argmax(normal @ np.array([np.sin(np.radians(30)), 0.0, -np.cos(np.radians(30))])))
    a, b, c = tri[k]
    fn = np.cross(b - a, c - a)
    fn = fn / np.linalg.norm(fn) * np.sign(fn @ centre[k])
    theta_i = float(np.arccos(-fn[2]))
    mesh.nodes = mesh.nodes - centre[k]  # the beam enters at the centre of that facet
    p, info = simulate_tets(mesh, (Layer(0.0, 0.0, 0.0, 0.0, 1.5),), 20_000, seed=1, n_all=1.0)
    r1 = fresnel_unpolarised(1.0, 1.5, np.array([np.cos(theta_i)]))[0][0]
    assert abs(info["r_specular"] - r1) < 1e-15
    theta_t = float(np.arcsin(np.sin(theta_i) / 1.5))
    r2 = fresnel_unpolarised(1.5, 1.0, np.array([np.cos(theta_t)]))[0][0]
    angle = np.arccos(np.clip(p.uz, -1, 1))
    direct = np.abs(angle - 2 * (theta_i - theta_t)) < 0.05
    expected = (1 - r1) * (1 - r2)
    assert abs(p.weight[direct].sum() - expected) < 4 * np.sqrt(expected * (1 - expected) / 20_000)
    assert abs(p.weight.sum() + info["r_specular"] - 1) < 1e-12
