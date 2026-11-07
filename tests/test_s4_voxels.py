"""S4a: voxel geometry (T6): invariance, layered media in voxels, open boundaries, inclusions."""

from __future__ import annotations

import numpy as np

from scattuccino.adding_doubling import slab
from scattuccino.layers import Layer, Stack
from scattuccino.random import CounterRandom
from scattuccino.transport import simulate
from scattuccino.voxels import VoxelMedium, simulate_voxels

TISSUE = Layer(0.0, 0.1, 10.0, 0.9, 1.4)


def _rt(sim: object) -> tuple[float, float]:
    p = sim.protocol  # type: ignore[attr-defined]
    return float(p.weight[p.side == -1].sum()), float(p.weight[p.side == 1].sum())


def test_counter_random_numbers() -> None:
    r = CounterRandom(7, 200_000)
    u = r.uniform(np.arange(200_000))
    v = r.uniform(np.arange(200_000))
    assert abs(u.mean() - 0.5) < 5 * np.sqrt(1 / 12 / u.size)
    assert abs(np.var(u) * 12 - 1) < 0.01
    assert abs(np.corrcoef(u, v)[0, 1]) < 5 / np.sqrt(u.size)
    assert 0 < u.min() and u.max() < 1
    # a packet's numbers do not depend on the other packets
    r2 = CounterRandom(7, 200_000)
    assert r2.uniform(np.array([12345]))[0] == u[12345]


def test_refining_voxels_and_batches_changes_nothing() -> None:
    coarse = VoxelMedium(np.zeros((1, 1, 1), dtype=np.int64), (1.0, 1.0, 1.0), (TISSUE,))
    fine = VoxelMedium(np.zeros((10, 10, 10), dtype=np.int64), (0.1, 0.1, 0.1), (TISSUE,))
    a = simulate_voxels(coarse, 20_000, seed=3)
    b = simulate_voxels(fine, 20_000, seed=3, batch=7_000)
    ra, ta = _rt(a)
    rb, tb = _rt(b)
    assert abs(ra - rb) < 1e-12 and abs(ta - tb) < 1e-12


def test_single_layer_against_adding_doubling() -> None:
    med = VoxelMedium(np.zeros((4, 4, 10), dtype=np.int64), (0.25, 0.25, 0.1), (TISSUE,))
    sim = simulate_voxels(med, 100_000, seed=5)
    ad = slab(1.0, 0.1, 10.0, 0.9, 1.4)
    p = sim.protocol
    for side, key in ((-1, "r_diffuse"), (1, "t_diffuse")):
        w = np.where(p.side == side, p.weight * p.packets, 0.0)
        mean = w.sum() / p.packets
        err = np.sqrt((np.sum(w**2) / p.packets - mean**2) / p.packets)
        assert abs(mean - ad[key]) < 4 * err


def test_two_layers_with_index_step_against_layered_transport() -> None:
    top = Layer(0.0, 0.2, 8.0, 0.8, 1.37)
    bottom = Layer(0.0, 0.05, 5.0, 0.9, 1.52)
    lab = np.zeros((2, 2, 15), dtype=np.int64)
    lab[:, :, 5:] = 1  # 0.5 mm over 1.0 mm
    vox = simulate_voxels(
        VoxelMedium(lab, (0.5, 0.5, 0.1), (top, bottom), 1.0, 1.33),
        60_000,
        seed=9,
        cos_incident=np.cos(0.5),
    )
    lay = simulate(
        Stack((Layer(0.5, 0.2, 8.0, 0.8, 1.37), Layer(1.0, 0.05, 5.0, 0.9, 1.52)), 1.0, 1.33),
        60_000,
        seed=10,
        cos_incident=np.cos(0.5),
    )
    assert abs(vox.r_specular - lay.r_specular) < 1e-15
    for got, ref in zip(_rt(vox), _rt(lay), strict=True):
        assert abs(got - ref) < 0.006  # two independent runs, 4σ of the difference ≈ 0.006
    assert np.allclose(vox.absorbed.sum() + sum(_rt(vox)) + vox.r_specular, 1.0, atol=1e-4)


def test_open_boundaries_and_absorbing_inclusion() -> None:
    n = 21
    lab = np.zeros((n, n, n), dtype=np.int64)
    c = np.arange(n) - n // 2
    xx, yy, zz = np.meshgrid(c, c, np.arange(n), indexing="ij")
    inclusion = (xx**2 + yy**2 + (zz - 6) ** 2) <= 9  # sphere r = 0.3 mm at 0.65 mm depth
    lab[inclusion] = 1
    blood = Layer(0.0, 2.0, 10.0, 0.98, 1.4)
    plain = simulate_voxels(
        VoxelMedium(np.zeros_like(lab), (0.1, 0.1, 0.1), (TISSUE, blood), lateral="open"),
        30_000,
        seed=2,
    )
    vessel = simulate_voxels(
        VoxelMedium(lab, (0.1, 0.1, 0.1), (TISSUE, blood), lateral="open"), 30_000, seed=2
    )
    for sim in (plain, vessel):
        total = sim.r_specular + sim.protocol.weight.sum() + sim.absorbed.sum()
        assert abs(total - 1) < 1e-4  # side losses are part of the protocol (side = 0)
        assert sim.protocol.weight[sim.protocol.side == 0].sum() > 0
    assert abs(vessel.absorbed_voxels[inclusion].sum() - vessel.absorbed[1]) < 1e-12
    assert _rt(vessel)[0] < _rt(plain)[0] - 0.005  # the absorber darkens the reflection
