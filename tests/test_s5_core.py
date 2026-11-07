"""S5a: C++ core against the Python reference, packet by packet (T13).

Skipped without the core unless SCATTUCCINO_REQUIRE_CORE=1 (CI).
"""

from __future__ import annotations

import os

import numpy as np
import pytest

from scattuccino import fast
from scattuccino.layers import Layer
from scattuccino.phase import Gegenbauer, HenyeyGreenstein, ModifiedHG, PhaseFunction, mie_phase
from scattuccino.random import CounterRandom
from scattuccino.voxels import VoxelMedium

if not fast.CORE_AVAILABLE:
    if os.environ.get("SCATTUCCINO_REQUIRE_CORE") == "1":
        raise RuntimeError("SCATTUCCINO_REQUIRE_CORE=1 but scattuccino._core is not built")
    pytest.skip("scattuccino._core not built", allow_module_level=True)

from scattuccino import _core  # type: ignore[attr-defined]

MIE = mie_phase(1.59 / 1.33 + 0j, 2 * np.pi * 0.5 * 1.33 / 0.633)[0]


def _inclusion() -> VoxelMedium:
    n = 21
    lab = np.zeros((n, n, n), dtype=np.int64)
    c = np.arange(n) - n // 2
    xx, yy, zz = np.meshgrid(c, c, np.arange(n), indexing="ij")
    lab[(xx**2 + yy**2 + (zz - 6) ** 2) <= 9] = 1
    return VoxelMedium(
        lab,
        (0.1, 0.1, 0.1),
        (
            Layer(0.0, 0.1, 10.0, 0.9, 1.4),
            Layer(0.0, 2.0, 10.0, 0.0, 1.4, phase=ModifiedHG(0.95, 0.9)),
        ),
        lateral="open",
    )


def _layers() -> VoxelMedium:
    lab = np.zeros((4, 4, 15), dtype=np.int64)
    lab[:, :, 5:] = 1
    return VoxelMedium(
        lab,
        (0.5, 0.5, 0.1),
        (
            Layer(0.0, 0.2, 8.0, 0.0, 1.37, phase=MIE),
            Layer(0.0, 0.05, 5.0, 0.0, 1.52, phase=Gegenbauer(0.8, 1.0)),
        ),
        1.0,
        1.33,
    )


def test_random_numbers_are_identical() -> None:
    r = CounterRandom(3, 100)
    for pid in (0, 5, 99):
        py = [r.uniform(np.array([pid]))[0] for _ in range(4)]
        assert np.array_equal(_core.uniforms(3, pid, 4), py)


@pytest.mark.parametrize(
    "phase", [HenyeyGreenstein(0.9), ModifiedHG(0.9, 0.8), Gegenbauer(0.8, 1.0), MIE]
)
def test_phase_sampling_is_identical(phase: PhaseFunction) -> None:
    xi = np.random.default_rng(1).random(10_000)
    c = fast._phase(phase)
    # identical up to the last bit of cbrt/pow in libm vs. NumPy
    got = np.array([c.sample(float(v)) for v in xi])
    assert np.max(np.abs(got - phase.sample(xi))) < 1e-15


@pytest.mark.parametrize(("medium", "cos_i"), [(_inclusion, 1.0), (_layers, float(np.cos(0.5)))])
def test_core_equals_python_packet_by_packet(medium: object, cos_i: float) -> None:
    med = medium()  # type: ignore[operator]
    a = fast.simulate_voxels(med, 10_000, seed=4, cos_incident=cos_i, backend="python")
    b = fast.simulate_voxels(med, 10_000, seed=4, cos_incident=cos_i, backend="cpu", threads=3)
    pa, pb = a.protocol, b.protocol
    assert pa.packet is not None and pb.packet is not None
    o = np.argsort(pa.packet)
    assert np.array_equal(pa.packet[o], pb.packet)
    # same decisions for every packet: sides, weights and numbers of events bit for bit;
    # positions, directions and paths up to the accumulated last bits of the elementary functions
    # (measured ≤ 7e-13 after hundreds of events)
    for name in ("side", "weight", "scatterings"):
        assert np.array_equal(getattr(pa, name)[o], getattr(pb, name)), name
    for name in ("x", "y", "ux", "uy", "uz", "path", "max_depth"):
        assert np.allclose(getattr(pa, name)[o], getattr(pb, name), rtol=1e-10, atol=1e-10), name
    assert np.allclose(a.absorbed, b.absorbed, rtol=1e-12)
    assert np.allclose(a.absorbed_voxels, b.absorbed_voxels, rtol=1e-12, atol=1e-18)


def test_threads_do_not_change_the_result() -> None:
    med = _inclusion()
    one = fast.simulate_voxels(med, 5_000, seed=6, backend="cpu", threads=1)
    many = fast.simulate_voxels(med, 5_000, seed=6, backend="cpu", threads=4)
    assert np.array_equal(one.protocol.weight, many.protocol.weight)
