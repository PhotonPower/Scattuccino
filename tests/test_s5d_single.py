"""S5d: single-precision option of the voxel transport (T15).

CPU in single is the reference for the GPU in single; the acceptance of single against double is
statistical (decision D7). Skipped without the core unless SCATTUCCINO_REQUIRE_CORE=1; the GPU test
needs CUDA (SCATTUCCINO_REQUIRE_CUDA=1 makes a missing backend fail).
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pytest

from scattuccino import fast
from scattuccino.layers import Layer
from scattuccino.voxels import VoxelMedium

if not fast.CORE_AVAILABLE:
    if os.environ.get("SCATTUCCINO_REQUIRE_CORE") == "1":
        raise RuntimeError("SCATTUCCINO_REQUIRE_CORE=1 but scattuccino._core is not built")
    pytest.skip("scattuccino._core not built", allow_module_level=True)

sys.path.insert(0, str(Path(__file__).parent))
from test_s5_core import _inclusion, _layers


def _weak() -> VoxelMedium:
    """Weakly absorbing slab: μa = 0.01, μs = 10 /mm, g = 0.9, n = 1.4, 5 mm, periodic."""
    return VoxelMedium(
        np.zeros((1, 1, 50), dtype=np.int64),
        (100.0, 100.0, 0.1),
        (Layer(0.0, 0.01, 10.0, 0.9, 1.4),),
    )


CASES = [(_inclusion, 1.0), (_layers, float(np.cos(0.5))), (_weak, 1.0)]
RING_EDGES = np.linspace(0.0, 4.0, 9)


def test_precision_argument() -> None:
    med = _inclusion()
    with pytest.raises(ValueError):
        fast.simulate_voxels(med, 10, backend="python", precision="single")
    with pytest.raises(ValueError):
        fast.simulate_voxels(med, 10, backend="cpu", precision="half")


@pytest.mark.parametrize(("medium", "cos_i"), CASES[:2])
def test_double_equals_python_reference(medium: object, cos_i: float) -> None:
    med = medium()  # type: ignore[operator]
    a = fast.simulate_voxels(med, 10_000, seed=4, cos_incident=cos_i, backend="python")
    b = fast.simulate_voxels(
        med, 10_000, seed=4, cos_incident=cos_i, backend="cpu", precision="double"
    )
    pa, pb = a.protocol, b.protocol
    assert pa.packet is not None and pb.packet is not None
    o = np.argsort(pa.packet)
    assert np.array_equal(pa.packet[o], pb.packet)
    for name in ("side", "weight", "scatterings"):
        assert np.array_equal(getattr(pa, name)[o], getattr(pb, name)), name


def _batch_stats(
    med: VoxelMedium, cos_i: float, precision: str, seed0: int, nb: int = 20, n: int = 10_000
) -> tuple[np.ndarray, np.ndarray]:  # type: ignore[type-arg]
    """Mean and standard error over nb independent batches of n packets (2·10⁵ in all).

    Quantities: R_d, T_d, lateral loss, absorption per material, R(ρ) on 8 rings 0–4 mm.
    """
    rows = []
    for b in range(nb):
        s = fast.simulate_voxels(
            med, n, seed=seed0 + b, cos_incident=cos_i, backend="cpu", precision=precision
        )
        assert s.totals is not None
        p = s.protocol
        top = p.side == -1
        ring, _ = np.histogram(p.rho[top], bins=RING_EDGES, weights=p.weight[top])
        rows.append(
            np.concatenate(
                [[s.totals["r_diffuse"], s.totals["t_diffuse"], s.totals["side"]], s.absorbed, ring]
            )
        )
    v = np.array(rows)
    return v.mean(axis=0), v.std(axis=0, ddof=1) / np.sqrt(nb)


@pytest.mark.parametrize(("medium", "cos_i"), CASES)
def test_single_agrees_with_double_statistically(medium: object, cos_i: float) -> None:
    med = medium()  # type: ignore[operator]
    m_d, e_d = _batch_stats(med, cos_i, "double", 1000)
    m_s, e_s = _batch_stats(med, cos_i, "single", 2000)
    sigma = np.hypot(e_d, e_s)
    dev = np.abs(m_s - m_d)
    # 4σ of the combined standard error; quantities that are zero in both (rings outside the
    # lateral extent) agree trivially
    worst = np.max(dev / np.where(sigma > 0, sigma, np.inf))
    print(f"max deviation {worst:.2f} σ")
    assert np.all(dev <= 4 * sigma + 1e-12), (worst, dev, sigma)


@pytest.mark.parametrize(("medium", "cos_i"), CASES)
def test_single_conserves_energy_and_loses_no_packet(medium: object, cos_i: float) -> None:
    med = medium()  # type: ignore[operator]
    n = 200_000
    s = fast.simulate_voxels(med, n, seed=7, cos_incident=cos_i, backend="cpu", precision="single")
    t = s.totals
    assert t is not None
    total = s.r_specular + t["r_diffuse"] + t["t_diffuse"] + t["side"] + s.absorbed.sum()
    assert abs(total - 1.0) < 1e-4, total
    # every packet either left or was killed by the roulette; the protocol holds the exits
    assert t["n_exit"] + t["n_roulette"] == n
    assert s.protocol.weight.size == t["n_exit"]


def test_gpu_single_equals_cpu_single_packet_by_packet() -> None:
    if not fast.cuda_available():
        if os.environ.get("SCATTUCCINO_REQUIRE_CUDA") == "1":
            raise RuntimeError("SCATTUCCINO_REQUIRE_CUDA=1 but the CUDA backend is not available")
        pytest.skip("CUDA backend not built or no GPU")
    for medium, cos_i in CASES:
        med = medium()
        kw = {"seed": 4, "cos_incident": cos_i, "precision": "single"}
        n = 50_000 if medium is not _weak else 20_000
        a = fast.simulate_voxels(med, n, backend="cpu", **kw)  # type: ignore[arg-type]
        b = fast.simulate_voxels(med, n, backend="cuda", **kw)  # type: ignore[arg-type]
        pa, pb = a.protocol, b.protocol
        assert pa.packet is not None and pb.packet is not None
        _, ia, ib = np.intersect1d(pa.packet, pb.packet, return_indices=True)
        same = (pa.side[ia] == pb.side[ib]) & (pa.scatterings[ia] == pb.scatterings[ib])
        frac = same.mean()
        sums = abs(pa.weight.sum() - pb.weight.sum()) / pa.weight.sum()
        print(f"{medium.__name__}: same packets {frac:.5f}, sum deviation {sums:.2e}")
        assert frac >= 0.99
        assert sums < 1e-3
        assert np.allclose(a.absorbed, b.absorbed, rtol=1e-3)
        assert a.totals is not None and b.totals is not None
        assert a.totals["n_exit"] + a.totals["n_roulette"] == n
        assert b.totals["n_exit"] + b.totals["n_roulette"] == n
