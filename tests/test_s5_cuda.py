"""S5b: CUDA kernel against the CPU core, packet by packet (T14).

Runs only where the core was built with CUDA and a GPU is present; with
SCATTUCCINO_REQUIRE_CUDA=1 a missing CUDA backend fails (large machine).
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pytest

from scattuccino import fast

if not fast.cuda_available():
    if os.environ.get("SCATTUCCINO_REQUIRE_CUDA") == "1":
        raise RuntimeError("SCATTUCCINO_REQUIRE_CUDA=1 but the CUDA backend is not available")
    pytest.skip("CUDA backend not built or no GPU", allow_module_level=True)

sys.path.insert(0, str(Path(__file__).parent))
from test_s5_core import _inclusion, _layers


@pytest.mark.parametrize(("medium", "cos_i"), [(_inclusion, 1.0), (_layers, float(np.cos(0.5)))])
def test_cuda_equals_cpu_packet_by_packet(medium: object, cos_i: float) -> None:
    med = medium()  # type: ignore[operator]
    a = fast.simulate_voxels(med, 50_000, seed=4, cos_incident=cos_i, backend="cpu")
    b = fast.simulate_voxels(med, 50_000, seed=4, cos_incident=cos_i, backend="cuda")
    pa, pb = a.protocol, b.protocol
    assert pa.packet is not None and pb.packet is not None
    # the GPU math library rounds log/cos/sin/pow/cbrt differently in the last bits; a decision
    # may flip in rare packets – count them instead of requiring all
    common, ia, ib = np.intersect1d(pa.packet, pb.packet, return_indices=True)
    same = (
        (pa.side[ia] == pb.side[ib])
        & (pa.scatterings[ia] == pb.scatterings[ib])
        & np.isclose(pa.weight[ia], pb.weight[ib], rtol=1e-9, atol=0)
    )
    assert same.mean() > 0.999
    assert abs(pa.weight.sum() - pb.weight.sum()) < 1e-3 * pa.weight.sum()
    assert np.allclose(a.absorbed, b.absorbed, rtol=1e-3)


def test_cuda_large_run_statistics() -> None:
    med = _inclusion()
    cpu = fast.simulate_voxels(med, 200_000, seed=8, backend="cpu")
    gpu = fast.simulate_voxels(med, 200_000, seed=8, backend="cuda")
    r_cpu = cpu.protocol.weight[cpu.protocol.side == -1].sum()
    r_gpu = gpu.protocol.weight[gpu.protocol.side == -1].sum()
    assert abs(r_cpu - r_gpu) < 1e-3 * r_cpu


def test_cuda_without_protocol_and_voxel_scoring() -> None:
    # S5c: sums without the protocol and without voxel scoring equal the full run
    med = _inclusion()
    full = fast.simulate_voxels(med, 50_000, seed=3, backend="cuda")
    lean = fast.simulate_voxels(
        med, 50_000, seed=3, backend="cuda", record=False, score_voxels=False
    )
    assert full.totals is not None and lean.totals is not None
    for key in ("r_diffuse", "t_diffuse", "side"):
        assert abs(full.totals[key] - lean.totals[key]) < 1e-12
    assert np.allclose(full.absorbed, lean.absorbed, rtol=1e-12)
    assert lean.protocol.weight.size == 0
