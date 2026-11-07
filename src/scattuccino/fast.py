"""Facade for the C++ core (stage S5): same inputs and results as the Python reference.

`simulate_voxels(medium, packets, …, backend="auto"|"python"|"cpu"|"cuda")` returns a
`scattuccino.voxels.VoxelSimulation`; with the core the packets run in parallel (OpenMP) and the
protocol is sorted by packet number (the Python reference records in exit order; compare by
`protocol.packet`).
"""

from __future__ import annotations

import os
from typing import Any

import numpy as np

from scattuccino.phase import Gegenbauer, HenyeyGreenstein, ModifiedHG, PhaseFunction, Tabulated
from scattuccino.transport import Protocol
from scattuccino.voxels import VoxelMedium, VoxelSimulation
from scattuccino.voxels import simulate_voxels as _python_voxels

try:
    from scattuccino import _core  # type: ignore[attr-defined]

    CORE_AVAILABLE = True
except ImportError:
    _core = None
    CORE_AVAILABLE = False


def core_module() -> Any:
    """The compiled module `scattuccino._core`, or None without the core."""
    return _core


def cuda_available() -> bool:
    """True if the core was built with SCATTUCCINO_ENABLE_CUDA=ON and a GPU is present."""
    return bool(_core is not None and _core.HAS_CUDA and _core.cuda_device_count() > 0)


def _phase(ph: PhaseFunction) -> Any:
    assert _core is not None
    c = _core.Phase()
    if isinstance(ph, HenyeyGreenstein):
        c.kind, c.p1 = _core.PhaseKind.HenyeyGreenstein, ph.g_hg
    elif isinstance(ph, ModifiedHG):
        c.kind, c.p1, c.p2 = _core.PhaseKind.ModifiedHG, ph.g_hg, ph.beta
    elif isinstance(ph, Gegenbauer):
        c.kind, c.p1, c.p2 = _core.PhaseKind.Gegenbauer, ph.g_gk, ph.alpha
    elif isinstance(ph, Tabulated):
        c.kind = _core.PhaseKind.Tabulated
        c.mu, c.p, c.cdf = list(ph.mu), list(ph.p), list(ph.cdf)
    else:
        raise TypeError(f"no C++ form for {type(ph).__name__}")
    return c


def simulate_voxels(
    medium: VoxelMedium,
    packets: int,
    seed: int = 1,
    cos_incident: float = 1.0,
    w_threshold: float = 1e-4,
    roulette_m: float = 10.0,
    backend: str = "auto",
    threads: int | None = None,
    device: int = 0,
) -> VoxelSimulation:
    """backend: "python", "cpu", "cuda" or "auto" (CUDA if available, else CPU, else Python)."""
    if backend not in ("auto", "python", "cpu", "cuda"):
        raise ValueError("backend must be 'auto', 'python', 'cpu' or 'cuda'")
    if backend == "python" or (backend == "auto" and _core is None):
        return _python_voxels(medium, packets, seed, cos_incident, w_threshold, roulette_m)
    if _core is None:
        raise RuntimeError("scattuccino._core is not built")
    lab = np.asarray(medium.labels, dtype=np.int32)
    inp = _core.VoxelInput()
    inp.nx, inp.ny, inp.nz = (int(v) for v in lab.shape)
    inp.dx, inp.dy, inp.dz = medium.voxel
    inp.labels = lab.ravel().tolist()
    mats = []
    for layer in medium.materials:
        m = _core.Material()
        m.mua, m.mus, m.n = layer.mua, layer.mus, layer.n
        m.phase = _phase(layer.phase_function)
        mats.append(m)
    inp.materials = mats
    inp.n_above, inp.n_below = medium.n_above, medium.n_below
    inp.periodic = medium.lateral == "periodic"
    inp.cos_incident, inp.seed, inp.packets = cos_incident, seed, packets
    inp.w_threshold, inp.roulette_m = w_threshold, roulette_m
    inp.threads = threads if threads and threads > 0 else (os.cpu_count() or 1)
    use_cuda = backend == "cuda" or (backend == "auto" and cuda_available())
    if use_cuda:
        if not cuda_available():
            raise RuntimeError("CUDA backend requested but not built or no GPU present")
        out = _core.simulate_voxels_cuda(inp, device)
    else:
        out = _core.simulate_voxels(inp)
    order = np.argsort(out["id"], kind="stable")
    nm = len(medium.materials)
    path = np.asarray(out["path"]).reshape(-1, nm)[order]
    proto = Protocol(
        np.asarray(out["side"], dtype=np.int64)[order],
        out["x"][order],
        out["y"][order],
        out["ux"][order],
        out["uy"][order],
        out["uz"][order],
        out["w"][order] / packets,
        path,
        np.asarray(out["scatterings"], dtype=np.int64)[order],
        out["zmax"][order],
        np.array([m.n for m in medium.materials]),
        packets,
        np.array([m.mua for m in medium.materials]),
        np.asarray(out["id"], dtype=np.int64)[order],
    )
    return VoxelSimulation(
        medium,
        float(out["r_specular"]),
        proto,
        out["absorbed"] / packets,
        np.asarray(out["absorbed_voxels"]).reshape(lab.shape) / packets,
        seed,
    )
