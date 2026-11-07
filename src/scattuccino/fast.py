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
    record: bool = True,
    score_voxels: bool = True,
    precision: str = "double",
) -> VoxelSimulation:
    """backend: "python", "cpu", "cuda" or "auto" (CUDA if available, else CPU, else Python).

    precision: "double" (default, reference) or "single" (float32 for direction, step lengths,
    phase function, Fresnel and position; weights, sums and random numbers stay double). The
    Python reference is double only.
    """
    if precision not in ("double", "single"):
        raise ValueError("precision must be 'double' or 'single'")
    if backend not in ("auto", "python", "cpu", "cuda"):
        raise ValueError("backend must be 'auto', 'python', 'cpu' or 'cuda'")
    if precision == "single" and (backend == "python" or (backend == "auto" and _core is None)):
        raise ValueError("precision='single' needs the C++ core (backend 'cpu' or 'cuda')")
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
    inp.record, inp.score_voxels = record, score_voxels
    inp.precision = 0 if precision == "double" else 1
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
    n = float(packets)
    totals = {
        "r_diffuse": out["sum_r"] / n,
        "t_diffuse": out["sum_t"] / n,
        "side": out["sum_side"] / n,
        "n_exit": float(out["n_exit"]),
        "n_roulette": float(out["n_roulette"]),
        "r_err": float(np.sqrt(max(out["sum_r2"] / n - (out["sum_r"] / n) ** 2, 0.0) / n)),
        "t_err": float(np.sqrt(max(out["sum_t2"] / n - (out["sum_t"] / n) ** 2, 0.0) / n)),
    }
    return VoxelSimulation(
        medium,
        float(out["r_specular"]),
        proto,
        out["absorbed"] / packets,
        np.asarray(out["absorbed_voxels"]).reshape(lab.shape) / packets,
        seed,
        totals,
    )


def simulate_polarised_slab(
    d: float,
    mua: float,
    mus: float,
    table: Any,
    stokes_in: Any,
    packets: int,
    seed: int = 1,
    n_slab: float = 1.0,
    n_above: float = 1.0,
    n_below: float = 1.0,
    cos_incident: float = 1.0,
    w_threshold: float = 1e-4,
    roulette_m: float = 10.0,
    threads: int | None = None,
) -> Any:
    """Polarised slab on the C++ core (S6c), same inputs and result as
    `scattuccino.polarisation.simulate_polarised` and the same packets (counter-based random
    numbers): a voxel grid of one voxel, laterally periodic and wide; d = np.inf → 10⁹ mm."""
    from scattuccino.polarisation import PolarisedResult

    if _core is None:
        raise RuntimeError("scattuccino._core is not built")
    depth = d if np.isfinite(d) else 1e9
    inp = _core.VoxelInput()
    inp.nx = inp.ny = inp.nz = 1
    inp.dx = inp.dy = 1e9
    inp.dz = depth
    inp.labels = [0]
    mat = _core.Material()
    mat.mua, mat.mus, mat.n = mua, mus, n_slab
    mat.phase = _phase(table.phase)
    inp.materials = [mat]
    inp.n_above, inp.n_below, inp.periodic = n_above, n_below, True
    inp.cos_incident, inp.seed, inp.packets = cos_incident, seed, packets
    inp.w_threshold, inp.roulette_m = w_threshold, roulette_m
    inp.threads = threads if threads and threads > 0 else (os.cpu_count() or 1)
    mt = _core.MuellerTable()
    mt.mu, mt.m12, mt.m33, mt.m34 = (
        list(np.asarray(a)) for a in (table.mu, table.m12, table.m33, table.m34)
    )
    mt.pure = bool(table.pure)
    pin = _core.PolInput()
    pin.voxels = inp
    pin.mueller = [mt]
    pin.stokes_in = [float(x) for x in stokes_in]
    out = _core.simulate_voxels_polarised(pin)
    order = np.argsort(out["id"], kind="stable")
    res = PolarisedResult(
        np.asarray(out["side"], dtype=np.int64)[order],
        np.asarray(out["u"]).reshape(-1, 3)[order],
        np.asarray(out["v"]).reshape(-1, 3)[order],
        np.asarray(out["stokes"]).reshape(-1, 4)[order],
        np.asarray(out["w"])[order] / packets,
        np.asarray(out["scatterings"], dtype=np.int64)[order],
        packets,
    )
    res.specular = (float(out["specular_r"]), np.asarray(out["specular_stokes"]))
    return res
