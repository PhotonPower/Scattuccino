"""Detectors: evaluations of the photon protocol (stage S2). Nothing is recomputed.

All results are per incident power 1, with standard errors from the per-packet contributions
(weights are stored divided by the number of packets, so sums are means).
"""

from __future__ import annotations

import numpy as np

from scattuccino._types import FloatArray
from scattuccino.transport import Protocol


def _select(p: Protocol, side: int, rho_lo: float, rho_hi: float) -> np.ndarray:  # type: ignore[type-arg]
    rho = p.rho
    return np.asarray((p.side == side) & (rho >= rho_lo) & (rho < rho_hi))


def reflectance_vs_rho(
    p: Protocol, edges: FloatArray, side: int = -1, weight: FloatArray | None = None
) -> tuple[FloatArray, FloatArray]:
    """R(ρ) in 1/mm² on annuli [edges_k, edges_k+1) and its standard error."""
    w = p.weight if weight is None else weight
    sel = p.side == side
    rho = p.rho[sel]
    ws = w[sel]
    area = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
    s1, _ = np.histogram(rho, bins=edges, weights=ws)
    s2, _ = np.histogram(rho, bins=edges, weights=(ws * p.packets) ** 2)
    mean = s1
    err = np.sqrt(np.maximum(s2 / p.packets - mean**2, 0.0) / p.packets)
    return np.asarray(mean / area), np.asarray(err / area)


def tpsf(
    p: Protocol, rho_lo: float, rho_hi: float, t_edges: FloatArray, side: int = -1
) -> FloatArray:
    """Time-resolved reflectance (per ps, per mm²) in the annulus, histogram over t_edges (ps)."""
    sel = _select(p, side, rho_lo, rho_hi)
    h, _ = np.histogram(p.time_ps[sel], bins=t_edges, weights=p.weight[sel])
    area = np.pi * (rho_hi**2 - rho_lo**2)
    return np.asarray(h / np.diff(t_edges) / area)


def frequency_domain(
    p: Protocol, rho_lo: float, rho_hi: float, f_ghz: float, side: int = -1
) -> complex:
    """Σ w e^(−iωt) per mm² in the annulus (modulation frequency f in GHz, t in ps)."""
    sel = _select(p, side, rho_lo, rho_hi)
    omega = 2 * np.pi * f_ghz * 1e-3  # rad/ps
    area = np.pi * (rho_hi**2 - rho_lo**2)
    return complex(np.sum(p.weight[sel] * np.exp(-1j * omega * p.time_ps[sel])) / area)


def mean_partial_paths(p: Protocol, rho_lo: float, rho_hi: float, side: int = -1) -> FloatArray:
    """⟨L_i⟩ per layer of the detected light (weight-averaged), mm. ⟨L⟩/ρ is the DPF."""
    sel = _select(p, side, rho_lo, rho_hi)
    w = p.weight[sel]
    return np.asarray((w[:, None] * p.path[sel]).sum(axis=0) / w.sum())
