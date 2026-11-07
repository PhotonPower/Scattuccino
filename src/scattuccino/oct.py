"""Optical coherence tomography from the photon protocol (stage S9a, intensity model).

Mean (speckle-free) OCT signal as in the Monte-Carlo OCT literature (Yao & Wang 1999): packets
leaving through the top surface into a detector (fibre or objective) at the source position –
radius `r_det` around the entry point, exit direction within the acceptance angle `theta_det`
of −ẑ – contribute at the depth z = OPL/2 (OPL = Σ n_i L_i, the optical path in the medium) with
their weight, convolved with the coherence function of the source (Gaussian, FWHM coherence
length l_c in the same round-trip depth scale).

Single-backscattering class (exactly one scattering event): for a half space with matched index
a packet scattered at depth z into polar angle θ (from −ẑ) leaves with weight μs e^(−μt z)
p(cos Θ = −cos θ) e^(−μt z/cos θ) and z_OPL = z(1 + 1/cos θ)/2 – the exact reference of the tests.

Focusing (Fieldtracing), the wave-optical start in the first mean free paths (BWPM) and the
extended Huygens–Fresnel model follow in S9b/S9c.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np

from scattuccino._types import FloatArray
from scattuccino.transport import Protocol


def detected(p: Protocol, r_det: float, theta_det: float) -> np.ndarray[Any, Any]:
    """Mask of the packets leaving through the top into the detector."""
    return np.asarray((p.side == -1) & (p.rho <= r_det) & (-p.uz >= np.cos(theta_det)))


def depth_opl(p: Protocol) -> FloatArray:
    """z_OPL = OPL/2 in mm."""
    return np.asarray(0.5 * (p.path @ p.n_layer))


def a_scan(
    p: Protocol,
    z: FloatArray,
    r_det: float,
    theta_det: float,
    l_c: float = 0.0,
    mask: np.ndarray[Any, Any] | None = None,
) -> FloatArray:
    """Mean OCT signal per mm of depth on the grid z (uniform). l_c = 0: histogram (no gating);
    else convolution with a Gaussian of FWHM l_c."""
    sel = detected(p, r_det, theta_det)
    if mask is not None:
        sel = sel & mask
    zo = depth_opl(p)[sel]
    w = p.weight[sel]
    dz = float(z[1] - z[0])
    if l_c <= 0:
        edges = np.concatenate([z - dz / 2, [z[-1] + dz / 2]])
        h, _ = np.histogram(zo, bins=edges, weights=w)
        return np.asarray(h / dz)
    sigma = l_c / (2 * np.sqrt(2 * np.log(2)))
    out = np.zeros(z.size)
    for i in range(0, zo.size, 20_000):
        d = (z[None, :] - zo[i : i + 20_000, None]) / sigma
        out += np.sum(w[i : i + 20_000, None] * np.exp(-0.5 * d * d), axis=0)
    return np.asarray(out / (sigma * np.sqrt(2 * np.pi)))


def single_scattering_a_scan(
    z_edges: FloatArray,
    mua: float,
    mus: float,
    pdf_back: Callable[[FloatArray], FloatArray],
    theta_det: float,
    n_z: int = 4000,
    n_theta: int = 400,
) -> FloatArray:
    """Expected single-backscattering signal per depth bin (matched index, half space, detector
    radius unlimited): ∫∫ μs e^(−μt z) p(−cos θ) e^(−μt z/cos θ) 2π sin θ dθ dz binned by
    z_OPL = z(1 + 1/cos θ)/2. pdf_back(μ) is the phase function per sr."""
    mut = mua + mus
    zmax = float(z_edges[-1])
    zz = (np.arange(n_z) + 0.5) * (zmax / n_z)
    tt = (np.arange(n_theta) + 0.5) * (theta_det / n_theta)
    zg, tg = np.meshgrid(zz, tt, indexing="ij")
    c = np.cos(tg)
    w = (
        mus
        * np.exp(-mut * zg)
        * pdf_back(-c)
        * np.exp(-mut * zg / c)
        * 2
        * np.pi
        * np.sin(tg)
        * (zmax / n_z)
        * (theta_det / n_theta)
    )
    zo = zg * (1 + 1 / c) / 2
    h, _ = np.histogram(zo.ravel(), bins=z_edges, weights=w.ravel())
    return np.asarray(h)
