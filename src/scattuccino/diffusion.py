"""Diffusion-theory references for stage S2 (approximations; valid for μs′ ≫ μa and far from the
source). Extrapolated-boundary dipole (Farrell, Patterson, Wilson 1992; Patterson, Chance, Wilson
1989; boundary coefficient after Groenhuis/Egan–Hilgeman as used by Farrell 1992).

Units: mm, 1/mm, ps.
"""

from __future__ import annotations

import numpy as np

from scattuccino._types import FloatArray
from scattuccino.transport import C_MM_PER_PS


def boundary_a(n_rel: float) -> float:
    """A = (1 + r_d)/(1 − r_d), r_d = −1.440/n² + 0.710/n + 0.668 + 0.0636 n, n = n_tissue/n_out."""
    if n_rel == 1.0:
        return 1.0
    r_d = -1.440 / n_rel**2 + 0.710 / n_rel + 0.668 + 0.0636 * n_rel
    return (1 + r_d) / (1 - r_d)


def _params(mua: float, musp: float, n_rel: float) -> tuple[float, float, float, float]:
    mutp = mua + musp
    d = 1.0 / (3.0 * mutp)
    mueff = float(np.sqrt(mua / d))
    z0 = 1.0 / mutp
    zb = 2.0 * boundary_a(n_rel) * d
    return d, mueff, z0, zb


def reflectance_rho(rho: FloatArray, mua: float, musp: float, n_rel: float = 1.0) -> FloatArray:
    """Farrell's R(ρ) (1/mm²): flux through the surface of the extrapolated-boundary dipole."""
    _, mueff, z0, zb = _params(mua, musp, n_rel)
    r1 = np.sqrt(z0**2 + rho**2)
    r2 = np.sqrt((z0 + 2 * zb) ** 2 + rho**2)
    return np.asarray(
        (
            z0 * (mueff + 1 / r1) * np.exp(-mueff * r1) / r1**2
            + (z0 + 2 * zb) * (mueff + 1 / r2) * np.exp(-mueff * r2) / r2**2
        )
        / (4 * np.pi)
    )


def fluence_surface(rho: FloatArray, mua: float, musp: float, n_rel: float = 1.0) -> FloatArray:
    """Φ(ρ, z = 0) (1/mm²) of the extrapolated-boundary dipole."""
    d, mueff, z0, zb = _params(mua, musp, n_rel)
    r1 = np.sqrt(z0**2 + rho**2)
    r2 = np.sqrt((z0 + 2 * zb) ** 2 + rho**2)
    return np.asarray((np.exp(-mueff * r1) / r1 - np.exp(-mueff * r2) / r2) / (4 * np.pi * d))


def fluence(
    rho: FloatArray, z: FloatArray, mua: float, musp: float, n_rel: float = 1.0
) -> FloatArray:
    """Φ(ρ, z) (1/mm²) inside the medium, extrapolated-boundary dipole (source at z0, image at
    −z0 − 2zb)."""
    d, mueff, z0, zb = _params(mua, musp, n_rel)
    r1 = np.sqrt((z - z0) ** 2 + rho**2)
    r2 = np.sqrt((z + z0 + 2 * zb) ** 2 + rho**2)
    return np.asarray((np.exp(-mueff * r1) / r1 - np.exp(-mueff * r2) / r2) / (4 * np.pi * d))


def fluence_laterally_integrated(
    z: FloatArray, mua: float, musp: float, n_rel: float = 1.0
) -> FloatArray:
    """∫ Φ(ρ, z) 2πρ dρ (1/mm) of the extrapolated-boundary dipole (plane-source form)."""
    d, mueff, z0, zb = _params(mua, musp, n_rel)
    return np.asarray(
        (np.exp(-mueff * np.abs(z - z0)) - np.exp(-mueff * (z + z0 + 2 * zb))) / (2 * d * mueff)
    )


def kp_coefficients(n_rel: float) -> tuple[float, float]:
    """C1, C2 of R = C1 Φ + C2 j (Kienle & Patterson 1997): radiance in P1 form leaving through
    the Fresnel boundary, C1 = ½∫(1 − R_F)μ dμ, C2 = (3/2)∫(1 − R_F)μ² dμ (from the tissue)."""
    from scattuccino.layers import fresnel_unpolarised

    x, w = np.polynomial.legendre.leggauss(400)
    if n_rel > 1:
        mu_c = np.sqrt(1 - 1 / n_rel**2)
        mu = mu_c + (1 - mu_c) * 0.5 * (x + 1)
        wq = (1 - mu_c) * 0.5 * w
    else:
        mu, wq = 0.5 * (x + 1), 0.5 * w
    rf, _ = fresnel_unpolarised(n_rel, 1.0, mu)
    return float(0.5 * np.sum(wq * (1 - rf) * mu)), float(1.5 * np.sum(wq * (1 - rf) * mu**2))


def reflectance_rho_kp(rho: FloatArray, mua: float, musp: float, n_rel: float = 1.0) -> FloatArray:
    """Kienle–Patterson R(ρ) = C1 Φ(ρ, 0) + C2 j(ρ) (1/mm²)."""
    c1, c2 = kp_coefficients(n_rel)
    return np.asarray(
        c1 * fluence_surface(rho, mua, musp, n_rel) + c2 * reflectance_rho(rho, mua, musp, n_rel)
    )


def tpsf(
    rho: float, t_ps: FloatArray, mua: float, musp: float, n: float, n_rel: float | None = None
) -> FloatArray:
    """Patterson's time-resolved reflectance R(ρ, t) (1/(mm² ps)), extrapolated boundary."""
    n_rel = n if n_rel is None else n_rel
    d, _, z0, zb = _params(mua, musp, n_rel)
    c = C_MM_PER_PS / n
    t = np.asarray(t_ps, dtype=np.float64)
    pre = (4 * np.pi * d * c) ** -1.5 * t**-2.5 * np.exp(-mua * c * t) / 2
    zp = z0 + 2 * zb
    return np.asarray(
        pre
        * (
            z0 * np.exp(-(rho**2 + z0**2) / (4 * d * c * t))
            + zp * np.exp(-(rho**2 + zp**2) / (4 * d * c * t))
        )
    )
