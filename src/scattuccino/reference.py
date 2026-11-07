"""Independent references for stage S1 (docs/spezifikation.md §4), not using the Monte Carlo.

- `slab_without_scattering`: Beer-Lambert slab with incoherent multiple reflections.
- `chandrasekhar_h`, `semi_infinite_isotropic_reflectance`: Chandrasekhar's H function and the
  total reflectance R = 1 − H(μ₀)√(1 − a) of an isotropically scattering half space with matched
  index (Chandrasekhar, Radiative Transfer, 1960).
"""

from __future__ import annotations

import numpy as np

from scattuccino._types import FloatArray


def slab_without_scattering(
    mua: float, d: float, n_slab: float, n_out: float = 1.0
) -> tuple[float, float]:
    """(R, T) at normal incidence, same medium on both sides, incoherent reflections."""
    r = ((n_out - n_slab) / (n_out + n_slab)) ** 2
    att = np.exp(-mua * d)
    den = 1 - r * r * att * att
    return float(r + (1 - r) ** 2 * r * att * att / den), float((1 - r) ** 2 * att / den)


def chandrasekhar_h(
    albedo: float, mu: FloatArray, n_quad: int = 200, tol: float = 1e-14
) -> FloatArray:
    """H(μ) for isotropic scattering, characteristic function ψ = a/2, by iteration of
    1/H(μ) = 1 − (a/2) μ ∫₀¹ H(μ′)/(μ + μ′) dμ′ (Gauss-Legendre in μ′; converges for a < 1)."""
    x, w = np.polynomial.legendre.leggauss(n_quad)
    mq, wq = 0.5 * (x + 1), 0.5 * w
    h = np.ones_like(mq)
    for _ in range(10_000):
        h_new = 1.0 / (
            1.0 - 0.5 * albedo * mq * np.sum(wq * h / (mq[:, None] + mq[None, :]), axis=1)
        )
        if np.max(np.abs(h_new - h)) < tol:
            h = h_new
            break
        h = h_new
    mu = np.atleast_1d(np.asarray(mu, dtype=np.float64))
    return np.asarray(
        1.0 / (1.0 - 0.5 * albedo * mu * np.sum(wq * h / (mu[:, None] + mq[None, :]), axis=1))
    )


def semi_infinite_isotropic_reflectance(albedo: float, mu0: float = 1.0) -> float:
    """Total (diffuse) reflectance of an isotropic half space, matched index: 1 − H(μ₀)√(1 − a)."""
    return float(1.0 - chandrasekhar_h(albedo, np.array([mu0]))[0] * np.sqrt(1.0 - albedo))
