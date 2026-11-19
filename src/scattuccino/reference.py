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


# S8b: laser Doppler model of Bonner and Nossal (bonner1981). Moving cells of radius a (mm) with
# the Rayleigh-Gans structure factor S(Qa) ≈ exp(−2ξ(Qa)²) (Eq. 13'), Q = 2k sin(θ/2), isotropic
# illumination, Gaussian speeds with rms ⟨V²⟩^½ = v (mm/s; Eqs. 12, 12'); m̄ mean number of
# collisions with moving cells (Poisson). Times in s, angular frequencies in rad/s.


def bonner_i1(tau: FloatArray | float, a: float, v: float, k: float, xi: float = 0.1) -> FloatArray:
    """Single-shift field correlation I₁(τ), bonner1981 Eq. 18 (Eq. 2 with Eqs. 12', 13'):
    I₁ = 2ξ/(2ξ + T²) [1 − e^(−2ξL²) e^(−L²T²)]/[1 − e^(−2ξL²)], L = 2ka, T = vτ/(√6 a) (Eq. 17);
    k (1/mm) in the medium."""
    big_l = 2 * k * a
    t2 = (v * np.asarray(tau, dtype=np.float64) / (np.sqrt(6) * a)) ** 2
    e = np.exp(-2 * xi * big_l**2)
    return np.asarray(2 * xi / (2 * xi + t2) * (1 - e * np.exp(-(big_l**2) * t2)) / (1 - e))


def bonner_i1_large(tau: FloatArray | float, a: float, v: float, xi: float = 0.1) -> FloatArray:
    """I₁ = 2ξ/(2ξ + T²) for large particles, 2ξL² ≫ 1 (bonner1981 Eq. 20)."""
    t2 = (v * np.asarray(tau, dtype=np.float64) / (np.sqrt(6) * a)) ** 2
    return np.asarray(2 * xi / (2 * xi + t2))


def bonner_g2(
    tau: FloatArray | float, mbar: float, a: float, v: float, xi: float = 0.1, beta: float = 1.0
) -> FloatArray:
    """g⁽²⁾(τ) = 1 + β[exp(2m̄(I₁ − 1)) − exp(−2m̄)] with I₁ from Eq. 20 (bonner1981 Eq. 24): the
    fluctuating part, the unshifted fraction e^(−m̄) beating with itself removed."""
    i1 = bonner_i1_large(tau, a, v, xi)
    return np.asarray(1 + beta * (np.exp(2 * mbar * (i1 - 1)) - np.exp(-2 * mbar)))


def bonner_s_j(j: int, w: float) -> float:
    """S_j(W) = (1/π) ∫₀^∞ cos(WZ) (1 + Z²)^(−j) dZ (bonner1981 Eq. 30 with Eq. 20) for j = 1–3.
    Eq. 31 prints S₁ = e^(−W); Eq. 30 gives ½ e^(−W), consistent with S₂, S₃ and Eq. A3."""
    if j == 1:
        return float(0.5 * np.exp(-w))
    if j == 2:
        return float(0.25 * (1 + w) * np.exp(-w))
    if j == 3:
        return float((w * w + 3 * w + 3) * np.exp(-w) / 16)
    raise ValueError("closed forms for j = 1, 2, 3 only")


def bonner_mean_frequency(
    mbar: float, a: float, v: float, xi: float = 0.1, beta: float = 1.0, j_max: int = 400
) -> float:
    """First moment ⟨ω⟩ = ∫ |ω| S(ω) dω (rad/s) of the spectrum of Eq. 24, bonner1981 Eq. 36/A4:
    ⟨ω⟩ = v/((12ξ)^½ a) β f(m̄), f = (2e^(−2m̄)/√π) Σ_(j≥1) (2m̄)^j Γ(j + ½)/(Γ(j + 1)Γ(j))."""
    from scipy.special import gammaln

    j = np.arange(1, j_max + 1)
    terms = np.exp(j * np.log(2 * mbar) + gammaln(j + 0.5) - gammaln(j + 1) - gammaln(j) - 2 * mbar)
    f = 2 / np.sqrt(np.pi) * float(np.sum(terms))
    return float(v / (np.sqrt(12 * xi) * a) * beta * f)
