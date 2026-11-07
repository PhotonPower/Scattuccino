"""Phase functions (stage S3): normalised ∫p dΩ = 1, sampling of cos θ, Legendre moments.

- `HenyeyGreenstein(g)`;
- `ModifiedHG(g_hg, beta)`: β·HG + (1 − β)·(3/4π) cos²θ (Bevilacqua & Depeursinge 1999), g = β g_hg;
- `Gegenbauer(g_gk, alpha)`: p ∝ (1 + g² − 2g cos θ)^(−(α + 1)) (Reynolds & McCormick 1980),
  α = ½ is HG; sampled by the inverse distribution in closed form;
- `Tabulated(mu, p)`: any p(cos θ) on a grid (e.g. Mie), sampled by the inverse distribution.

Moments g_n = ⟨P_n(cos θ)⟩; the subdiffusive parameters γ = (1 − g₂)/(1 − g₁) and
δ = (1 − g₃)/(1 − g₁) (Bevilacqua & Depeursinge 1999) describe reflectance near the source.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.polynomial import legendre

from scattuccino._types import FloatArray

_MU, _W = legendre.leggauss(2000)


class PhaseFunction:
    """Base class: pdf(μ) per steradian, sample(ξ) → μ."""

    def pdf(self, mu: FloatArray) -> FloatArray:
        raise NotImplementedError

    def sample(self, xi: FloatArray) -> FloatArray:
        raise NotImplementedError

    def moment(self, n: int) -> float:
        """g_n = 2π ∫ p(μ) P_n(μ) dμ."""
        coef = np.zeros(n + 1)
        coef[n] = 1.0
        return float(2 * np.pi * np.sum(_W * self.pdf(_MU) * legendre.legval(_MU, coef)))

    @property
    def g(self) -> float:
        return self.moment(1)

    @property
    def gamma(self) -> float:
        return (1 - self.moment(2)) / (1 - self.moment(1))

    @property
    def delta(self) -> float:
        return (1 - self.moment(3)) / (1 - self.moment(1))


@dataclass(frozen=True)
class HenyeyGreenstein(PhaseFunction):
    g_hg: float

    def pdf(self, mu: FloatArray) -> FloatArray:
        g = self.g_hg
        return np.asarray((1 - g * g) / (4 * np.pi * (1 + g * g - 2 * g * mu) ** 1.5))

    def sample(self, xi: FloatArray) -> FloatArray:
        g = self.g_hg
        if abs(g) < 1e-12:
            return np.asarray(2 * xi - 1)
        t = (1 - g * g) / (1 - g + 2 * g * xi)
        return np.asarray(np.clip((1 + g * g - t * t) / (2 * g), -1.0, 1.0))


@dataclass(frozen=True)
class ModifiedHG(PhaseFunction):
    g_hg: float
    beta: float

    def pdf(self, mu: FloatArray) -> FloatArray:
        return np.asarray(
            self.beta * HenyeyGreenstein(self.g_hg).pdf(mu)
            + (1 - self.beta) * 3 / (4 * np.pi) * mu**2
        )

    def sample(self, xi: FloatArray) -> FloatArray:
        # one uniform number: choose the branch, rescale it within the branch
        hg = xi < self.beta
        x_hg = np.where(hg, xi / max(self.beta, 1e-300), 0.5)
        x_c2 = np.where(hg, 0.5, (xi - self.beta) / max(1 - self.beta, 1e-300))
        return np.asarray(
            np.where(hg, HenyeyGreenstein(self.g_hg).sample(x_hg), np.cbrt(2 * x_c2 - 1))
        )


@dataclass(frozen=True)
class Gegenbauer(PhaseFunction):
    g_gk: float
    alpha: float

    def _norm(self) -> float:
        g, a = self.g_gk, self.alpha
        return float((1 - g) ** (-2 * a) - (1 + g) ** (-2 * a))

    def pdf(self, mu: FloatArray) -> FloatArray:
        g, a = self.g_gk, self.alpha
        u = 1 + g * g - 2 * g * mu
        return np.asarray(2 * g * a * u ** (-(a + 1)) / self._norm() / (2 * np.pi))

    def sample(self, xi: FloatArray) -> FloatArray:
        g, a = self.g_gk, self.alpha
        u = ((1 + g) ** (-2 * a) + xi * self._norm()) ** (-1 / a)
        return np.asarray(np.clip((1 + g * g - u) / (2 * g), -1.0, 1.0))


class Tabulated(PhaseFunction):
    """p on an ascending grid of cos θ (linear interpolation), normalised on construction."""

    def __init__(self, mu: FloatArray, p: FloatArray) -> None:
        mu = np.asarray(mu, dtype=np.float64)
        p = np.asarray(p, dtype=np.float64)
        if mu[0] > -1 + 1e-12 or mu[-1] < 1 - 1e-12 or np.any(np.diff(mu) <= 0):
            raise ValueError("grid must be ascending and cover [-1, 1]")
        cdf = np.concatenate([[0.0], np.cumsum(0.5 * (p[1:] + p[:-1]) * np.diff(mu))])
        self.mu, self.p = mu, p * (1.0 / (2 * np.pi * cdf[-1]))
        self.cdf = cdf / cdf[-1]

    def pdf(self, mu: FloatArray) -> FloatArray:
        return np.asarray(np.interp(mu, self.mu, self.p))

    def sample(self, xi: FloatArray) -> FloatArray:
        # exact inverse of the piecewise-linear pdf within each interval
        k = np.clip(np.searchsorted(self.cdf, xi, side="right") - 1, 0, len(self.mu) - 2)
        m0, m1 = self.mu[k], self.mu[k + 1]
        p0, p1 = self.p[k], self.p[k + 1]
        slope = (p1 - p0) / (m1 - m0)
        # 2π ∫_{m0}^{μ} (p0 + slope (t − m0)) dt = ξ − cdf_k (p normalised with ∫ p dΩ = 1)
        c = xi - self.cdf[k]
        lin = np.abs(slope) < 1e-12 * np.maximum(p0, 1e-300)
        with np.errstate(divide="ignore", invalid="ignore"):
            dq = np.where(
                lin,
                c / (2 * np.pi * p0),
                (-p0 + np.sqrt(np.maximum(p0 * p0 + 2 * slope * c / (2 * np.pi), 0))) / slope,
            )
        return np.asarray(np.clip(m0 + dq, -1.0, 1.0))


def mie_phase(m: complex, x: float, n_mu: int = 4001) -> tuple[Tabulated, float, float]:
    """Unpolarised Mie phase function of a sphere (relative index m, size parameter x = 2πr n_med/λ)
    as `Tabulated`, with Q_sca and Q_ext (Bohren & Huffman 1983, ch. 4). The grid is dense near the
    forward direction."""
    from scattuccino.mie import qeff, s12

    t = np.linspace(0, 1, n_mu)
    theta = np.pi * t**2  # dense at θ = 0
    mu = np.unique(np.cos(theta))  # ascending, without coincident nodes near θ = 0
    s1, s2 = s12(m, x, np.arccos(mu))  # noqa: RUF100
    p = 0.5 * (np.abs(s1) ** 2 + np.abs(s2) ** 2)
    q_ext, q_sca = qeff(m, x)
    return Tabulated(mu, p), q_sca, q_ext
