"""Dynamic light scattering in tissue (stage S8a): DCS, laser speckle contrast, laser Doppler.

The photon protocol carries the momentum transfer per layer, Y_t = Σ (1 − cos θ_j) over the
scattering events in layer t (`Protocol.momentum`). With q_j² = 2k²(1 − cos θ_j), k = 2πn/λ, and
Gaussian displacements of the moving scatterers, a path contributes
exp(−(1/3) k_t² α_t Y_t ⟨Δr_t²(τ)⟩) to the field autocorrelation, so

    g1(τ) = Σ_paths w · exp(−Σ_t (1/3) k_t² α_t Y_t ⟨Δr_t²(τ)⟩) / Σ_paths w

(Boas & Yodh 1997; α_t: fraction of scattering events from moving scatterers, used as a mean
factor). Motion models: Brownian ⟨Δr²⟩ = 6 D_B τ, random flow ⟨Δr²⟩ = v² τ².

g2(τ) = 1 + β |g1(τ)|² (Siegert), speckle contrast for exposure T:
K² = (2β/T) ∫₀ᵀ (1 − τ/T) |g1(τ)|² dτ (Bandyopadhyay et al. 2005).

Units: mm, µm (wavelength), s; D_B in mm²/s, v in mm/s.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from scattuccino._types import FloatArray
from scattuccino.transport import Protocol


@dataclass(frozen=True)
class Motion:
    """Motion of the moving scatterers in one layer: alpha (fraction of events), Brownian
    diffusion coefficient d_b (mm²/s) and/or random flow speed v (mm/s)."""

    alpha: float = 0.0
    d_b: float = 0.0
    v: float = 0.0

    def msd(self, tau: FloatArray) -> FloatArray:
        """⟨Δr²(τ)⟩ in mm²."""
        return np.asarray(6 * self.d_b * tau + (self.v * tau) ** 2)


def g1(
    p: Protocol,
    motions: tuple[Motion, ...],
    wavelength_um: float,
    tau: FloatArray,
    rho_lo: float,
    rho_hi: float,
    side: int = -1,
) -> FloatArray:
    """Normalised field autocorrelation of the light leaving through `side` in the annulus."""
    if p.momentum is None:
        raise ValueError("the protocol carries no momentum transfer")
    sel = (p.side == side) & (p.rho >= rho_lo) & (p.rho < rho_hi)
    w = p.weight[sel]
    y = p.momentum[sel]
    k2 = (2 * np.pi * p.n_layer / (wavelength_um * 1e-3)) ** 2  # 1/mm²
    tau = np.atleast_1d(np.asarray(tau, dtype=np.float64))
    out = np.empty(tau.size)
    for i, t in enumerate(tau):
        expo = np.zeros(w.size)
        for layer, m in enumerate(motions):
            if m.alpha > 0:
                expo += k2[layer] / 3 * m.alpha * y[:, layer] * float(m.msd(np.array([t]))[0])
        out[i] = np.sum(w * np.exp(-expo)) / np.sum(w)
    return out


def g2(g1_values: FloatArray, beta: float = 1.0) -> FloatArray:
    return np.asarray(1 + beta * np.abs(g1_values) ** 2)


def speckle_contrast(
    tau: FloatArray, g1_values: FloatArray, exposure: float, beta: float = 1.0
) -> float:
    """K for exposure T from g1 on a grid tau covering [0, T] (trapezoidal rule)."""
    sel = tau <= exposure * (1 + 1e-12)
    t, g = tau[sel], np.abs(g1_values[sel]) ** 2
    return float(np.sqrt(2 * beta / exposure * np.trapezoid((1 - t / exposure) * g, t)))


def speckle_contrast_single_exponential(tau_c: float, exposure: float, beta: float = 1.0) -> float:
    """K for g1 = exp(−τ/τ_c): K² = β[τ_c/T + (τ_c²/2T²)(e^(−2T/τ_c) − 1)] (Bandyopadhyay 2005)."""
    x = exposure / tau_c
    return float(np.sqrt(beta * (1 / x + (np.exp(-2 * x) - 1) / (2 * x * x))))


# S8b: laser Doppler with tagged moving scatterers (`Layer.moving_fraction`). A path with momentum
# transfer Y_t over its moving events in layer t contributes exp(−Σ_t (1/3) k_t² Y_t ⟨Δr_t²(τ)⟩)
# exactly (Gaussian displacements); paths without moving events keep their phase (unshifted
# light), so g1 → the unshifted fraction for τ → ∞ (bonner1981 Eqs. 21–23).


def g1_tagged(
    p: Protocol,
    motions: tuple[Motion, ...],
    wavelength_um: float,
    tau: FloatArray,
    rho_lo: float,
    rho_hi: float,
    side: int = -1,
) -> FloatArray:
    """Normalised field autocorrelation from the tagged moving events (`Motion.alpha` is not
    used: the fraction is in the transport)."""
    if p.momentum_moving is None:
        raise ValueError("the protocol carries no tagged moving events (Layer.moving_fraction)")
    sel = (p.side == side) & (p.rho >= rho_lo) & (p.rho < rho_hi)
    w = p.weight[sel]
    y = p.momentum_moving[sel]
    k2 = (2 * np.pi * p.n_layer / (wavelength_um * 1e-3)) ** 2
    tau = np.atleast_1d(np.asarray(tau, dtype=np.float64))
    out = np.empty(tau.size)
    for i, t in enumerate(tau):
        expo = np.zeros(w.size)
        for layer, m in enumerate(motions):
            expo += k2[layer] / 3 * y[:, layer] * float(m.msd(np.array([t]))[0])
        out[i] = np.sum(w * np.exp(-expo)) / np.sum(w)
    return out


def unshifted_fraction(p: Protocol, rho_lo: float, rho_hi: float, side: int = -1) -> float:
    """Weight fraction of the detected light without any moving event (limit of g1 for τ → ∞)."""
    if p.n_moving is None:
        raise ValueError("the protocol carries no tagged moving events (Layer.moving_fraction)")
    sel = (p.side == side) & (p.rho >= rho_lo) & (p.rho < rho_hi)
    w = p.weight[sel]
    return float(np.sum(w[p.n_moving[sel].sum(axis=1) == 0]) / np.sum(w))


def doppler_spectrum(tau: FloatArray, g: FloatArray, omega: FloatArray) -> FloatArray:
    """S(ω) = (1/π) ∫₀^∞ cos(ωτ) G(τ) dτ for the fluctuating part G of the intensity correlation,
    e.g. G = g2 − 1 − β f₀² with the unshifted fraction f₀ (bonner1981 Eq. 26); two-sided,
    ∫ S dω = G(0). Trapezoidal rule on `tau` (starting at 0, G negligible at its end)."""
    tau = np.asarray(tau, dtype=np.float64)
    g = np.asarray(g, dtype=np.float64)
    omega = np.atleast_1d(np.asarray(omega, dtype=np.float64))
    out = np.empty(omega.size)
    for i0 in range(0, omega.size, 256):
        c = np.cos(omega[i0 : i0 + 256, None] * tau[None, :])
        out[i0 : i0 + 256] = np.trapezoid(c * g[None, :], tau, axis=1) / np.pi
    return out


def mean_frequency(tau: FloatArray, g: FloatArray) -> float:
    """First moment ⟨ω⟩ = ∫ |ω| S(ω) dω (rad/s, bonner1981 Eq. 32) directly from G:
    ⟨ω⟩ = (2/π) ∫₀^∞ (G(0) − G(τ))/τ² dτ (G even and smooth at 0; the part beyond the grid as
    G(0)/τ_max with G → 0)."""
    tau = np.asarray(tau, dtype=np.float64)
    g = np.asarray(g, dtype=np.float64)
    t, d = tau[1:], g[0] - g[1:]
    head = float(d[0] / t[0] ** 2 * t[0])  # (G(0) − G)/τ² is finite at 0: −G''(0)/2
    return float(2 / np.pi * (head + np.trapezoid(d / t**2, t) + g[0] / tau[-1]))
