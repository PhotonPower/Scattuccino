"""S9c-2: general extended Huygens-Fresnel model with the phase function of the medium (T10b).

Ψ(z) = ∫ M_s(ϱ, z)² e^(−ϱ²/A²) d²ϱ/(πA²): the signal is the overlap of the mean irradiances of the
illumination and of the (reciprocal) detection mode, each from lutomirski1971 Eq. 10, by
Parseval; M_s is the spherical-wave MTF (lutomirski1971 Eq. 11, J₀ → ∫₀¹ J₀(·u) du). With the
Gaussian approximation of the MTF (thrane2000 Eq. 20) it must reproduce the closed form (S9c-1).
"""

from __future__ import annotations

import numpy as np
from scipy.special import erf

from scattuccino._types import FloatArray
from scattuccino.ehf import (
    beam_radii,
    general_efficiency,
    heterodyne_efficiency,
    lateral_coherence_length,
    overlap_width,
    phase_correlation,
    spherical_mtf,
)
from scattuccino.phase import Tabulated

LAM = 0.814


def test_overlap_width_against_thrane_beam_radius() -> None:
    # aperture at the surface (d = 0, n = 1): w_H² = 2 z²/(k² A²) with w_H from thrane2000 Eq. 13
    k = 2 * np.pi / (LAM * 1e-3)
    w0, f = 0.05, 1.0
    a = 1 / (2 * w0**2) + 1j * k / (2 * f)  # thrane2000 Eq. 6: U = exp(−p²(1/w0² + ik/f)/2)
    for z in (0.2, 0.9, 1.0, 1.4):
        w_h, _ = beam_radii(z, 1.0, 0.3, 1.0, w0, f, 0.0, LAM)
        big_a = overlap_width(z, a, k)
        assert abs(2 * z**2 / (k * big_a) ** 2 / w_h**2 - 1) < 1e-12


def test_gaussian_mtf_reproduces_the_closed_form() -> None:
    # thrane2000 Eq. 20 with ρ₀ from S9c-1 (n = 1, aperture at the surface, d = 0)
    k = 2 * np.pi / (LAM * 1e-3)
    w0, f, mus, th = 0.05, 1.0, 8.0, 0.1
    a = 1 / (2 * w0**2) + 1j * k / (2 * f)
    for z in (0.1, 0.3, 0.9):
        rho0 = float(lateral_coherence_length(z, mus, th, 1.0, 0.0, LAM))
        e = np.exp(-mus * z)

        def mtf(rho: FloatArray, e: float = e, rho0: float = rho0) -> FloatArray:
            return np.asarray(e + (1 - e) * np.exp(-(rho**2) / rho0**2))

        psi = general_efficiency(z, a, k, mtf)
        ref = heterodyne_efficiency(z, mus, th, 1.0, w0, f, 0.0, LAM)
        assert abs(psi / ref - 1) < 1e-6, z


def test_phase_correlation_of_a_narrow_gaussian_lobe() -> None:
    # p ∝ exp(−θ²/θ_r²) with θ_r ≪ 1: b(ρ) = exp(−k² θ_r² ρ²/4) (small-angle Fourier pair)
    k, th = 10.0, 0.05
    mu = np.cos(np.linspace(np.pi, 0, 200_001))
    phase = Tabulated(mu, np.exp(-(np.arccos(np.clip(mu, -1, 1)) ** 2) / th**2))
    rho = np.array([0.0, 1.0, 2.0, 4.0, 6.0])
    b = phase_correlation(phase, k, rho_max=8.0)(rho)
    assert np.allclose(b, np.exp(-(k**2) * th**2 * rho**2 / 4), atol=2e-3)


def test_spherical_mtf_of_a_gaussian_correlation() -> None:
    # M_s = exp{−μs z [1 − ∫₀¹ b(ρu) du]}, b = exp(−ρ²/ρ_φ²): ∫₀¹ = √π erf(c)/(2c), c = ρ/ρ_φ
    mus, z, rho_phi = 5.0, 0.4, 0.01
    rho = np.linspace(1e-6, 0.05, 50).astype(np.float64)
    m = spherical_mtf(lambda r: np.exp(-(r**2) / rho_phi**2), mus, z)(rho)
    c = rho / rho_phi
    ref = np.exp(-mus * z * (1 - np.sqrt(np.pi) * erf(c) / (2 * c)))
    assert np.allclose(m, ref, rtol=1e-6)


def test_no_scattering() -> None:
    k = 2 * np.pi / (LAM * 1e-3)
    a = 1 / (2 * 0.05**2) + 1j * k / 2
    assert abs(general_efficiency(0.5, a, k, lambda r: np.ones_like(r)) - 1) < 1e-6
