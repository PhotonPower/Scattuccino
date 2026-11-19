"""S9c-1: closed-form extended Huygens-Fresnel model of the OCT signal (T10b, thrane2000).

The closed form Eq. 22 against the quadrature of Eq. 21 (cross term 4, levitz2004 Eq. 2), the
lateral coherence length against the quadrature of Eq. A1, limits, and the curves of thrane2000
Figs. 2 and 3 digitised from the paper.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad

from scattuccino.ehf import Reflection, beam_radii, heterodyne_efficiency, lateral_coherence_length

LAM = 0.814  # µm, thrane2000 Figs. 2 and 3


def test_eq22_against_quadrature_of_eq21() -> None:
    # Eq. 21: ⟨i²⟩ ∝ ∫ [a e^(−r²/w_H²)/w_H² + b e^(−r²/w_S²)/w_S²]² d²r with a = e^(−μs z),
    # b = 1 − a; Eq. 22 divides by ⟨i²⟩₀, i.e. Ψ = (2 w_H²/π) × the integral
    z, mus, th, n, w0, f = 0.05, 4.0, 0.1, 1.4, 0.5, 5.0
    d = f - z / n  # focus on the discontinuity: the halo is comparable to w_H
    w_h, w_s = beam_radii(z, mus, th, n, w0, f, d, LAM)
    a = np.exp(-mus * z)
    b = 1 - a

    def integrand(r: float) -> float:
        h = a * np.exp(-(r * r) / w_h**2) / w_h**2 + b * np.exp(-(r * r) / w_s**2) / w_s**2
        return float(2 * np.pi * r * h * h)

    ref = 2 * w_h**2 / np.pi * quad(integrand, 0, 30 * w_s, limit=500)[0]
    assert 0.05 < w_h**2 / w_s**2 < 0.95  # all three terms matter
    assert abs(heterodyne_efficiency(z, mus, th, n, w0, f, d, LAM) / ref - 1) < 1e-8


def test_coherence_length_against_quadrature_of_a1() -> None:
    # Eq. A1: Γ_pt(ρ) = exp{−ρ²/4 ∫₀ᴸ k² θ² μs (B(s)/B(L))² ds}, B(s) = s/n in the tissue
    # (0 ≤ s ≤ z), no scattering between the surface and the lens; Γ = exp(−ρ²/ρ₀²)
    z, mus, th, n, d = 0.4, 10.0, 0.3, 1.4, 3.0
    k = 2 * np.pi / (LAM * 1e-3)
    b_l = d + z / n
    integral = quad(lambda s: k**2 * th**2 * mus * (s / n / b_l) ** 2, 0, z)[0]
    rho0 = np.sqrt(4 / integral)
    assert abs(lateral_coherence_length(z, mus, th, n, d, LAM) / rho0 - 1) < 1e-10


def test_limits() -> None:
    z, n, w0, f, d = 0.5, 1.4, 0.5, 5.0, 5.0 - 0.5 / 1.4
    refls: tuple[Reflection, ...] = ("diffuse", "specular")
    for refl in refls:
        # no scattering: the signal of the focused beam without loss
        psi = heterodyne_efficiency(z, 1e-9, 0.3, n, w0, f, d, LAM, reflection=refl)
        assert abs(psi - 1) < 1e-6
        # scattered light without lateral coherence (θ_rms → ∞): only the unscattered round trip
        mus = 4.0
        psi = heterodyne_efficiency(z, mus, 1e6, n, w0, f, d, LAM, reflection=refl)
        assert abs(psi / np.exp(-2 * mus * z) - 1) < 1e-6


def test_thrane_fig2() -> None:
    # log10 Ψ digitised from thrane2000 Fig. 2 (300 dpi, plot frame as axes; reading ≈ 0.1
    # decades): curve 1 diffuse with shower curtain, 2 diffuse without, 3 specular; focus on the
    # discontinuity (d = f − z/n); λ = 814 nm, μs = 20/mm, θ_rms = 0.3, n = 1.4, f = 5 mm,
    # w0 = 0.5 mm
    fig = {
        0.5: (-4.50, -6.77, -7.49),
        1.0: (-5.40, -7.05, -7.81),
        1.5: (-5.93, -7.23, -7.98),
        1.75: (-6.13, -7.29, -8.03),
    }
    n, f = 1.4, 5.0
    for z, curves in fig.items():
        d = f - z / n
        cases = (
            {"reflection": "diffuse"},
            {"reflection": "diffuse", "shower_curtain": False},
            {"reflection": "specular"},
        )
        for ref, kw in zip(curves, cases, strict=True):
            psi = heterodyne_efficiency(z, 20.0, 0.3, n, 0.5, f, d, LAM, **kw)  # type: ignore[arg-type]
            assert abs(np.log10(psi) - ref) < 0.25, (z, kw)


def test_thrane_fig3() -> None:
    # log10 Ψ against μs digitised from thrane2000 Fig. 3 (diffuse, specular): z = 0.5 mm,
    # θ_rms = 0.38, n = 1.33, f = 16 mm, w0 = 0.125 mm, glass plate δ = 1.256 mm, n_g = 1.5 in
    # d(z) = f − z/n + δ(1 − 1/n_g) (Eq. 27)
    fig = {2.0: (-0.88, -0.99), 4.5: (-1.61, -2.02), 6.5: (-1.91, -2.84), 7.5: (-1.99, -3.25)}
    z, n, f = 0.5, 1.33, 16.0
    d = f - z / n + 1.256 * (1 - 1 / 1.5)
    for mus, (diffuse, specular) in fig.items():
        cases: tuple[tuple[float, Reflection], ...] = ((diffuse, "diffuse"), (specular, "specular"))
        for ref, refl in cases:
            psi = heterodyne_efficiency(z, mus, 0.38, n, 0.125, f, d, LAM, reflection=refl)
            assert abs(np.log10(psi) - ref) < 0.25, (mus, refl)
