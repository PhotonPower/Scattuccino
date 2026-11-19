"""Extended Huygens-Fresnel model of the OCT signal: closed form (S9c-1, D11 A) and general form
with the phase function of the medium (S9c-2, D11 B).

Thrane, Yura and Andersen (thrane2000): mean square heterodyne signal from a discontinuity at
depth z (diffuse or specular reflection) behind tissue with small-angle forward scattering,
⟨i²(z)⟩ = ⟨i²⟩₀ Ψ(z). The tissue is described by μs, the rms scattering angle θ_rms (levitz2004:
e⁻¹ half width of a Gaussian fit to the forward lobe, g_eff = cos θ_rms) and the mean index n;
absorption and bulk backscattering are neglected. A lens of focal length f at distance d from the
surface focuses a Gaussian beam (1/e intensity radius w0 in the lens plane); the mixing plane is
the lens plane. Units: lengths in mm, wavelength in µm (vacuum), angles in rad.

The cross term of Ψ is 4 e^(−μs z)(1 − e^(−μs z))/(1 + w_S²/w_H²) as obtained by integrating
thrane2000 Eq. 21 (and printed in levitz2004 Eq. 2); thrane2000 Eq. 22 prints 2 (issue #36).
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Literal

import numpy as np
from scipy.special import j0

from scattuccino._types import FloatArray

Reflection = Literal["diffuse", "specular"]


def lateral_coherence_length(
    z: FloatArray | float,
    mus: float,
    theta_rms: float,
    n: float,
    d: FloatArray | float,
    wavelength_um: float,
    reflection: Reflection = "diffuse",
    shower_curtain: bool = True,
) -> FloatArray:
    """ρ₀(z) in the lens plane (mm).

    Diffuse reflection: spherical wave from a point in the discontinuity plane,
    ρ₀ = √(3/(μs z)) λ/(π θ_rms) (n d + z)/z (thrane2000 Eqs. A6, A7, 24); the factor (n d + z)/z
    is the shower-curtain effect, `shower_curtain=False` drops it (thrane2000 Fig. 2, curve 2).
    Specular reflection: plane wave, ρ₀ = √(1/(2μs z)) λ/(π θ_rms) (Eq. 26)."""
    lam = wavelength_um * 1e-3
    z = np.asarray(z, dtype=np.float64)
    if reflection == "specular":
        return np.asarray(np.sqrt(1 / (2 * mus * z)) * lam / (np.pi * theta_rms))
    rho0 = np.sqrt(3 / (mus * z)) * lam / (np.pi * theta_rms)
    if shower_curtain:
        rho0 = rho0 * (n * np.asarray(d) + z) / z
    return np.asarray(rho0)


def beam_radii(
    z: FloatArray | float,
    mus: float,
    theta_rms: float,
    n: float,
    w0: float,
    f: float,
    d: FloatArray | float,
    wavelength_um: float,
    reflection: Reflection = "diffuse",
    shower_curtain: bool = True,
) -> tuple[FloatArray, FloatArray]:
    """(w_H, w_S): 1/e intensity radii in the discontinuity plane without and with scattering
    (thrane2000 Eqs. 13, 14), ray-matrix elements A = 1, B = d + z/n, k = 2π/λ (vacuum)."""
    k = 2 * np.pi / (wavelength_um * 1e-3)
    b = np.asarray(d) + np.asarray(z) / n
    w_h2 = w0**2 * (1 - b / f) ** 2 + (b / (k * w0)) ** 2
    rho0 = lateral_coherence_length(
        z, mus, theta_rms, n, d, wavelength_um, reflection, shower_curtain
    )
    w_s2 = w_h2 + (2 * b / (k * rho0)) ** 2
    return np.asarray(np.sqrt(w_h2)), np.asarray(np.sqrt(w_s2))


def heterodyne_efficiency(
    z: FloatArray | float,
    mus: float,
    theta_rms: float,
    n: float,
    w0: float,
    f: float,
    d: FloatArray | float,
    wavelength_um: float,
    reflection: Reflection = "diffuse",
    shower_curtain: bool = True,
) -> FloatArray:
    """Heterodyne efficiency factor Ψ(z) = ⟨i²(z)⟩/⟨i²⟩₀ (dimensionless).

    Diffuse: Ψ = e^(−2μs z) + 4 e^(−μs z)(1 − e^(−μs z))/(1 + w_S²/w_H²)
    + (1 − e^(−μs z))² w_H²/w_S² (thrane2000 Eq. 21 integrated; levitz2004 Eq. 2).
    Specular: Ψ = e^(−2μs z) + (1 − e^(−2μs z)) w_H²/w_S² (thrane2000 Eq. 25).
    For the focus on the discontinuity pass d = f − z/n (Eq. 23)."""
    w_h, w_s = beam_radii(z, mus, theta_rms, n, w0, f, d, wavelength_um, reflection, shower_curtain)
    q = (w_h / w_s) ** 2
    e = np.exp(-mus * np.asarray(z, dtype=np.float64))
    if reflection == "specular":
        return np.asarray(e * e + (1 - e * e) * q)
    return np.asarray(e * e + 4 * e * (1 - e) / (1 + 1 / q) + (1 - e) ** 2 * q)


# S9c-2 (D11 B): general form with the phase function of the medium.
#
# The mean irradiance of a beam U(p) = exp(−a p²) launched at the surface (index-matched medium,
# wavenumber k in the medium) at depth z is lutomirski1971 Eq. 10, the Fourier transform of
# M_s(ϱ, z) W(ϱ) with W(ϱ) = ∫ U(r + ϱ/2) U*(r − ϱ/2) e^(ikϱ·r/z) d²r ∝ exp(−ϱ²/(2A²)). The
# heterodyne signal from a δ-correlated reflector at z is the overlap of the mean irradiances of
# the illumination and of the detection mode (thrane2000 Eqs. 11, 15, 16; the same beam by
# reciprocity), by Parseval ∝ ∫ M_s² |W|² d²ϱ; Ψ normalises it to M_s = 1.


def overlap_width(z: float, a: complex, k: float) -> float:
    """A (mm) of |W(ϱ)|² ∝ exp(−ϱ²/A²) for U = exp(−a p²) at the surface (Re a > 0):
    1/A² = α + (k/z − 2β)²/(4α), a = α + iβ (Gaussian integral over r)."""
    alpha, beta = a.real, a.imag
    return float(1 / np.sqrt(alpha + (k / z - 2 * beta) ** 2 / (4 * alpha)))


def phase_correlation(
    phase: object, k: float, rho_max: float, n_rho: int = 4001, n_theta: int = 20_001
) -> Callable[[FloatArray], FloatArray]:
    """b(ρ): small-angle phase correlation of a medium with phase function p (per sr) and
    wavenumber k (1/mm) in the medium, b(ρ) = 2π ∫₀^(π/2) p(cos θ) J₀(kρ sin θ) sin θ dθ –
    lutomirski1971 Eq. 11 with the scattered wave vector k sin θ in place of the spectrum (the
    forward half only: backscattered light leaves the beam, b(0) = forward fraction). Tabulated on
    0 ≤ ρ ≤ rho_max (mm), zero beyond."""
    theta = np.linspace(0.0, np.pi / 2, n_theta)
    w = phase.pdf(np.cos(theta)) * np.sin(theta)  # type: ignore[attr-defined]
    w = w * np.gradient(theta) * 2 * np.pi
    rho = np.linspace(0.0, rho_max, n_rho)
    b = np.empty(n_rho)
    for i0 in range(0, n_rho, 64):
        arg = k * rho[i0 : i0 + 64, None] * np.sin(theta)[None, :]
        b[i0 : i0 + 64] = j0(arg) @ w

    def table(r: FloatArray) -> FloatArray:
        return np.asarray(np.interp(r, rho, b, right=0.0))

    return table


def spherical_mtf(
    b: Callable[[FloatArray], FloatArray], mus: float, z: float, n_u: int = 200
) -> Callable[[FloatArray], FloatArray]:
    """M_s(ϱ, z) = exp{−μs z [1 − ∫₀¹ b(ϱu) du]} (lutomirski1971 Eq. 11 for spherical waves; the
    scatterer at depth s sees the separation ϱ(z − s)/z; μs z is the phase variance, thrane2000
    Eq. 18). Gauss-Legendre in u."""
    x, w = np.polynomial.legendre.leggauss(n_u)
    u, wq = 0.5 * (x + 1), 0.5 * w

    def mtf(rho: FloatArray) -> FloatArray:
        r = np.asarray(rho, dtype=np.float64)
        mean_b = b((r[..., None] * u).reshape(-1)).reshape(r.shape + u.shape) @ wq
        return np.asarray(np.exp(-mus * z * (1 - mean_b)))

    return mtf


def general_efficiency(
    z: float,
    a: complex,
    k: float,
    mtf: Callable[[FloatArray], FloatArray],
    n_rho: int = 20_000,
) -> float:
    """Ψ(z) = ∫ M_s(ϱ)² e^(−ϱ²/A²) d²ϱ/(πA²) (dimensionless), A from `overlap_width`. With the
    Gaussian approximation of M_s (thrane2000 Eq. 20) this is the closed form (Eq. 21 integrated)
    for an aperture at the surface. Quadrature on a logarithmic grid in ϱ (10⁻⁷ A … 10 A)."""
    big_a = overlap_width(z, a, k)
    rho = np.geomspace(1e-7 * big_a, 10 * big_a, n_rho).astype(np.float64)
    m = mtf(rho)
    f = m * m * 2 * (rho / big_a) ** 2 * np.exp(-((rho / big_a) ** 2))
    head = float(m[0] ** 2 * (rho[0] / big_a) ** 2)  # 0 ≤ ϱ < ϱ_min
    return head + float(np.trapezoid(f, np.log(rho)))
