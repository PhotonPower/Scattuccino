"""Extended Huygens-Fresnel model of the OCT signal, closed form (stage S9c-1, D11 A).

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

from typing import Literal

import numpy as np

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
