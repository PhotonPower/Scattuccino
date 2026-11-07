"""Adding-doubling for a homogeneous slab (stage S1b), independent of the Monte Carlo.

Matrix-operator form of the adding-doubling method (van de Hulst; Prahl, van Gemert, Welch 1993):
radiance at quadrature cosines μ_i ∈ (0, 1] for up- and downward directions, operators that map
incident to reflected/transmitted radiance with the quadrature weights built in,

    R̃_ij = 2π a Δτ p̄(−μ_i, μ_j) w_j / μ_i,   T̃_ij = δ_ij e^(−Δτ/μ_i) + 2π a Δτ p̄(μ_i, μ_j) w_j / μ_i

for a thin start layer Δτ = τ/2ⁿ (first order in Δτ), doubled n times. A collimated beam at μ₀
is carried separately (direct attenuation e^(−τ/μ₀) and its diffuse sources r₀, t₀). Boundaries
with an index step reflect diffuse light by Fresnel (diagonal, TIR beyond the critical cosine; the
quadrature is split there) and the collimated beam by its Fresnel factor; the coupled system of the
slab and both boundaries is solved directly. p̄ is the azimuthal average of the Henyey-Greenstein
phase function (normalised ∫p dΩ = 1), renormalised on the quadrature so that scattering conserves
energy exactly.

Returns the same quantities as the Monte Carlo of `scattuccino.layers.run`: specular reflection at
the first surface, all other reflected and transmitted power, absorption.
"""

from __future__ import annotations

import itertools

import numpy as np

from scattuccino._types import FloatArray
from scattuccino.layers import fresnel_unpolarised


def _quadrature(n: int, breaks: list[float]) -> tuple[FloatArray, FloatArray]:
    """Gauss-Legendre on (0, 1] split at the given cosines (critical angles)."""
    edges = sorted({0.0, 1.0, *[b for b in breaks if 0.0 < b < 1.0]})
    xs, ws = [], []
    per = max(n // (len(edges) - 1), 4)
    x, w = np.polynomial.legendre.leggauss(per)
    for a, b in itertools.pairwise(edges):
        xs.append(a + (b - a) * 0.5 * (x + 1))
        ws.append((b - a) * 0.5 * w)
    return np.concatenate(xs), np.concatenate(ws)


def _hg_azimuthal(g: float, mu_out: FloatArray, mu_in: FloatArray, n_phi: int = 1024) -> FloatArray:
    """p̄(μ, μ′) = (1/2π) ∫ p_HG(μμ′ + √(1−μ²)√(1−μ′²) cos φ) dφ (rows μ, columns μ′)."""
    phi = 2 * np.pi * (np.arange(n_phi) + 0.5) / n_phi
    s_out = np.sqrt(np.clip(1 - mu_out**2, 0, 1))[:, None, None]
    s_in = np.sqrt(np.clip(1 - mu_in**2, 0, 1))[None, :, None]
    cos_t = mu_out[:, None, None] * mu_in[None, :, None] + s_out * s_in * np.cos(phi)[None, None, :]
    p = (1 - g * g) / (4 * np.pi * (1 + g * g - 2 * g * cos_t) ** 1.5)
    return np.asarray(p.mean(axis=2))


def slab(
    d: float,
    mua: float,
    mus: float,
    g: float,
    n_slab: float,
    n_above: float = 1.0,
    n_below: float = 1.0,
    cos_incident: float = 1.0,
    n_quad: int = 48,
    doublings: int = 24,
) -> dict[str, float]:
    """Specular and diffuse reflectance, transmittance, absorption (d = np.inf: half space)."""
    mut = mua + mus
    a = mus / mut if mut > 0 else 0.0
    breaks = []
    for n_out in (n_above, n_below):
        if n_out < n_slab:
            breaks.append(float(np.sqrt(1 - (n_out / n_slab) ** 2)))
    mu, w = _quadrature(n_quad, breaks)
    r_sp, mu0_arr = fresnel_unpolarised(n_above, n_slab, np.array([cos_incident]))
    mu0 = float(mu0_arr[0])
    power = 1.0 - float(r_sp[0])
    # phase function on the nodes, plus the collimated direction
    p_up = _hg_azimuthal(g, -mu, mu)  # into upward directions from downward ones
    p_dn = _hg_azimuthal(g, mu, mu)
    norm = 2 * np.pi * ((p_up + p_dn) * w[:, None]).sum(axis=0)  # per incident direction
    p_up, p_dn = p_up / norm, p_dn / norm
    p0_up = _hg_azimuthal(g, -mu, np.array([mu0]))[:, 0]
    p0_dn = _hg_azimuthal(g, mu, np.array([mu0]))[:, 0]
    norm0 = 2 * np.pi * ((p0_up + p0_dn) * w).sum()
    p0_up, p0_dn = p0_up / norm0, p0_dn / norm0
    tau = mut * d
    infinite = not np.isfinite(tau)
    n_dbl = doublings if not infinite else 60
    dtau = (tau if not infinite else 1.0) / 2.0**doublings if not infinite else 1e-6
    big_r = 2 * np.pi * a * dtau * p_up * w[None, :] / mu[:, None]
    big_t = np.diag(np.exp(-dtau / mu)) + 2 * np.pi * a * dtau * p_dn * w[None, :] / mu[:, None]
    r0 = a * dtau * p0_up / (mu * mu0)
    t0 = a * dtau * p0_dn / (mu * mu0)
    e0 = np.exp(-dtau / mu0)
    eye = np.eye(len(mu))
    for _ in range(n_dbl):
        inv = np.linalg.inv(eye - big_r @ big_r)
        u = inv @ (big_r @ t0 + e0 * r0)  # upward at the middle (B reflects)
        dvec = np.linalg.solve(eye - big_r @ big_r, t0 + e0 * (big_r @ r0))
        r0, t0 = r0 + big_t @ u, big_t @ dvec + e0 * t0
        big_r, big_t = big_r + big_t @ inv @ big_r @ big_t, big_t @ inv @ big_t
        e0 = e0 * e0
        if infinite and e0 < 1e-300 and np.max(np.abs(big_t)) < 1e-14:
            break
    # boundaries: diagonal Fresnel reflection of diffuse light, collimated Fresnel factors
    f_top, _ = fresnel_unpolarised(n_slab, n_above, mu)
    f_bot, _ = fresnel_unpolarised(n_slab, n_below, mu)
    rt0 = float(fresnel_unpolarised(n_slab, n_above, np.array([mu0]))[0][0])
    rb0 = float(fresnel_unpolarised(n_slab, n_below, np.array([mu0]))[0][0])
    q = e0 * e0 * rb0 * rt0
    pd = power / (1 - q)  # all downward collimated passes
    pu = power * e0 * rb0 / (1 - q)  # all upward collimated passes
    n = len(mu)
    ft, fb = np.diag(f_top), np.diag(f_bot)
    # unknowns U_top (up, at the top) and D_bot (down, at the bottom)
    a_mat = np.block([[eye - big_r @ ft, -big_t @ fb], [-big_t @ ft, eye - big_r @ fb]])
    rhs = np.concatenate([r0 * pd + t0 * pu, t0 * pd + r0 * pu])
    sol = np.linalg.solve(a_mat, rhs)
    u_top, d_bot = sol[:n], sol[n:]
    flux = 2 * np.pi * w * mu
    # escaping diffuse light, plus collimated light reflected at the bottom and leaving the top
    r_diff = float(np.sum(flux * (1 - f_top) * u_top)) + pu * e0 * (1 - rt0)
    t_diff = float(np.sum(flux * (1 - f_bot) * d_bot)) + pd * e0 * (1 - rb0)
    total_in = 1.0
    absorbed = total_in - float(r_sp[0]) - r_diff - t_diff
    return {
        "r_specular": float(r_sp[0]),
        "r_diffuse": r_diff,
        "t_diffuse": t_diff,
        "absorbed": absorbed,
    }
