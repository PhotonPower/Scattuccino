"""Polarised photon transport (stage S6a): Stokes vectors and Mueller matrices.

Conventions (docs/konventionen.md): time factor e^(−iωt), Stokes vector (I, Q, U, V) with respect
to a reference frame (v, w) perpendicular to the propagation direction u (w = u × v):
Q = |E_v|² − |E_w|², U = 2 Re(E_v E_w*), V = −2 Im(E_v E_w*) (Bohren & Huffman 1983, §2.11), so
that e^(−iωt) fields with E_w = i E_v (… V = +1 …) follow the book's sign of circular polarisation.

Scattering (Bohren & Huffman ch. 3/4): in the frame whose v lies in the scattering plane
(parallel) the Mueller matrix of a sphere is

    [[S11, S12, 0, 0], [S12, S11, 0, 0], [0, 0, S33, S34], [0, 0, −S34, S33]]

with S11 = (|S2|² + |S1|²)/2, S12 = (|S2|² − |S1|²)/2, S33 = Re(S2 S1*), S34 = Im(S2 S1*).

Algorithm (frame vectors, as Ramella-Roman, Prahl, Jacques 2005, method with the reference vector
carried along): θ from S11 (the φ-averaged phase function), φ by rejection from
S11 + S12 (Q cos 2φ + U sin 2φ)/I, then rotate the frame by φ about u (Stokes rotation
Q′ = Q cos 2φ + U sin 2φ, U′ = −Q sin 2φ + U cos 2φ), scatter by θ in the plane (u, v₁):
u′ = cos θ u + sin θ v₁, v′ = −sin θ u + cos θ v₁, and S′ = M(θ) S′ normalised to I = 1 (the
weight is unchanged because θ, φ are sampled from the exact distribution).

Boundaries (S6b): plane interfaces with the polarised Fresnel event `interface` (amplitudes
r_p, r_s, t_p, t_s in the plane-of-incidence frame p̂ = ŝ × û, including total internal
reflection with its phase).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from scattuccino._types import FloatArray, IntArray
from scattuccino.mie import s12
from scattuccino.phase import Tabulated
from scattuccino.random import CounterRandom


@dataclass
class MuellerTable:
    """Mueller ratios over cos θ and the phase function (S11 normalised) for sampling θ."""

    mu: FloatArray
    m12: FloatArray  # S12/S11
    m33: FloatArray  # S33/S11
    m34: FloatArray  # S34/S11
    phase: Tabulated
    pure: bool = False  # from amplitudes S1, S2: m12² + m33² + m34² = 1 (pure states stay pure)

    def ratios(self, mu: FloatArray) -> tuple[FloatArray, FloatArray, FloatArray]:
        m12 = np.asarray(np.interp(mu, self.mu, self.m12))
        m33 = np.asarray(np.interp(mu, self.mu, self.m33))
        m34 = np.asarray(np.interp(mu, self.mu, self.m34))
        if self.pure:  # restore the identity that linear interpolation breaks
            scale = np.sqrt(
                np.clip(1 - m12 * m12, 0, None) / np.maximum(m33 * m33 + m34 * m34, 1e-300)
            )
            m33, m34 = m33 * scale, m34 * scale
        return m12, m33, m34


def mueller_from_amplitudes(mu: FloatArray, s1: np.ndarray, s2: np.ndarray) -> MuellerTable:  # type: ignore[type-arg]
    s11 = 0.5 * (np.abs(s2) ** 2 + np.abs(s1) ** 2)
    return MuellerTable(
        mu,
        0.5 * (np.abs(s2) ** 2 - np.abs(s1) ** 2) / s11,
        np.real(s2 * np.conj(s1)) / s11,
        np.imag(s2 * np.conj(s1)) / s11,
        Tabulated(mu, s11),
        pure=True,
    )


def mie_mueller(m: complex, x: float, n_mu: int = 4001) -> MuellerTable:
    t = np.linspace(0, 1, n_mu)
    mu = np.unique(np.cos(np.pi * t**2))
    s1, s2 = s12(m, x, np.arccos(mu))
    return mueller_from_amplitudes(mu.astype(np.float64), s1, s2)


def rayleigh_mueller(n_mu: int = 2001) -> MuellerTable:
    mu = np.linspace(-1.0, 1.0, n_mu).astype(np.float64)
    return mueller_from_amplitudes(mu, np.ones_like(mu, dtype=complex), mu.astype(complex))


def rotate_stokes(q: FloatArray, u: FloatArray, phi: FloatArray) -> tuple[FloatArray, FloatArray]:
    c, s = np.cos(2 * phi), np.sin(2 * phi)
    return np.asarray(q * c + u * s), np.asarray(-q * s + u * c)


def scatter_step(
    table: MuellerTable,
    u: FloatArray,
    v: FloatArray,
    stokes: FloatArray,
    xi_theta: FloatArray,
    rng: CounterRandom,
    ids: IntArray,
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """One scattering for N packets: directions u (N,3), frames v (N,3), Stokes (N,4) with I = 1."""
    mu = table.phase.sample(xi_theta)
    m12, m33, m34 = table.ratios(mu)
    q, uu = stokes[:, 1], stokes[:, 2]
    bound = 1.0 + np.abs(m12) * np.sqrt(q * q + uu * uu)
    phi = np.empty(len(mu))
    todo = np.arange(len(mu))
    while todo.size:  # rejection: p(φ) ∝ 1 + m12 (Q cos 2φ + U sin 2φ)
        trial = 2 * np.pi * rng.uniform(ids[todo])
        f = 1.0 + m12[todo] * (q[todo] * np.cos(2 * trial) + uu[todo] * np.sin(2 * trial))
        ok = rng.uniform(ids[todo]) * bound[todo] <= f
        phi[todo[ok]] = trial[ok]
        todo = todo[~ok]
    return apply_scattering(table, u, v, stokes, mu, phi)


def apply_scattering(
    table: MuellerTable,
    u: FloatArray,
    v: FloatArray,
    stokes: FloatArray,
    mu: FloatArray,
    phi: FloatArray,
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Deterministic part of a scattering: rotate the frame by φ about u, scatter by θ = acos μ in
    the plane (u, v₁), apply the Mueller matrix, normalise I to 1."""
    q, uu = stokes[:, 1], stokes[:, 2]
    m12, m33, m34 = table.ratios(mu)
    w = np.cross(u, v)
    cp, sp = np.cos(phi)[:, None], np.sin(phi)[:, None]
    v1 = cp * v + sp * w
    q1, u1 = rotate_stokes(q, uu, phi)
    st, ct = np.sqrt(np.clip(1 - mu * mu, 0, 1))[:, None], mu[:, None]
    u_new = ct * u + st * v1
    v_new = -st * u + ct * v1
    i2 = 1.0 + m12 * q1
    q2 = m12 + q1
    u2 = m33 * u1 + m34 * stokes[:, 3]
    v2 = -m34 * u1 + m33 * stokes[:, 3]
    s_new = np.column_stack([np.ones_like(i2), q2 / i2, u2 / i2, v2 / i2])
    # renormalise the frame (rounding)
    u_new /= np.linalg.norm(u_new, axis=1)[:, None]
    v_new -= np.sum(v_new * u_new, axis=1)[:, None] * u_new
    v_new /= np.linalg.norm(v_new, axis=1)[:, None]
    return u_new, v_new, s_new


def fresnel_amplitudes(
    n1: float, n2: float, cos_i: FloatArray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:  # type: ignore[type-arg]
    """Complex Fresnel amplitudes (r_p, r_s, t_p, t_s) and cos θ_t for fields in the bases
    p̂ = ŝ × û, ŝ = (û × n̂)/|û × n̂| (e^(−iωt); beyond the critical angle cos θ_t = +i√…, the
    evanescent wave decays). With these bases r_p = r_s at normal incidence up to the sign of the
    reflected p̂ (p̂_r = −p̂_i), as in Raytatouille's Convention A."""
    cos_i = np.abs(np.asarray(cos_i, dtype=np.float64))
    sin2_t = (n1 / n2) ** 2 * (1 - cos_i**2)
    cos_t = np.where(
        sin2_t <= 1,
        np.sqrt(np.clip(1 - sin2_t, 0, None)) + 0j,
        1j * np.sqrt(np.clip(sin2_t - 1, 0, None)),
    )
    r_s = (n1 * cos_i - n2 * cos_t) / (n1 * cos_i + n2 * cos_t)
    r_p = (n2 * cos_i - n1 * cos_t) / (n2 * cos_i + n1 * cos_t)
    t_s = 2 * n1 * cos_i / (n1 * cos_i + n2 * cos_t)
    t_p = 2 * n1 * cos_i / (n2 * cos_i + n1 * cos_t)
    return r_p, r_s, t_p, t_s, cos_t


def apply_amplitudes(a_p: np.ndarray, a_s: np.ndarray, stokes: FloatArray) -> FloatArray:  # type: ignore[type-arg]
    """Stokes (in the frame v = p̂, w = ŝ) after E_p → a_p E_p, E_s → a_s E_s; not normalised."""
    i, q, u, v = stokes[:, 0], stokes[:, 1], stokes[:, 2], stokes[:, 3]
    pp, ss = np.abs(a_p) ** 2, np.abs(a_s) ** 2
    c = a_p * np.conj(a_s)
    return np.column_stack(
        [
            0.5 * (pp + ss) * i + 0.5 * (pp - ss) * q,
            0.5 * (pp - ss) * i + 0.5 * (pp + ss) * q,
            c.real * u + c.imag * v,
            -c.imag * u + c.real * v,
        ]
    )


def frame_to_plane(
    u: FloatArray, v: FloatArray, stokes: FloatArray, normal: FloatArray
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Rotate the Stokes frame to (p̂, ŝ) of the plane of incidence with the surface normal.
    Returns (p̂, ŝ, rotated Stokes). Along the normal any ŝ ⊥ û is chosen (the current w)."""
    w = np.cross(u, v)
    s_hat = np.cross(u, normal)
    norm = np.linalg.norm(s_hat, axis=1)
    along = norm < 1e-12
    s_hat = np.where(along[:, None], w, s_hat / np.where(along, 1.0, norm)[:, None])
    p_hat = np.cross(s_hat, u)
    phi = np.arctan2(np.sum(w * p_hat, axis=1), np.sum(v * p_hat, axis=1))
    q, uu = rotate_stokes(stokes[:, 1], stokes[:, 2], phi)
    out = np.column_stack([stokes[:, 0], q, uu, stokes[:, 3]])
    return (
        np.asarray(p_hat, dtype=np.float64),
        np.asarray(s_hat, dtype=np.float64),
        np.asarray(out, dtype=np.float64),
    )


def interface(
    u: FloatArray,
    v: FloatArray,
    stokes: FloatArray,
    n1: float,
    n2: float,
    normal: FloatArray,
    xi: FloatArray,
) -> tuple[FloatArray, FloatArray, FloatArray, np.ndarray]:  # type: ignore[type-arg]
    """Polarised Fresnel event at a plane surface (normal pointing from medium 1 into medium 2):
    the reflectance R(S) = ½(|r_p|² + |r_s|²) + ½(|r_p|² − |r_s|²) Q/I (in the p/s frame) decides
    reflection (ξ < R) or transmission; the Stokes vector follows the Fresnel amplitudes and is
    normalised to I = 1. Returns new u, v (= p̂ of the new direction), Stokes, reflected mask."""
    if n1 == n2:
        return u, v, stokes, np.zeros(len(u), dtype=bool)
    p_hat, s_hat, st = frame_to_plane(u, v, stokes, normal)
    cos_i = np.abs(np.sum(u * normal, axis=1))
    r_p, r_s, t_p, t_s, cos_t = fresnel_amplitudes(n1, n2, cos_i)
    big_r = (
        0.5 * (np.abs(r_p) ** 2 + np.abs(r_s) ** 2)
        + 0.5 * (np.abs(r_p) ** 2 - np.abs(r_s) ** 2) * st[:, 1]
    )
    refl = xi < big_r
    u_new = u.copy()
    u_new[refl] = u[refl] - 2 * (np.sum(u[refl] * normal[refl], axis=1))[:, None] * normal[refl]
    tr = ~refl
    ratio = n1 / n2
    ct = np.real(cos_t[tr])
    u_new[tr] = ratio * u[tr] + (ct - ratio * cos_i[tr])[:, None] * normal[tr]
    s_new = np.where(refl[:, None], apply_amplitudes(r_p, r_s, st), apply_amplitudes(t_p, t_s, st))
    s_new = s_new / s_new[:, 0:1]
    u_new /= np.linalg.norm(u_new, axis=1)[:, None]
    v_new = np.cross(s_hat, u_new)
    return (
        np.asarray(u_new),
        np.asarray(v_new, dtype=np.float64),
        np.asarray(s_new, dtype=np.float64),
        refl,
    )


@dataclass
class PolarisedResult:
    side: np.ndarray  # type: ignore[type-arg]
    u: FloatArray  # exit directions (N, 3)
    v: FloatArray  # exit reference vectors (N, 3)
    stokes: FloatArray  # (N, 4), I = 1 (weight separately)
    weight: FloatArray  # per incident power
    scatterings: np.ndarray  # type: ignore[type-arg]
    packets: int
    specular: tuple[float, FloatArray] | None = None  # (reflectance, Stokes) at entry


def simulate_polarised(
    d: float,
    mua: float,
    mus: float,
    table: MuellerTable,
    stokes_in: FloatArray,
    packets: int,
    seed: int = 1,
    w_threshold: float = 1e-4,
    roulette_m: float = 10.0,
    max_events: int = 100_000,
    n_slab: float = 1.0,
    n_above: float = 1.0,
    n_below: float = 1.0,
    cos_incident: float = 1.0,
) -> PolarisedResult:
    """Slab 0 < z < d (np.inf: half space) of index n_slab between n_above and n_below; incidence
    from above in the x–z plane at θ = acos(cos_incident) with the Stokes vector `stokes_in` in the
    frame v = p̂ (in the plane of incidence; at normal incidence v = x̂), w = ŝ. Fresnel events with
    the polarised reflectance (S6b); the specular reflection at entry is returned separately in
    `PolarisedResult.specular` (Stokes in the reflected p/s frame, weight = reflectance)."""
    rng = CounterRandom(seed, packets)
    mut = mua + mus
    s_in = np.asarray(stokes_in, dtype=np.float64) / stokes_in[0]
    sin_i = np.sqrt(1 - cos_incident**2)
    u0 = np.array([[sin_i, 0.0, cos_incident]])
    # p̂ = ŝ × û with ŝ = (û × ẑ)/|…| = −ŷ for incidence in the x–z plane; at normal incidence x̂
    v0 = np.array([[cos_incident, 0.0, -sin_i]]) if sin_i > 0 else np.array([[1.0, 0.0, 0.0]])
    zhat = np.array([[0.0, 0.0, 1.0]])
    # entry: deterministic split into specular reflection and transmitted part
    p_hat, s_hat, st = frame_to_plane(u0, v0, s_in[None, :], zhat)
    r_p, r_s, t_p, t_s, cos_t = fresnel_amplitudes(n_above, n_slab, np.array([cos_incident]))
    s_r = apply_amplitudes(r_p, r_s, st)[0]
    s_t = apply_amplitudes(t_p, t_s, st)[0]
    big_r = float(s_r[0])  # reflectance for I = 1
    ratio = n_above / n_slab
    ut = ratio * u0[0] + (float(np.real(cos_t[0])) - ratio * cos_incident) * zhat[0]
    ut /= np.linalg.norm(ut)
    vt = np.cross(s_hat[0], ut)
    specular = (big_r, s_r / s_r[0] if big_r > 0 else s_r)
    u = np.tile(ut, (packets, 1))
    v = np.tile(vt, (packets, 1))
    s = np.tile(s_t / s_t[0], (packets, 1))
    z = np.zeros(packets)
    w = np.full(packets, 1.0 - big_r)
    alive = np.ones(packets, dtype=bool)
    nsc = np.zeros(packets, dtype=np.int64)
    side = np.zeros(packets, dtype=np.int64)
    s_left = np.zeros(packets)
    for _ in range(max_events):
        idx = np.nonzero(alive)[0]
        if idx.size == 0:
            break
        new = s_left[idx] <= 0
        s_left[idx[new]] = -np.log(rng.uniform(idx[new]))
        step = s_left[idx] / mut if mus > 0 else np.full(idx.size, np.inf)
        uz = u[idx, 2]
        with np.errstate(divide="ignore", invalid="ignore"):
            d_b = np.where(uz < 0, -z[idx] / uz, np.where(uz > 0, (d - z[idx]) / uz, np.inf))
        hit = d_b <= step
        dist = np.where(hit, d_b, step)
        z[idx] = z[idx] + uz * dist
        if mua > 0:
            w[idx] *= np.exp(-mua * dist) if mus == 0 else 1.0  # without scattering: Beer-Lambert
        s_left[idx] = np.where(hit, s_left[idx] - d_b * mut, 0.0)
        hb = idx[hit]
        if hb.size:
            up = u[hb, 2] < 0
            for sel, n_out, nrm in ((up, n_above, -1.0), (~up, n_below, 1.0)):
                ids = hb[sel]
                if ids.size == 0:
                    continue
                normal = np.tile([0.0, 0.0, nrm], (ids.size, 1))
                un, vn, sn, refl = interface(
                    u[ids],
                    v[ids],
                    s[ids],
                    n_slab,
                    n_out,
                    normal,
                    rng.uniform(ids) if n_slab != n_out else np.ones(ids.size),
                )
                u[ids], v[ids], s[ids] = un, vn, sn
                out = ids[~refl]
                side[out] = -1 if nrm < 0 else 1
                alive[out] = False
        ib = idx[~hit]
        if ib.size and mus > 0:
            w[ib] *= mus / mut
            u[ib], v[ib], s[ib] = scatter_step(table, u[ib], v[ib], s[ib], rng.uniform(ib), rng, ib)
            nsc[ib] += 1
            low = ib[w[ib] < w_threshold]
            if low.size:
                survive = rng.uniform(low) < 1.0 / roulette_m
                w[low[survive]] *= roulette_m
                w[low[~survive]] = 0.0
                alive[low[~survive]] = False
                side[low[~survive]] = 2
        elif ib.size:  # no scattering: absorption only (handled along the path)
            pass
    leave = (side == -1) | (side == 1)
    res = PolarisedResult(
        side[leave], u[leave], v[leave], s[leave], w[leave] / packets, nsc[leave], packets
    )
    res.specular = specular
    return res


def stokes_lab(res: PolarisedResult) -> FloatArray:
    """Stokes vectors of the exits in the lab frame: reference = projection of x̂ onto the plane
    ⊥ u (fine for exits not along ±x)."""
    u, v = res.u, res.v
    ref = np.array([1.0, 0.0, 0.0])[None, :] - u[:, 0:1] * u
    ref /= np.linalg.norm(ref, axis=1)[:, None]
    w = np.cross(u, v)
    phi = np.arctan2(np.sum(ref * w, axis=1), np.sum(ref * v, axis=1))
    q, uu = rotate_stokes(res.stokes[:, 1], res.stokes[:, 2], phi)
    return np.column_stack([res.stokes[:, 0], q, uu, res.stokes[:, 3]])
