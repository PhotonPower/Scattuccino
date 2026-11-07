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

Stage S6a: layers with matched refractive indices (no Fresnel at the boundaries); Fresnel Mueller
matrices at index steps follow in S6b.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from scattuccino._types import FloatArray
from scattuccino.mie import s12
from scattuccino.phase import Tabulated


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
    rng: np.random.Generator,
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """One scattering for N packets: directions u (N,3), frames v (N,3), Stokes (N,4) with I = 1."""
    mu = table.phase.sample(xi_theta)
    m12, m33, m34 = table.ratios(mu)
    q, uu = stokes[:, 1], stokes[:, 2]
    bound = 1.0 + np.abs(m12) * np.sqrt(q * q + uu * uu)
    phi = np.empty(len(mu))
    todo = np.arange(len(mu))
    while todo.size:  # rejection: p(φ) ∝ 1 + m12 (Q cos 2φ + U sin 2φ)
        trial = 2 * np.pi * rng.random(todo.size)
        f = 1.0 + m12[todo] * (q[todo] * np.cos(2 * trial) + uu[todo] * np.sin(2 * trial))
        ok = rng.random(todo.size) * bound[todo] <= f
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


@dataclass
class PolarisedResult:
    side: np.ndarray  # type: ignore[type-arg]
    u: FloatArray  # exit directions (N, 3)
    v: FloatArray  # exit reference vectors (N, 3)
    stokes: FloatArray  # (N, 4), I = 1 (weight separately)
    weight: FloatArray  # per incident power
    scatterings: np.ndarray  # type: ignore[type-arg]
    packets: int


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
) -> PolarisedResult:
    """Slab 0 < z < d (np.inf: half space), index-matched, normal incidence along +z with the
    Stokes vector `stokes_in` in the frame v = x̂ (w = ŷ)."""
    rng = np.random.default_rng(seed)
    mut = mua + mus
    s_in = np.asarray(stokes_in, dtype=np.float64) / stokes_in[0]
    u = np.tile([0.0, 0.0, 1.0], (packets, 1))
    v = np.tile([1.0, 0.0, 0.0], (packets, 1))
    s = np.tile(s_in, (packets, 1))
    z = np.zeros(packets)
    w = np.ones(packets)
    alive = np.ones(packets, dtype=bool)
    nsc = np.zeros(packets, dtype=np.int64)
    side = np.zeros(packets, dtype=np.int64)
    for _ in range(max_events):
        idx = np.nonzero(alive)[0]
        if idx.size == 0:
            break
        step = -np.log(rng.random(idx.size)) / mut
        z_new = z[idx] + u[idx, 2] * step
        out_top, out_bot = z_new < 0, z_new > d
        side[idx[out_top]], side[idx[out_bot]] = -1, 1
        alive[idx[out_top | out_bot]] = False
        stay = idx[~(out_top | out_bot)]
        z[stay] = z_new[~(out_top | out_bot)]
        w[stay] *= mus / mut
        u[stay], v[stay], s[stay] = scatter_step(
            table, u[stay], v[stay], s[stay], rng.random(stay.size), rng
        )
        nsc[stay] += 1
        low = stay[w[stay] < w_threshold]
        if low.size:
            survive = rng.random(low.size) < 1.0 / roulette_m
            w[low[survive]] *= roulette_m
            w[low[~survive]] = 0.0
            alive[low[~survive]] = False
            side[low[~survive]] = 2
    leave = (side == -1) | (side == 1)
    return PolarisedResult(
        side[leave], u[leave], v[leave], s[leave], w[leave] / packets, nsc[leave], packets
    )


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
