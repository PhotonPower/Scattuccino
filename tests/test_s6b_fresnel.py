"""S6b: polarised Fresnel events at index steps (T7c)."""

from __future__ import annotations

import numpy as np

from scattuccino.polarisation import (
    apply_amplitudes,
    fresnel_amplitudes,
    interface,
    mie_mueller,
    rayleigh_mueller,
    simulate_polarised,
)


def test_fresnel_rhomb_turns_linear_into_circular() -> None:
    # two total internal reflections at 54.6° in glass n = 1.51: 45° retardation each
    n, th = 1.51, np.radians(54.6)
    u = np.array([[np.sin(th), 0.0, np.cos(th)]])
    v = np.array([[np.cos(th), 0.0, -np.sin(th)]])  # p̂ of the plane of incidence
    s = np.array([[1.0, 0.0, 1.0, 0.0]])  # linear at 45° between p and s
    u, v, s, r1 = interface(u, v, s, n, 1.0, np.array([[0.0, 0.0, 1.0]]), np.array([0.5]))
    u, v, s, r2 = interface(u, v, s, n, 1.0, np.array([[0.0, 0.0, -1.0]]), np.array([0.5]))
    assert r1[0] and r2[0]  # total reflection
    assert abs(abs(s[0, 3]) - 1) < 1e-6 and abs(s[0, 2]) < 1e-3


def test_brewster_reflection_is_s_polarised() -> None:
    res = simulate_polarised(
        1.0,
        0.0,
        1.0,
        rayleigh_mueller(),
        np.array([1.0, 0, 0, 0]),
        10,
        seed=1,
        n_slab=1.33,
        cos_incident=float(np.cos(np.arctan(1.33))),
    )
    assert res.specular is not None
    assert np.allclose(res.specular[1], [1.0, -1.0, 0.0, 0.0], atol=1e-12)


def test_slab_without_scattering_against_incoherent_sum() -> None:
    n, th, mua, d = 1.5, np.radians(50), 0.3, 1.0
    s_in = np.array([[1.0, 0.0, 1.0, 0.0]])
    res = simulate_polarised(
        d,
        mua,
        0.0,
        rayleigh_mueller(),
        s_in[0],
        200_000,
        seed=1,
        n_slab=n,
        cos_incident=float(np.cos(th)),
    )
    ci = np.array([np.cos(th)])
    rp1, rs1, tp1, ts1, ct = fresnel_amplitudes(1.0, n, ci)
    rp2, rs2, tp2, ts2, _ = fresnel_amplitudes(n, 1.0, np.real(ct))
    a = np.exp(-mua * d / np.real(ct[0]))
    t_ref, r_ref = np.zeros(4), np.zeros(4)
    for k in range(60):  # incoherent sum over bounces; each path keeps its p–s coherence
        t_ref += apply_amplitudes(tp1 * (rp2 * rp2) ** k * tp2, ts1 * (rs2 * rs2) ** k * ts2, s_in)[
            0
        ] * a ** (2 * k + 1)
        r_ref += apply_amplitudes(
            tp1 * rp2 ** (2 * k + 1) * tp2, ts1 * rs2 ** (2 * k + 1) * ts2, s_in
        )[0] * a ** (2 * k + 2)
    for side, ref in ((1, t_ref), (-1, r_ref)):
        sel = res.side == side
        contrib = (res.weight[sel] * res.packets)[:, None] * res.stokes[sel]
        full = np.zeros((res.packets, 4))
        full[: contrib.shape[0]] = contrib
        mean = full.mean(axis=0)
        err = full.std(axis=0) / np.sqrt(res.packets)
        assert np.all(np.abs(mean - ref) < 4 * err + 1e-12)
    spec, s_spec = res.specular  # type: ignore[misc]
    s_r = apply_amplitudes(rp1, rs1, s_in)[0]
    assert abs(spec - s_r[0]) < 1e-15 and np.allclose(s_spec, s_r / s_r[0], atol=1e-15)


def test_energy_with_scattering_and_index_steps() -> None:
    res = simulate_polarised(
        1.0,
        0.0,
        5.0,
        mie_mueller(1.59 / 1.33 + 0j, 2.0),
        np.array([1.0, 0.3, 0.5, 0.6]),
        20_000,
        seed=2,
        n_slab=1.4,
        n_below=1.5,
        cos_incident=float(np.cos(0.6)),
    )
    assert res.specular is not None
    assert abs(res.weight.sum() + res.specular[0] - 1) < 1e-12  # no absorption: everything leaves
    assert np.all(np.sum(res.stokes[:, 1:] ** 2, axis=1) <= 1 + 1e-12)
