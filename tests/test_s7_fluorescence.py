"""S7a: fluorescence and Raman (T8).

1. Non-scattering medium (both wavelengths), matched index: escaping fluorescence
   F = Q ∫ μa,f e^(−μa,x z) · ½ E₂(μa,m z) dz (exact), absorbed fraction μa,f/μa,x.
2. Energy: emitted = escaping + re-absorbed.
3. Forward against adjoint: the two-stage result equals Q ∫ A_f(z) E_m(z) dz with the escape
   function E_m(z) from separate isotropic point sources at fixed depths (independent runs).
4. Spatially offset detection (SORS): Raman-active only in a deep layer – its share of the
   detected Raman light grows with the source–detector offset relative to a superficial layer.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad
from scipy.special import expn

from scattuccino import detectors as det
from scattuccino.fluorescence import simulate_fluorescence
from scattuccino.layers import Layer, Stack
from scattuccino.transport import Sources, simulate


def test_non_scattering_against_closed_form() -> None:
    mua_t, mua_f, mua_m, q = 0.5, 0.3, 0.4, 0.2
    sx = Stack((Layer(40.0, mua_t + mua_f, 0.0, 0.0, 1.0),))
    sm = Stack((Layer(40.0, mua_m, 0.0, 0.0, 1.0),))
    res = simulate_fluorescence(sx, sm, np.array([mua_f]), q, 200_000, 200_000, seed=3)
    mux = mua_t + mua_f
    assert abs(res.absorbed_fluorophore - mua_f / mux) < 1e-12  # all light absorbed (thick slab)
    ref = (
        q
        * quad(lambda z: mua_f * np.exp(-mux * z) * 0.5 * expn(2, mua_m * z), 0, np.inf, limit=200)[
            0
        ]
    )
    p = res.emission.protocol
    top = p.weight[p.side == -1].sum()
    frac = ref / res.emitted
    sigma = res.emitted * np.sqrt(frac * (1 - frac) / 200_000)
    assert abs(top - ref) < 4 * sigma


def test_energy_of_the_emission() -> None:
    sx = Stack((Layer(2.0, 0.15, 10.0, 0.9, 1.4),))
    sm = Stack((Layer(2.0, 0.06, 9.0, 0.9, 1.4),))
    res = simulate_fluorescence(sx, sm, np.array([0.05]), 0.1, 50_000, 50_000, seed=4)
    p = res.emission.protocol
    assert abs(p.weight.sum() + res.emission.absorbed.sum() - res.emitted) < 1e-6 * res.emitted


def test_forward_against_adjoint() -> None:
    d, dz = 2.0, 0.1
    sx = Stack((Layer(d, 0.15, 10.0, 0.9, 1.4),))
    sm = Stack((Layer(d, 0.06, 9.0, 0.9, 1.4),))
    q, mua_f = 0.1, 0.05
    fwd = simulate_fluorescence(sx, sm, np.array([mua_f]), q, 100_000, 100_000, seed=5)
    p = fwd.emission.protocol
    f_forward = p.weight[p.side == -1].sum()
    # adjoint: fluorophore absorption per depth bin × escape probability of an isotropic source
    ev = fwd.excitation.events
    assert ev is not None
    edges = np.arange(0.0, d + dz / 2, dz)
    a_f, _ = np.histogram(ev["z"], bins=edges, weights=ev["dw"] * mua_f / 0.15)
    n = 6_000
    escape = []
    for k, z0 in enumerate(0.5 * (edges[1:] + edges[:-1])):
        src = Sources(
            np.zeros(n),
            np.zeros(n),
            np.full(n, z0),
            np.zeros(n, dtype=np.int64),
            np.full(n, 1.0 / n),
        )
        e = simulate(sm, n, seed=100 + k, sources=src).protocol
        escape.append(e.weight[e.side == -1].sum())
    f_adjoint = q * np.sum(a_f * np.array(escape))
    assert abs(f_forward / f_adjoint - 1) < 0.03  # two Monte Carlo estimates (≈1 % each)


def test_spatially_offset_raman() -> None:
    # Raman-active top layer (0.5 mm) or bottom layer (below), same optics otherwise
    def run(active_top: bool) -> tuple[float, float]:
        mua_r = (0.02, 0.0) if active_top else (0.0, 0.02)
        sx = Stack(
            (
                Layer(0.5, 0.05 + mua_r[0], 10.0, 0.9, 1.4),
                Layer(5.0, 0.05 + mua_r[1], 10.0, 0.9, 1.4),
            )
        )
        sm = Stack((Layer(0.5, 0.05, 10.0, 0.9, 1.4), Layer(5.0, 0.05, 10.0, 0.9, 1.4)))
        res = simulate_fluorescence(sx, sm, np.array(mua_r), 1.0, 40_000, 40_000, seed=6)
        r, _ = det.reflectance_vs_rho(res.emission.protocol, np.array([0.0, 0.5, 4.0, 6.0]))
        return float(r[0]), float(r[2])

    top_near, top_far = run(True)
    deep_near, deep_far = run(False)
    # the deep layer's share of the Raman signal is larger at a 4–6 mm offset than at 0–0.5 mm
    assert deep_far / top_far > 2 * deep_near / top_near
