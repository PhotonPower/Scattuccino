"""S2: photon protocol, detectors and volume scoring (T3).

Exact checks (bookkeeping) and comparisons with diffusion theory in its range of validity
(several transport lengths from source and boundary; approximations, documented tolerances).
"""

from __future__ import annotations

import numpy as np
import pytest

from scattuccino import detectors as det
from scattuccino import diffusion as dif
from scattuccino.adding_doubling import slab
from scattuccino.layers import Layer, Stack
from scattuccino.transport import Grid, Simulation, reweight, simulate

HALF = Stack((Layer(np.inf, 0.05, 10.0, 0.9, 1.0),), 1.0, 1.0)


@pytest.fixture(scope="module")
def half_space() -> Simulation:
    return simulate(HALF, 60_000, seed=5, grid=Grid(0.25, 40, 0.25, 60))


def test_bookkeeping(half_space: Simulation) -> None:
    p = half_space.protocol
    # scored absorption equals the absorption per layer; energy balance
    assert abs(half_space.absorbed_grid.sum() - half_space.absorbed.sum()) < 1e-10  # type: ignore[union-attr]
    assert abs(half_space.r_specular + p.weight.sum() + half_space.absorbed.sum() - 1) < 1e-4
    # R(ρ) integrates to the total diffuse reflectance
    edges = np.linspace(0, p.rho.max() + 1, 2001)
    r, _ = det.reflectance_vs_rho(p, edges)
    assert abs(np.sum(r * np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)) - p.weight.sum()) < 1e-12
    # zero modulation frequency = CW; time integral of the TPSF = CW in the same annulus
    cw = det.reflectance_vs_rho(p, np.array([5.0, 6.0]))[0][0]
    assert abs(det.frequency_domain(p, 5.0, 6.0, 0.0).real - cw) < 1e-15
    t_edges = np.linspace(0, 1e5, 20001).astype(np.float64)
    assert abs(np.sum(det.tpsf(p, 5.0, 6.0, t_edges) * np.diff(t_edges)) - cw) < 1e-12
    # partial paths add up; a single layer: ⟨L_1⟩ = ⟨L⟩
    assert np.allclose(p.path.sum(axis=1), p.total_path)


def test_total_reflectance_against_adding_doubling() -> None:
    stack = Stack((Layer(2.0, 0.2, 8.0, 0.8, 1.4),), 1.0, 1.0)
    ad = slab(2.0, 0.2, 8.0, 0.8, 1.4)
    for scheme in ("implicit", "mbl"):
        sim = simulate(stack, 40_000, seed=8, absorption=scheme)
        top = sim.protocol.weight[sim.protocol.side == -1]
        err = np.sqrt(np.var(top * sim.protocol.packets) * top.size) / sim.protocol.packets
        assert abs(top.sum() - ad["r_diffuse"]) < 4 * err + 2e-3


def test_microscopic_beer_lambert_reweighting_is_exact() -> None:
    stack_a = Stack((Layer(1.5, 0.05, 8.0, 0.85, 1.37), Layer(1.0, 0.2, 5.0, 0.9, 1.4)))
    stack_b = Stack((Layer(1.5, 0.12, 8.0, 0.85, 1.37), Layer(1.0, 0.03, 5.0, 0.9, 1.4)))
    a = simulate(stack_a, 5_000, seed=4, absorption="mbl", w_threshold=0.0)
    b = simulate(stack_b, 5_000, seed=4, absorption="mbl", w_threshold=0.0)
    # same seed, no roulette: identical paths; re-weighting a to the absorption of b gives b
    assert np.array_equal(a.protocol.path, b.protocol.path)
    w = reweight(a, np.array([0.12, 0.03]))
    assert np.max(np.abs(w - b.protocol.weight)) < 1e-15
    # ⟨L_i⟩ of the detected light = −∂ ln R/∂μa,i (microscopic Beer–Lambert)
    top = a.protocol.side == -1
    eps = 1e-6
    for i in range(2):
        dmu = np.zeros(2)
        dmu[i] = eps
        r0 = reweight(a, a.protocol.mua_used)[top].sum()
        r1 = reweight(a, a.protocol.mua_used + dmu)[top].sum()
        lw = (a.protocol.weight[top] * a.protocol.path[top, i]).sum() / a.protocol.weight[top].sum()
        assert abs(-(np.log(r1) - np.log(r0)) / eps - lw) < 1e-4 * lw


@pytest.mark.parametrize("n", [1.0, 1.4])
def test_reflectance_rho_against_diffusion(n: float) -> None:
    sim = simulate(Stack((Layer(np.inf, 0.05, 10.0, 0.9, n),)), 60_000, seed=2)
    edges = np.array([4.0, 5.0, 6.0, 7.0, 8.0])
    r, e = det.reflectance_vs_rho(sim.protocol, edges)
    kp = dif.reflectance_rho_kp(0.5 * (edges[1:] + edges[:-1]), 0.05, 1.0, n)
    # Kienle–Patterson in 4–8 transport lengths: within 10 % (measured 5–9 %)
    assert np.all(np.abs(r / kp - 1) < 0.10 + 4 * e / kp)


def test_fluence_and_tpsf_against_diffusion(half_space: Simulation) -> None:
    # laterally integrated fluence (all radial bins, the last collects everything outside)
    z = ((np.arange(60) + 0.5) * 0.25).astype(np.float64)
    phi_1d = half_space.absorbed_grid.sum(axis=0) / (0.05 * 0.25)  # type: ignore[union-attr]
    ref = dif.fluence_laterally_integrated(z, 0.05, 1.0, 1.0)
    depth = (z > 3.0) & (z < 10.0)
    assert np.all(np.abs(phi_1d[depth] / ref[depth] - 1) < 0.08)
    # mean time of flight in the annulus 4–7 mm against the diffusion TPSF (measured +11 %:
    # diffusion is early at a few transport lengths)
    p = half_space.protocol
    sel = (p.side == -1) & (p.rho >= 4) & (p.rho < 7)
    t_mc = float(np.sum(p.weight[sel] * p.time_ps[sel]) / np.sum(p.weight[sel]))
    tt = np.linspace(1.0, 5000.0, 50_000).astype(np.float64)
    num = den = 0.0
    for r in np.linspace(4.0, 7.0, 31):
        f = dif.tpsf(float(r), tt, 0.05, 1.0, 1.0)
        num += float(r * np.trapezoid(f * tt, tt))
        den += float(r * np.trapezoid(f, tt))
    assert abs(t_mc / (num / den) - 1) < 0.15
