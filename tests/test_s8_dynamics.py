"""S8a: dynamic light scattering – momentum transfer, DCS, speckle contrast (T9)."""

from __future__ import annotations

import numpy as np
import pytest

from scattuccino.diffusion import g1_correlation_diffusion
from scattuccino.dynamics import (
    Motion,
    g1,
    g2,
    speckle_contrast,
    speckle_contrast_single_exponential,
)
from scattuccino.layers import Layer, Stack
from scattuccino.transport import Simulation, simulate

LAM = 0.785


@pytest.fixture(scope="module")
def half_space() -> Simulation:
    return simulate(Stack((Layer(np.inf, 0.05, 10.0, 0.9, 1.4),)), 60_000, seed=3)


def test_momentum_transfer_bookkeeping(half_space: Simulation) -> None:
    p = half_space.protocol
    assert p.momentum is not None
    # ⟨1 − cos θ⟩ = 1 − g for Henyey–Greenstein
    n = p.scatterings.sum()
    assert abs(p.momentum.sum() / n - 0.1) < 4 * np.sqrt(0.81 / 3 / n)  # var(cos) ≤ (1−g²)/3
    assert np.all(p.momentum >= 0) and np.all(p.momentum <= 2 * p.scatterings[:, None] + 1e-12)
    # static scatterers: g1 = 1
    tau = np.logspace(-6, -2, 5).astype(np.float64)
    assert np.allclose(g1(p, (Motion(),), LAM, tau, 2.0, 3.0), 1.0)


def test_dcs_against_correlation_diffusion(half_space: Simulation) -> None:
    tau = np.logspace(-6, -2, 13).astype(np.float64)
    motion = Motion(alpha=1.0, d_b=1e-6)  # αD_B = 1e-8 cm²/s
    mc = g1(half_space.protocol, (motion,), LAM, tau, 2.0, 3.0)
    ref = g1_correlation_diffusion(2.5, tau, 0.05, 1.0, 1.4, LAM, 1.0, motion.msd(tau))
    assert np.all(np.diff(mc) <= 1e-12) and mc[0] > 0.99
    valid = ref > 0.1
    assert np.max(np.abs(mc[valid] - ref[valid])) < 0.04  # diffusion at 2–3 transport lengths
    assert np.allclose(g2(mc, 0.5), 1 + 0.5 * mc**2)


def test_speckle_contrast_single_exponential() -> None:
    tau_c, exposure = 1e-3, 5e-3
    t = np.linspace(0.0, exposure, 20_001).astype(np.float64)
    k = speckle_contrast(t, np.exp(-t / tau_c), exposure)
    assert abs(k - speckle_contrast_single_exponential(tau_c, exposure)) < 1e-7
    assert (
        abs(speckle_contrast_single_exponential(tau_c, 1e-9) - 1) < 1e-5
    )  # short exposure: K → √β


def test_depth_sensitivity_two_layers() -> None:
    # flow only in the deep layer (below 1 mm): decorrelation grows with the distance
    stack = Stack((Layer(1.0, 0.05, 10.0, 0.9, 1.4), Layer(np.inf, 0.05, 10.0, 0.9, 1.4)))
    sim = simulate(stack, 40_000, seed=4)
    mot = (Motion(), Motion(alpha=1.0, d_b=1e-6))
    tau = np.array([3e-4])
    near = g1(sim.protocol, mot, LAM, tau, 0.5, 1.0)[0]
    far = g1(sim.protocol, mot, LAM, tau, 3.0, 4.0)[0]
    assert far < near - 0.1
