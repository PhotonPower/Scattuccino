"""S8b: laser Doppler – tagged moving scatterers, Doppler spectrum (T9b, bonner1981).

The reference (bonner1981) is checked against its own definitions by quadrature (Eqs. 2, 18, 20,
30, 31, 36/A4); the Monte Carlo with tagged moving scatterers against Eq. 22 with the measured
distribution of moving collisions and I₁ from Eq. 2 by quadrature (independent of the transport).
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import quad

from scattuccino.dynamics import (
    Motion,
    doppler_spectrum,
    g1_tagged,
    mean_frequency,
    unshifted_fraction,
)
from scattuccino.layers import Layer, Stack
from scattuccino.phase import Tabulated
from scattuccino.reference import (
    bonner_g2,
    bonner_i1,
    bonner_i1_large,
    bonner_mean_frequency,
    bonner_s_j,
)
from scattuccino.transport import Simulation, simulate

LAM = 0.785
XI = 0.1  # bonner1981 Eq. 13'


def test_i1_against_quadrature_of_eq2() -> None:
    # Eq. 2 with the structure factor of Eq. 13' and Gaussian speeds (Eq. 12'):
    # I₁ = ∫ S(Q) exp(−Q²⟨V²⟩τ²/6) sin θ dθ / ∫ S(Q) sin θ dθ, Q = 2k sin(θ/2)
    k, a, v = 2 * np.pi * 1.33 / (LAM * 1e-3), 1e-3, 1.0  # 1/mm, mm, mm/s
    for tau in (1e-4, 5e-4, 2e-3):

        def num(th: float, tau: float = tau, w: bool = True) -> float:
            q = 2 * k * np.sin(th / 2)
            s = np.exp(-2 * XI * (q * a) ** 2)
            return float(s * (np.exp(-(q**2) * v**2 * tau**2 / 6) if w else 1.0) * np.sin(th))

        ref = quad(num, 0, np.pi, limit=400)[0] / quad(lambda t: num(t, w=False), 0, np.pi)[0]
        assert abs(bonner_i1(tau, a, v, k, XI) - ref) < 1e-8
    # large particles (2ξL² ≫ 1): Eq. 20
    t = np.array([1e-4, 1e-3])
    assert np.allclose(bonner_i1(t, 3e-3, v, k, XI), bonner_i1_large(t, 3e-3, v, XI), rtol=1e-6)


def test_spectra_of_powers_of_i1() -> None:
    # Eq. 30 with I₁(Z) = 1/(1 + Z²): S_j(W) = (1/π) ∫₀^∞ cos(WZ) I₁^j dZ by quadrature; Eq. 31
    # prints S₁ = e^(−W), Eq. 30 gives ½ e^(−W) (consistent with S₂, S₃ and Eq. A3)
    for j in (1, 2, 3):
        for w in (0.3, 1.0, 3.0):
            ref = quad(lambda z, j=j: (1 + z * z) ** -j, 0, np.inf, weight="cos", wvar=w)[0]
            assert abs(bonner_s_j(j, w) - ref / np.pi) < 1e-8


def test_mean_frequency_of_the_poisson_model() -> None:
    # Eq. 36/A4 against the first moment (Eq. 32) of the spectrum of Eq. 24 (Eq. 26), numerically
    a, v, beta = 2.75e-3, 1.0, 0.5
    for mbar in (0.2, 1.0, 3.0):
        tau = np.linspace(0, 0.5, 50_001).astype(np.float64)
        g2m1 = bonner_g2(tau, mbar, a, v, XI, beta) - 1
        omega = np.linspace(0, 1.5e4, 3001).astype(np.float64)
        s = doppler_spectrum(tau, g2m1, omega)
        ref = bonner_mean_frequency(mbar, a, v, XI, beta)
        assert abs(2 * float(np.trapezoid(omega * s, omega)) / ref - 1) < 2e-3, mbar
        assert abs(mean_frequency(tau, g2m1) / ref - 1) < 2e-3, mbar


@pytest.fixture(scope="module")
def tissue() -> tuple[Simulation, Tabulated, float]:
    # static tissue (HG 0.8) with a small fraction of moving cells (Eq. 13' structure factor)
    k = 2 * np.pi * 1.4 / (LAM * 1e-3)
    a = 2.75e-3
    fine = 1 - np.geomspace(1e-9, 2, 40_001)[::-1]  # resolves the forward peak
    mu = np.unique(np.concatenate([[-1.0], fine, [1.0]])).astype(np.float64)
    cell = Tabulated(mu, np.exp(-2 * XI * 2 * (k * a) ** 2 * (1 - mu)).astype(np.float64))
    layer = Layer(np.inf, 0.01, 10.0, 0.8, 1.4, moving_fraction=0.02, phase_moving=cell)
    return simulate(Stack((layer,)), 40_000, seed=5), cell, k


def test_tagging_bookkeeping(tissue: tuple[Simulation, Tabulated, float]) -> None:
    sim, cell, _ = tissue
    p = sim.protocol
    assert p.n_moving is not None and p.momentum_moving is not None and p.momentum is not None
    frac = p.n_moving.sum() / p.scatterings.sum()
    n = p.scatterings.sum()
    assert abs(frac - 0.02) < 4 * np.sqrt(0.02 * 0.98 / n)
    assert np.all(p.momentum_moving <= p.momentum + 1e-12)
    # ⟨1 − cos θ⟩ of the moving events = that of the cell phase function
    m1 = 1 - cell.g
    mean = p.momentum_moving.sum() / p.n_moving.sum()
    assert abs(mean / m1 - 1) < 0.05


def test_monte_carlo_against_eq22(tissue: tuple[Simulation, Tabulated, float]) -> None:
    # g1 = Σ_m P_m I₁^m (Eq. 22 with P₀ kept: unshifted light) with P_m from the protocol and I₁
    # from Eq. 2 for the sampled phase function by quadrature (Gaussian speeds, v = 1 mm/s)
    sim, cell, k = tissue
    p = sim.protocol
    assert p.n_moving is not None
    sel = (p.side == -1) & (p.rho >= 1.0) & (p.rho < 3.0)
    w = p.weight[sel]
    m = p.n_moving[sel, 0].astype(np.int64)
    pm = np.bincount(m, weights=w) / w.sum()
    tau = np.array([2e-4, 1e-3, 5e-3])
    motion = (Motion(v=1.0),)
    g = g1_tagged(p, motion, LAM, tau, 1.0, 3.0)
    for t, gi in zip(tau, g, strict=True):
        # I₁ = ∫ p(μ) exp(−k²(1 − μ)v²τ²/3) dΩ (Eq. 2 with Eq. 12', q² = 2k²(1 − μ))
        i1 = quad(
            lambda mu, t=t: 2
            * np.pi
            * cell.pdf(np.array([mu]))[0]
            * np.exp(-(k**2) * (1 - mu) * t * t / 3),
            -1,
            1,
            points=[0.99, 0.999, 0.9999],
            limit=500,
        )[0]
        ref = float(np.sum(pm * i1 ** np.arange(pm.size)))
        assert abs(gi - ref) < 0.01, (t, gi, ref)
    assert abs(unshifted_fraction(p, 1.0, 3.0) - pm[0]) < 1e-12
    assert 0.05 < pm[0] < 0.95  # both shifted and unshifted light present
