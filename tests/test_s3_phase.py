"""S3: phase functions, Mie and the subdiffusive regime (T4, T5)."""

from __future__ import annotations

import numpy as np
import pytest

from scattuccino import detectors as det
from scattuccino.adding_doubling import slab
from scattuccino.layers import Layer, Stack
from scattuccino.mie import qeff, s12
from scattuccino.phase import Gegenbauer, HenyeyGreenstein, ModifiedHG, Tabulated, mie_phase
from scattuccino.transport import simulate

PS_SPHERE = (1.59 / 1.33 + 0j, 2 * np.pi * 0.5 * 1.33 / 0.633)  # 1 µm polystyrene in water, 633 nm
PHASES = {
    "hg": HenyeyGreenstein(0.9),
    "mhg": ModifiedHG(0.9, 0.8),
    "gk": Gegenbauer(0.8, 1.0),
    "mie": mie_phase(*PS_SPHERE)[0],
}


@pytest.mark.parametrize("name", list(PHASES))
def test_normalisation_and_sampled_moments(name: str) -> None:
    ph = PHASES[name]
    x, w = np.polynomial.legendre.leggauss(2000)
    assert abs(2 * np.pi * np.sum(w * ph.pdf(x.astype(np.float64))) - 1) < 1e-6
    mu = ph.sample(np.random.default_rng(1).random(400_000).astype(np.float64))
    for n, sample in ((1, mu), (2, 1.5 * mu**2 - 0.5)):
        assert abs(sample.mean() - ph.moment(n)) < 5 * sample.std() / np.sqrt(mu.size)


def test_closed_forms() -> None:
    assert abs(ModifiedHG(0.9, 0.8).g - 0.72) < 1e-9  # g = β g_HG
    mu = np.linspace(-1.0, 1.0, 7).astype(np.float64)
    assert np.allclose(Gegenbauer(0.9, 0.5).pdf(mu), HenyeyGreenstein(0.9).pdf(mu))  # α = ½ is HG
    xi = np.linspace(0.01, 0.99, 9).astype(np.float64)
    assert np.allclose(Gegenbauer(0.9, 0.5).sample(xi), HenyeyGreenstein(0.9).sample(xi))
    # a tabulated HG (grid dense in the forward direction) has the moment of the analytic one
    grid = np.unique(np.cos(np.pi * np.linspace(0.0, 1.0, 20_001) ** 2)).astype(np.float64)
    tab = Tabulated(grid, HenyeyGreenstein(0.8).pdf(grid))
    assert abs(tab.g - 0.8) < 1e-4


def test_mie_limits() -> None:
    m, x = 1.5 + 0j, 0.05
    q_ext, q_sca = qeff(m, x)
    rayleigh = 8 / 3 * x**4 * abs((m * m - 1) / (m * m + 2)) ** 2
    assert abs(q_sca / rayleigh - 1) < 1e-3  # measured 1.8e-4 (x⁶ terms)
    s1, s2 = s12(m, x, np.array([0.0, np.pi / 2], dtype=np.float64))
    p = np.abs(s1) ** 2 + np.abs(s2) ** 2  # Rayleigh: ∝ 1 + cos²θ
    assert abs(p[1] / p[0] - 0.5) < 1e-2
    assert abs(qeff(1.33 + 0.01j, 300.0)[0] - 2.0) < 0.06  # extinction paradox, 2.044
    tab, q_s, q_e = mie_phase(*PS_SPHERE)
    assert q_s == pytest.approx(q_e, rel=1e-12)  # non-absorbing
    assert 0.90 < tab.g < 0.93  # 1 µm polystyrene in water: 0.917


@pytest.mark.parametrize("name", ["mhg", "gk", "mie"])
def test_monte_carlo_against_adding_doubling(name: str) -> None:
    ph = PHASES[name]
    ad = slab(1.0, 0.1, 10.0, 0.0, 1.4, phase=ph)
    sim = simulate(Stack((Layer(1.0, 0.1, 10.0, 0.0, 1.4, phase=ph),)), 100_000, seed=11)
    p = sim.protocol
    for side, key in ((-1, "r_diffuse"), (1, "t_diffuse")):
        w = np.where(p.side == side, p.weight * p.packets, 0.0)
        mean = w.sum() / p.packets
        err = np.sqrt((np.sum(w**2) / p.packets - mean**2) / p.packets)
        assert abs(mean - ad[key]) < 4 * err


def test_subdiffusive_reflectance_depends_on_gamma() -> None:
    # same g = 0.9 (same μs′), different γ: near the source the lower γ reflects much more,
    # far away the profiles approach each other (Bevilacqua & Depeursinge 1999)
    edges = np.array([0.0, 0.1, 0.4, 0.8, 6.0, 8.0])
    r = {}
    for name, ph in (("hg", HenyeyGreenstein(0.9)), ("mhg", ModifiedHG(0.9 / 0.92, 0.92))):
        assert abs(ph.g - 0.9) < 1e-9
        sim = simulate(Stack((Layer(np.inf, 0.05, 10.0, 0.0, 1.4, phase=ph),)), 40_000, seed=4)
        r[name] = det.reflectance_vs_rho(sim.protocol, edges)[0]
    assert ModifiedHG(0.9 / 0.92, 0.92).gamma < HenyeyGreenstein(0.9).gamma
    ratio = r["mhg"] / r["hg"]
    assert ratio[0] > 3.0  # 0–0.1 mm: about 5×
    assert abs(ratio[-1] - 1) < abs(ratio[0] - 1) / 5
