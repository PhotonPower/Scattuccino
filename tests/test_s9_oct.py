"""S9a: OCT A-scans from the photon protocol (T10)."""

from __future__ import annotations

import numpy as np
import pytest

from scattuccino import oct as o
from scattuccino.layers import Layer, Stack
from scattuccino.phase import HenyeyGreenstein
from scattuccino.transport import Simulation, simulate

MUA, MUS, G = 0.1, 2.0, 0.3


@pytest.fixture(scope="module")
def half_space() -> Simulation:
    return simulate(Stack((Layer(np.inf, MUA, MUS, G, 1.0),)), 400_000, seed=1)


def test_single_backscattering_against_exact_integral(half_space: Simulation) -> None:
    p = half_space.protocol
    theta = np.radians(30)
    edges = np.linspace(0.0, 1.0, 11).astype(np.float64)
    sel = o.detected(p, 1e9, theta) & (p.scatterings == 1)
    zo, w = o.depth_opl(p)[sel], p.weight[sel]
    h, _ = np.histogram(zo, bins=edges, weights=w)
    e2, _ = np.histogram(zo, bins=edges, weights=(w * p.packets) ** 2)
    err = np.sqrt(e2) / p.packets
    ref = o.single_scattering_a_scan(edges, MUA, MUS, HenyeyGreenstein(G).pdf, theta)
    assert np.all(np.abs(h - ref) < 4 * err)


def test_coherence_gating_conserves_the_signal(half_space: Simulation) -> None:
    p = half_space.protocol
    z = np.linspace(-0.5, 4.0, 901).astype(np.float64)
    hist = o.a_scan(p, z, 0.1, np.radians(10))
    gated = o.a_scan(p, z, 0.1, np.radians(10), l_c=0.01)
    dz = z[1] - z[0]
    sel = o.detected(p, 0.1, np.radians(10)) & (o.depth_opl(p) < 3.9)
    total = p.weight[sel].sum()
    assert abs(np.sum(gated[z < 3.95]) * dz - total) < 1e-3 * total
    assert abs(np.sum(hist[z < 3.95]) * dz - total) < 2e-2 * total  # binning at the edge


def test_multiple_scattering_grows_with_depth(half_space: Simulation) -> None:
    p = half_space.protocol
    z = np.linspace(0.05, 1.95, 20).astype(np.float64)
    all_ = o.a_scan(p, z, 1.0, np.radians(30))
    single = o.a_scan(p, z, 1.0, np.radians(30), mask=p.scatterings == 1)
    share = 1 - single / np.maximum(all_, 1e-300)
    assert share[0] < 0.4 and share[-1] > 0.9  # shallow: mostly single, deep: multiple
    # the single-scattering part decays like e^(−2 μt z) (Beer law of OCT)
    slope = np.polyfit(z[:8], np.log(np.maximum(single[:8], 1e-300)), 1)[0]
    assert abs(slope / (-2 * (MUA + MUS)) - 1) < 0.05  # measured 1.5 %
