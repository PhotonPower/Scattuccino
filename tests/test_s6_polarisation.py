"""S6a: polarised transport (T7).

1. The scattering step against an independent field calculation: Stokes vectors from the
   amplitudes S1, S2 (Bohren & Huffman) for Mie, and from the dipole field −û′ × (û′ × E) for
   Rayleigh – checks frame rotation, Mueller matrix and the sign of V.
2. Unpolarised light: intensities equal the scalar transport with the same phase function.
3. Pure states stay pure (sphere Mueller matrices map pure to pure).
4. Polarisation memory (MacKintosh et al. 1989): in backscattering from many small spheres linear
   polarisation survives better than circular (and the helicity flips); for large spheres
   circular survives better (helicity kept).
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from scattuccino._types import FloatArray
from scattuccino.layers import Layer, Stack
from scattuccino.mie import s12
from scattuccino.polarisation import (
    apply_scattering,
    mie_mueller,
    rayleigh_mueller,
    simulate_polarised,
    stokes_lab,
)
from scattuccino.transport import simulate

PS = (1.59 / 1.33 + 0j, 2 * np.pi * 0.5 * 1.33 / 0.633)


def _stokes(ev: Any, ew: Any) -> FloatArray:
    i = np.abs(ev) ** 2 + np.abs(ew) ** 2
    q = np.abs(ev) ** 2 - np.abs(ew) ** 2
    u = 2 * np.real(ev * np.conj(ew))
    v = -2 * np.imag(ev * np.conj(ew))
    return np.asarray(np.column_stack([i, q, u, v]) / i[:, None], dtype=np.float64)


def test_scattering_step_against_fields() -> None:
    rng = np.random.default_rng(3)
    n = 200
    u = rng.normal(size=(n, 3))
    u = (u / np.linalg.norm(u, axis=1)[:, None]).astype(np.float64)
    v = np.cross(u, rng.normal(size=(n, 3))).astype(np.float64)
    v /= np.linalg.norm(v, axis=1)[:, None]
    w = np.cross(u, v)
    ev = rng.normal(size=n) + 1j * rng.normal(size=n)
    ew = rng.normal(size=n) + 1j * rng.normal(size=n)
    s_in = _stokes(ev, ew)
    mu = rng.uniform(-1, 1, n)
    phi = rng.uniform(0, 2 * np.pi, n)
    ev1 = np.cos(phi) * ev + np.sin(phi) * ew
    ew1 = -np.sin(phi) * ev + np.cos(phi) * ew
    # Mie through the amplitudes
    m, x = 1.2 + 0.01j, 3.0
    _, _, s_out = apply_scattering(mie_mueller(m, x, 20_001), u, v, s_in, mu, phi)
    a1, a2 = s12(m, x, np.arccos(mu))
    assert np.max(np.abs(s_out - _stokes(a2 * ev1, a1 * ew1))) < 1e-6  # interpolation of the table
    # Rayleigh through the dipole field
    u2, v2, s_out = apply_scattering(rayleigh_mueller(20_001), u, v, s_in, mu, phi)
    e_in = ev[:, None] * v + ew[:, None] * w
    e_s = e_in - np.sum(u2 * e_in, axis=1)[:, None] * u2
    w2 = np.cross(u2, v2)
    ref = _stokes(np.sum(e_s * v2, axis=1), np.sum(e_s * w2, axis=1))
    assert np.max(np.abs(s_out - ref)) < 1e-6
    assert np.allclose(np.sum(u2 * v2, axis=1), 0, atol=1e-12)


def test_unpolarised_intensities_equal_scalar_transport() -> None:
    table = mie_mueller(*PS)
    pol = simulate_polarised(1.0, 0.1, 10.0, table, np.array([1.0, 0, 0, 0]), 40_000, seed=2)
    sca = simulate(Stack((Layer(1.0, 0.1, 10.0, 0.0, 1.0, phase=table.phase),)), 40_000, seed=3)
    p = sca.protocol
    for side in (-1, 1):
        a = pol.weight[pol.side == side].sum()
        b = p.weight[p.side == side].sum()
        assert abs(a - b) < 0.006  # two independent runs (4σ of the difference)


def test_pure_states_stay_pure() -> None:
    res = simulate_polarised(
        1.0, 0.1, 10.0, mie_mueller(*PS), np.array([1.0, 0.6, 0.0, 0.8]), 5_000, seed=4
    )
    assert np.max(np.abs(np.sum(res.stokes[:, 1:] ** 2, axis=1) - 1)) < 1e-12


@pytest.mark.parametrize(("x", "linear_wins"), [(0.5, True), (10.0, False)])
def test_polarisation_memory(x: float, linear_wins: bool) -> None:
    table = mie_mueller(1.59 / 1.33 + 0j, x)
    dop = {}
    for name, s in (("lin", np.array([1.0, 1.0, 0, 0])), ("circ", np.array([1.0, 0, 0, 1.0]))):
        res = simulate_polarised(np.inf, 0.0, 10.0, table, s, 15_000, seed=5, max_events=3000)
        top = res.side == -1
        st = stokes_lab(res)[top]
        w = res.weight[top]
        dop[name] = np.sum(w * (st[:, 1] if name == "lin" else st[:, 3])) / w.sum()
    if linear_wins:
        assert dop["lin"] > abs(dop["circ"]) + 0.05 and dop["circ"] < 0  # helicity flips
    else:
        assert dop["circ"] > abs(dop["lin"]) + 0.1 and dop["circ"] > 0  # helicity kept
