"""S1b: adding-doubling (independent of the Monte Carlo) and Monte Carlo against it (T1c)."""

from __future__ import annotations

import numpy as np
import pytest

from scattuccino.adding_doubling import slab
from scattuccino.layers import Layer, Stack, run
from scattuccino.reference import semi_infinite_isotropic_reflectance, slab_without_scattering


@pytest.mark.parametrize("albedo", [0.5, 0.9, 0.99])
def test_adding_doubling_against_h_function(albedo: float) -> None:
    r = slab(np.inf, 1 - albedo, albedo, 0.0, 1.0)["r_diffuse"]
    assert abs(r - semi_infinite_isotropic_reflectance(albedo)) < 1e-4  # measured ≤ 5e-5


def test_adding_doubling_without_scattering_is_exact() -> None:
    res = slab(1.0, 0.5, 0.0, 0.0, 1.5)
    r_ref, t_ref = slab_without_scattering(0.5, 1.0, 1.5)
    assert abs(res["r_specular"] + res["r_diffuse"] - r_ref) < 1e-12
    assert abs(res["t_diffuse"] - t_ref) < 1e-11


def test_adding_doubling_literature_case() -> None:
    # μa = 1, μs = 9 /mm, g = 0.75, d = 0.2 mm, n = 1: van de Hulst's values as quoted for the
    # validation of MCML (wang1995; status Z in docs/quellen.md)
    res = slab(0.2, 1.0, 9.0, 0.75, 1.0)
    assert abs(res["r_diffuse"] - 0.09734) < 2e-4
    assert abs(res["t_diffuse"] - 0.66096) < 2e-4


CASES = [
    (1.0, 0.1, 10.0, 0.9, 1.4, 1.0, 1.0, 1.0),
    (np.inf, 0.1, 0.9, 0.75, 1.33, 1.0, 1.0, 1.0),
    (0.5, 0.3, 6.0, 0.7, 1.45, 1.0, 1.52, 0.7),
    (0.2, 0.05, 20.0, 0.95, 1.37, 1.33, 1.0, 1.0),
]


@pytest.mark.parametrize("case", CASES)
def test_monte_carlo_against_adding_doubling(case: tuple[float, ...]) -> None:
    d, mua, mus, g, n, n_above, n_below, cos_i = case
    ad = slab(d, mua, mus, g, n, n_above, n_below, cos_i)
    mc = run(
        Stack((Layer(d, mua, mus, g, n),), n_above, n_below), 100_000, seed=7, cos_incident=cos_i
    )
    assert abs(mc.r_specular - ad["r_specular"]) < 1e-15
    assert abs(mc.r_diffuse - ad["r_diffuse"]) < 4 * mc.r_diffuse_err + 1e-4
    assert abs(mc.t_diffuse - ad["t_diffuse"]) < 4 * mc.t_diffuse_err + 1e-4
