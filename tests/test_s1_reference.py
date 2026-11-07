"""S1: Monte-Carlo reference in layers against independent solutions (docs/testfaelle.md).

Statistical tests use 4 standard errors; seeds are fixed.
"""

from __future__ import annotations

import numpy as np
import pytest

from scattuccino.layers import Layer, Stack, fresnel_unpolarised, run, sample_hg
from scattuccino.reference import semi_infinite_isotropic_reflectance, slab_without_scattering


def test_slab_without_scattering() -> None:
    res = run(Stack((Layer(1.0, 0.5, 0.0, 0.0, 1.5),)), 200_000, seed=3)
    r_ref, t_ref = slab_without_scattering(0.5, 1.0, 1.5)
    assert abs(res.r_specular + res.r_diffuse - r_ref) < 4 * res.r_diffuse_err
    assert abs(res.t_diffuse - t_ref) < 4 * res.t_diffuse_err


@pytest.mark.parametrize("albedo", [0.5, 0.9, 0.99])
def test_semi_infinite_isotropic_against_h_function(albedo: float) -> None:
    res = run(Stack((Layer(np.inf, 1 - albedo, albedo, 0.0, 1.0),)), 100_000, seed=4)
    ref = semi_infinite_isotropic_reflectance(albedo)
    assert abs(res.r_diffuse - ref) < 4 * res.r_diffuse_err


def test_energy_conservation_two_layers_with_index_steps() -> None:
    stack = Stack((Layer(0.3, 0.2, 8.0, 0.8, 1.4), Layer(2.0, 0.05, 5.0, 0.9, 1.37)), 1.0, 1.52)
    res = run(stack, 50_000, seed=5, cos_incident=np.cos(0.4))
    assert abs(res.total - 1.0) < 1e-4  # roulette conserves energy in expectation only
    assert res.absorbed.shape == (2,) and np.all(res.absorbed > 0)


def test_henyey_greenstein_sampling_mean() -> None:
    rng = np.random.default_rng(1)
    for g in (0.0, 0.5, 0.9, -0.3):
        c = sample_hg(np.full(400_000, g), rng.random(400_000))
        assert abs(c.mean() - g) < 4 * c.std() / np.sqrt(c.size)


def test_fresnel_normal_and_total_internal_reflection() -> None:
    r, _ = fresnel_unpolarised(1.0, 1.5, np.array([1.0]))
    assert abs(r[0] - 0.04) < 1e-15
    r_tir, _ = fresnel_unpolarised(1.5, 1.0, np.array([np.cos(np.radians(60))]))
    assert r_tir[0] == 1.0


def test_determinism() -> None:
    a = run(Stack((Layer(1.0, 0.1, 2.0, 0.7, 1.33),)), 5_000, seed=9)
    b = run(Stack((Layer(1.0, 0.1, 2.0, 0.7, 1.33),)), 5_000, seed=9)
    assert a.r_diffuse == b.r_diffuse and a.t_diffuse == b.t_diffuse
