"""S7c: fluorescence emission spectra (T8c).

The emission is split into wavelength bins with the normalised emission spectrum f_b; each bin is
transported with its own optical properties. Checks:

1. Non-scattering medium: per bin F_b = Q f_b ∫ μa,f e^(−μa,x z) ½E₂(μa,m,b z) dz (closed form of
   T8 per bin) – the detected spectrum shows the inner-filter distortion by re-absorption.
2. One excitation for all bins (equal to `simulate_fluorescence`), energy per bin.
3. Equal properties in all bins: the detected total equals the single-wavelength run.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad
from scipy.special import expn

from scattuccino.fluorescence import simulate_fluorescence, simulate_fluorescence_spectrum
from scattuccino.layers import Layer, Stack


def test_non_scattering_spectrum_against_closed_form() -> None:
    mua_t, mua_f, q = 0.5, 0.3, 0.2
    mux = mua_t + mua_f
    sx = Stack((Layer(40.0, mux, 0.0, 0.0, 1.0),))
    mua_m = np.array([0.8, 0.4, 0.1])  # strong re-absorption in the blue bin
    f_b = np.array([0.2, 0.5, 0.3])
    sm = [Stack((Layer(40.0, float(m), 0.0, 0.0, 1.0),)) for m in mua_m]
    n_m = 120_000
    res = simulate_fluorescence_spectrum(sx, sm, f_b, np.array([mua_f]), q, 100_000, n_m, seed=3)
    assert abs(res.absorbed_fluorophore - mua_f / mux) < 1e-12
    assert np.allclose(res.emitted, q * f_b * res.absorbed_fluorophore, rtol=1e-12)
    detected = res.detected()
    for b, m in enumerate(mua_m):
        ref = (
            q
            * f_b[b]
            * quad(
                lambda z, m=m: mua_f * np.exp(-mux * z) * 0.5 * expn(2, m * z),
                0,
                np.inf,
                limit=200,
            )[0]
        )
        n_b = res.packets[b]
        frac = ref / res.emitted[b]
        sigma = res.emitted[b] * np.sqrt(frac * (1 - frac) / n_b)
        assert abs(detected[b] - ref) < 4 * sigma, b
    # inner filter: the detected spectrum is shifted to the weakly re-absorbed bins
    shape = detected / detected.sum()
    assert shape[0] < f_b[0] and shape[2] > f_b[2]


def test_one_excitation_and_energy_per_bin() -> None:
    sx = Stack((Layer(2.0, 0.15, 10.0, 0.9, 1.4),))
    sm = [Stack((Layer(2.0, m, 9.0, 0.9, 1.4),)) for m in (0.1, 0.06)]
    f_b = np.array([0.4, 0.6])
    res = simulate_fluorescence_spectrum(sx, sm, f_b, np.array([0.05]), 0.1, 30_000, 40_000, seed=4)
    single = simulate_fluorescence(sx, sm[0], np.array([0.05]), 0.1, 30_000, 1000, seed=4)
    assert res.absorbed_fluorophore == single.absorbed_fluorophore
    for b, emi in enumerate(res.emission):
        p = emi.protocol
        total = p.weight.sum() + emi.absorbed.sum()
        # Russian roulette conserves energy in expectation only (measured 1e-5)
        assert abs(total - res.emitted[b]) < 1e-4 * res.emitted[b]


def test_equal_bins_equal_the_single_wavelength_run() -> None:
    sx = Stack((Layer(np.inf, 0.15, 10.0, 0.9, 1.4),))
    sm = Stack((Layer(np.inf, 0.06, 9.0, 0.9, 1.4),))
    n = 60_000
    res = simulate_fluorescence_spectrum(
        sx, [sm, sm, sm], np.array([0.2, 0.3, 0.5]), np.array([0.05]), 0.1, 40_000, n, seed=6
    )
    single = simulate_fluorescence(sx, sm, np.array([0.05]), 0.1, 40_000, n, seed=6)
    p = single.emission.protocol
    top = p.weight[p.side == -1]
    ref = top.sum()
    frac = ref / single.emitted
    sigma = single.emitted * np.sqrt(frac * (1 - frac) / n)  # binomial, both runs independent
    assert abs(res.detected().sum() - ref) < 4 * np.sqrt(2) * sigma
