"""S9b-2: random index fields and their scattering (T16).

Born integrals against closed forms; the generated field has the prescribed variance and
correlation on the grid; with BWPM (optional dependency, skipped otherwise) one small realisation
reproduces the Born μs and the forward g.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import quad

from scattuccino.microstructure import (
    WhittleMatern,
    born,
    calibrate_with_bwpm,
    quantise,
    random_field,
)

LAM, N0 = 0.8, 1.38


def test_born_against_closed_form_for_exponential_correlation() -> None:
    med = WhittleMatern(0.02, 0.5e-3, 4.0)
    b = born(med, N0, LAM)
    k = 2 * np.pi * N0 / (LAM * 1e-3)
    a = 2 * k**2 * med.l_c**2
    mus = 8 * k**4 * med.sigma_n**2 * med.l_c**3 / (N0**2 * (1 + 2 * a))
    assert abs(b.mus / mus - 1) < 1e-6
    g_num = quad(lambda m: m / (1 + a * (1 - m)) ** 2, -1, 1, limit=500)[0]
    g_den = quad(lambda m: 1 / (1 + a * (1 - m)) ** 2, -1, 1, limit=500)[0]
    assert abs(b.g - g_num / g_den) < 1e-6


def test_random_field_statistics() -> None:
    med = WhittleMatern(0.05, 0.5e-3, 4.0)
    h, shape = 1e-4, (64, 64, 64)
    dn = random_field(med, shape, h, np.random.default_rng(1))
    variances = [
        random_field(med, shape, h, np.random.default_rng(seed)).var() for seed in range(8)
    ]
    qs = [2 * np.pi * np.fft.fftfreq(n, d=h) for n in shape]
    qz, qy, qx = np.meshgrid(*qs, indexing="ij")
    psd = med.psd(np.sqrt(qx**2 + qy**2 + qz**2))
    psd[0, 0, 0] = 0.0  # the sample variance removes the mean (q = 0)
    expected = np.sum(psd) * (2 * np.pi) ** 3 / (h**3 * dn.size)
    # low-q modes make one realisation fluctuate by ~6 %; eight realisations average it down
    assert abs(np.mean(variances) / expected - 1) < 0.08
    corr = [np.mean(dn * np.roll(dn, s, axis=2)) / np.mean(dn * dn) for s in (2, 5)]
    assert (
        abs(corr[0] - np.exp(-2 * h / med.l_c)) < 0.06
        and abs(corr[1] - np.exp(-5 * h / med.l_c)) < 0.06
    )
    lab, idx = quantise(dn, N0, 64)
    n = np.array(idx)[lab].real - N0
    assert np.max(np.abs(n - dn)) <= (dn.max() - dn.min()) / 126 + 1e-15


def test_bwpm_realisation_reproduces_born() -> None:
    pytest.importorskip("bwpm")
    med = WhittleMatern(0.0572, 0.5e-3, 4.0)  # Born μs ≈ 200/mm, g ≈ 0.935
    cal = calibrate_with_bwpm(
        med, N0, LAM, (60, 96, 96), 1e-4, np.random.default_rng(1), keep_every=20
    )
    assert abs(cal.mus_fit / cal.born.mus - 1) < 0.2  # one small realisation
    assert abs(cal.g_forward - cal.g_forward_born) < 0.01
    assert np.all(np.abs(cal.total - 1) < 0.02)  # forward FPWPM conserves power
