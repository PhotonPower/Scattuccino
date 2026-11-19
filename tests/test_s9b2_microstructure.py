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
    bwpm_medium,
    calibrate_with_bwpm,
    plane_wave_power,
    quantise,
    random_field,
)

LAM, N0 = 0.8, 1.38


def test_born_against_closed_form_for_exponential_correlation() -> None:
    med = WhittleMatern(0.02, 0.5e-3, 4.0)
    b = born(med, N0, LAM, vector=False)
    k = 2 * np.pi * N0 / (LAM * 1e-3)
    a = 2 * k**2 * med.l_c**2
    mus = 8 * k**4 * med.sigma_n**2 * med.l_c**3 / (N0**2 * (1 + 2 * a))
    assert abs(b.mus / mus - 1) < 1e-6
    g_num = quad(lambda m: m / (1 + a * (1 - m)) ** 2, -1, 1, limit=500)[0]
    g_den = quad(lambda m: 1 / (1 + a * (1 - m)) ** 2, -1, 1, limit=500)[0]
    assert abs(b.g - g_num / g_den) < 1e-6


def test_vector_born_has_the_dipole_factor() -> None:
    # S9b-7: σ_vec = σ_scalar (1 + cos²θ)/2 (unpolarised); exponential correlation, closed form of
    # σ_scalar ∝ 1/(1 + a(1 − μ))², integrals by quadrature
    med = WhittleMatern(0.02, 0.5e-3, 4.0)
    b = born(med, N0, LAM)
    k = 2 * np.pi * N0 / (LAM * 1e-3)
    a = 2 * k**2 * med.l_c**2
    c = 2 * np.pi * k**4 / N0**2 * med.sigma_n**2 * med.l_c**3 / np.pi**2  # σ at a(1 − μ) = 0

    def shape(m: float) -> float:
        return (1 + m * m) / 2 / (1 + a * (1 - m)) ** 2

    mus = 2 * np.pi * c * quad(shape, -1, 1, limit=500)[0]
    g = quad(lambda m: m * shape(m), -1, 1, limit=500)[0] / quad(shape, -1, 1, limit=500)[0]
    assert abs(b.mus / mus - 1) < 1e-6
    assert abs(b.g - g) < 1e-6
    mu = np.array([-1.0, 0.0, 1.0])  # backward and forward unchanged, 90° halved
    scalar = born(med, N0, LAM, vector=False)
    ratio = b.phase.pdf(mu) * b.mus / (scalar.phase.pdf(mu) * scalar.mus)
    assert np.allclose(ratio, [1.0, 0.5, 1.0], rtol=1e-6)


def test_plane_wave_power_includes_e_z() -> None:
    # a plane wave at 40° in the x–z plane, E ⊥ k (p-polarised): |E|² = |E_x|²/cos²θ; s-polarised
    # (E along y): |E|² = |E_y|²
    n, h = 64, 1e-4
    k = 2 * np.pi * N0 / (LAM * 1e-3)
    kx = 2 * np.pi * np.fft.fftfreq(n, h)
    kxx, kyy = np.meshgrid(kx, kx)
    m = round(float(k * np.sin(np.radians(40)) * n * h / (2 * np.pi)))
    cos_t = np.sqrt(1 - (kx[m] / k) ** 2)
    x = np.arange(n) * h
    wave = np.exp(1j * kx[m] * x)[None, :] * np.ones((n, 1))
    for comp, expected in ((0, 1 / cos_t**2), (1, 1.0)):
        f = np.zeros((n, n, 2), dtype=np.complex128)
        f[:, :, comp] = wave
        s = plane_wave_power(f, kxx, kyy, k)
        assert abs(s[0, m] / n**4 - expected) < 1e-9 and abs(s.sum() - s[0, m]) < 1e-9 * s[0, m]


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
    med = WhittleMatern(0.0572, 0.5e-3, 4.0)  # vector Born μs ≈ 190/mm, g ≈ 0.948
    cal = calibrate_with_bwpm(
        med, N0, LAM, (60, 96, 96), 1e-4, np.random.default_rng(1), keep_every=20
    )
    # one small realisation; measured 1.010 against vector Born (0.960 against scalar Born). The
    # FPWPM conserves the transverse norm, so its extinction is ∫ σ r²/cos θ dΩ (BWPM K17 point 8):
    # 1.028 × vector Born – the closeness to vector Born is no proof of vector scattering
    assert abs(cal.mus_fit / cal.born.mus - 1) < 0.1
    assert np.all(np.abs(cal.total - 1) < 0.02)  # forward FPWPM conserves power


@pytest.mark.xfail(
    strict=True,
    reason="BWPM E9 (K17): the FPWPM scatters the p part too strongly at large angles "
    "(r/cos²θ); full field against vector Born gives g_forward 0.913 vs 0.955",
)
def test_bwpm_realisation_reproduces_born_g() -> None:
    pytest.importorskip("bwpm")
    med = WhittleMatern(0.0572, 0.5e-3, 4.0)
    cal = calibrate_with_bwpm(
        med, N0, LAM, (60, 96, 96), 1e-4, np.random.default_rng(1), keep_every=20
    )
    assert abs(cal.g_forward - cal.g_forward_born) < 0.02  # μs z = 0.4: some multiple scattering


def test_measured_phase_function_for_the_monte_carlo() -> None:
    pytest.importorskip("bwpm")
    med = WhittleMatern(0.0572, 0.5e-3, 4.0)
    cal = calibrate_with_bwpm(
        med, N0, LAM, (60, 96, 96), 1e-4, np.random.default_rng(1), keep_every=20
    )
    mus, phase = cal.scattering()
    x, w = np.polynomial.legendre.leggauss(4000)
    assert abs(2 * np.pi * np.sum(w * phase.pdf(x.astype(np.float64))) - 1) < 1e-3
    assert mus == cal.mus_fit
    # the backward half (not seen by the forward propagation) is Born in absolute terms
    mu_b = np.array([-1.0, -0.5])
    sigma = phase.pdf(mu_b) * mus
    sigma_born = cal.born.phase.pdf(mu_b) * cal.born.mus
    assert np.all(np.abs(sigma / sigma_born - 1) < 0.05)


def test_bwpm_medium_models() -> None:
    # S9b-8: δn as index steps (default) or as the fluctuation on one host index
    dn = random_field(WhittleMatern(0.05, 0.5e-3, 4.0), (4, 16, 16), 1e-4, np.random.default_rng(3))
    lab, idx, fl = bwpm_medium(dn, N0)
    lab_q, idx_q = quantise(dn, N0)
    assert fl is None and np.array_equal(lab, lab_q) and idx == idx_q
    lab, idx, fl = bwpm_medium(dn, N0, "increment")
    assert idx == [complex(N0)] and not lab.any() and fl is dn
    with pytest.raises(ValueError):
        bwpm_medium(dn, N0, "other")


@pytest.mark.xfail(
    strict=True,
    reason="BWPM B7: for σ_n ≈ 0.06 the Born increment pumps power into near-grazing waves "
    "(|E⊥|² at cos θ < 0.3: 0.03 / 0.20 / 0.87 after 2 / 4 / 6 µm; BWPM #29)",
)
def test_bwpm_born_increment_conserves_power() -> None:
    pytest.importorskip("bwpm.fluctuation")
    med = WhittleMatern(0.0572, 0.5e-3, 4.0)
    cal = calibrate_with_bwpm(
        med,
        N0,
        LAM,
        (60, 96, 96),
        1e-4,
        np.random.default_rng(1),
        keep_every=20,
        model="increment",
    )
    assert np.all(np.abs(cal.total - 1) < 0.02)  # as the index steps above
