"""S9b-6: coherent hand-off without Husimi bias (T18).

The field at z_h is split into its projection onto the freely propagated incident beam (coherent,
packets from the exact Wigner function of the Gaussian beam) and the orthogonal rest (scattered,
Husimi packets). Reference values are analytic: second moments of a Gaussian beam, the overlap of
two displaced Gaussians, and the single-scattering OCT signal of a focused beam.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("rtt_fusion")

from scattuccino._types import FloatArray  # noqa: E402
from scattuccino.hybrid import (  # noqa: E402
    Beam,
    free_field,
    gaussian_beam_field,
    gaussian_beam_packets,
    hybrid_oct,
    split_coherent,
)
from scattuccino.microstructure import WhittleMatern  # noqa: E402
from scattuccino.phase import HenyeyGreenstein  # noqa: E402

N0, LAM = 1.38, 0.8
K = 2 * np.pi * N0 / (LAM * 1e-3)  # 1/mm


def test_wigner_packets_have_the_moments_of_the_beam() -> None:
    # paraxial Gaussian beam: Var x = w(z)²/4 with w(z)² = w_f² (1 + (Δz/z_R)²), z_R = k w_f²/2;
    # Var k_x = 1/w_f²; Cov(x, u_x) = Δz Var(u_x) (rays from the waist)
    beam = Beam(2e-3, 0.012)
    n = 400_000
    z_r = K * beam.w_f**2 / 2
    for z_p in (0.0, 0.012, 0.03):
        x, y, u = gaussian_beam_packets(beam, N0, LAM, z_p, n, np.random.default_rng(5))
        dz = z_p - beam.z_f
        var_x = beam.w_f**2 * (1 + (dz / z_r) ** 2) / 4
        assert abs(np.var(x) / var_x - 1) < 0.015  # 5σ of the sample variance: 1.1 %
        assert abs(np.var(y) / var_x - 1) < 0.015
        assert abs(np.var(u[:, 0] * K) * beam.w_f**2 - 1) < 0.015
        cov = float(np.mean(x * u[:, 0]) / (dz * np.var(u[:, 0]))) if dz != 0 else 1.0
        assert abs(cov - 1) < 0.02
        assert np.allclose(np.linalg.norm(u, axis=1), 1.0) and np.all(u[:, 2] > 0)


def test_free_field_focuses_the_beam() -> None:
    # the beam field at z = 0 propagated to z_f: rms radius w_f/2 (as T17a, through free_field)
    beam, nx, dx = Beam(2e-3, 0.01), 128, 2.5e-4
    e = free_field(gaussian_beam_field(beam, N0, LAM, nx, dx), N0, LAM, dx, beam.z_f)
    x = (np.arange(nx) - nx // 2) * dx
    xx, _ = np.meshgrid(x, x)
    inten = np.sum(np.abs(e) ** 2, axis=2)
    assert abs(np.sqrt(np.sum(inten * xx**2) / inten.sum()) / (beam.w_f / 2) - 1) < 0.02


def test_projection_splits_the_power_exactly() -> None:
    # beam at its waist (z_f = 0) and a copy displaced by Δ: |c|² = exp(−Δ²/w²) (overlap of two
    # Gaussians E ∝ exp(−r²/w²)); the rest is orthogonal, the power splits without cross term
    beam, nx, dx = Beam(2e-3, 0.0), 128, 1.5e-4
    e_g = gaussian_beam_field(beam, N0, LAM, nx, dx)
    shift = 10
    e = 0.7 * e_g + 0.4 * np.roll(e_g, shift, axis=1)
    c, rest = split_coherent(e, e_g)
    overlap = np.exp(-((shift * dx) ** 2) / (2 * beam.w_f**2))
    assert abs(c - (0.7 + 0.4 * overlap)) < 1e-9

    def p(f: np.ndarray) -> float:  # type: ignore[type-arg]
        return float(np.sum(np.abs(f) ** 2))

    assert abs(p(e) - (abs(c) ** 2 * p(e_g) + p(rest))) < 1e-12 * p(e)
    assert abs(np.vdot(e_g, rest)) < 1e-12 * p(e_g)


def single_scattering(
    d: FloatArray, beam: Beam, mus: float, mua: float, p_back: float, l_c: float
) -> FloatArray:
    """First-order OCT signal per optical depth d = n₀ z (index-matched surface, paraxial): events
    μs e^(−μt z) per depth with the beam's intensity density exp(−ρ²/2V)/(2πV), each coupling
    w p(−1) e^(−μt z) (2π/k²) exp(−ρ²/2V)/V into the mode; ∫ d²ρ gives 1/(2V):
    S₁(z) = μs p(−1) (2π/k²) e^(−2μt z)/(2V(z)), V = w_f²/4 + ((z_f − z)/k)²/w_f², per physical
    depth; the Gaussian coherence gate (FWHM l_c) in the optical depth n₀z, integrated over z,
    gives the signal per optical depth (Jacobian 1/n₀ implicit)."""
    mut = mus + mua
    zz = np.linspace(0.0, 0.2, 40_001)  # physical depth of the event
    v = beam.w_f**2 / 4 + ((beam.z_f - zz) / K) ** 2 / beam.w_f**2
    s1 = mus * p_back * (2 * np.pi / K**2) * np.exp(-2 * mut * zz) / (2 * v)
    sig = l_c / (2 * np.sqrt(2 * np.log(2)))
    g = np.exp(-0.5 * ((d[:, None] - N0 * zz[None, :]) / sig) ** 2) / (sig * np.sqrt(2 * np.pi))
    return np.asarray(np.trapezoid(s1[None, :] * g, zz, axis=1))


def test_single_scattering_in_a_tight_focus() -> None:
    # z_h = 0, weakly scattering (μs z ≤ 0.04 over the window: higher orders ≲ 2 %): the hybrid
    # signal is the analytic single-scattering signal of the focused beam. Isotropic scattering:
    # the reference needs only p(−1), and a forward peak would let rare paths (backscattered,
    # then forward into the mode, p(1)/p(−1) ≈ 7000 for g = 0.9) dominate the error.
    # Measured 0.994–1.016 ± 0.4–2 %; with coherent="husimi" (σ = 0.6 µm) 0.85 → 0.37.
    beam, mus, l_c = Beam(2e-3, 0.012), 1.0, 0.004
    hg = HenyeyGreenstein(0.0)
    d = np.linspace(0.005, 0.05, 10).astype(np.float64)
    runs = [
        hybrid_oct(
            WhittleMatern(0.0572, 0.5e-3),
            N0,
            0.0,
            LAM,
            beam,
            0.0,
            d,
            l_c,
            400_000,
            seed=s,
            nx=128,
            dx=1.5e-4,
            max_path=0.12,
            scattering=(mus, hg),
        ).signal
        for s in (1, 2, 3, 4, 5, 6)
    ]
    mean = np.mean(runs, axis=0)
    err = np.std(runs, axis=0, ddof=1) / np.sqrt(len(runs))
    ref = single_scattering(d, beam, mus, 0.0, float(hg.pdf(np.array(-1.0))), l_c)
    assert np.all(np.abs(mean - ref) < 4 * err + 0.02 * ref)


def test_coherent_share_behind_the_random_medium() -> None:
    # z_h = 3 µm of the Whittle–Matérn medium (Born μs = 200 /mm): the coherent share is the
    # ballistic power e^(−μs z_h) = 0.55 within the BWPM calibration's deviation near the surface
    # (coherent power 0.63 at 3 µm, issue 29: fit μs 0.95 × Born, curvature at small depth)
    pytest.importorskip("bwpm")
    from scattuccino.microstructure import born

    med = WhittleMatern(0.0572, 0.5e-3)
    z = np.linspace(0.0, 0.03, 16).astype(np.float64)
    r = hybrid_oct(
        med, N0, 1.0, LAM, Beam(2e-3, 0.004), 0.003, z, 0.004, 20_000, seed=3, nx=96, dx=1.5e-4,
        max_path=0.1, sigma_husimi=0.4e-3,
    )  # fmt: skip
    ballistic = np.exp(-born(med, N0, LAM).mus * 0.003)
    assert abs(r.power_coherent / ballistic - 1) < 0.25
    assert 0 < r.power_coherent < r.power_handoff < 1.02
    assert np.all(np.isfinite(r.signal)) and r.signal.sum() > 0
