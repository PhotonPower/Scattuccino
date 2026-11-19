"""S9c-4a: back scattering of BWPM at a random δn against first-order Born (T10c).

`microstructure.born_reflection` (first Born for δn constant within each slice, Weyl expansion)
is checked against the exact reflection of a weak homogeneous layer (Airy formula); then the
reflected field of BWPM (B8: back scattering at `fluctuation`, iterations = 1) for the same
realisation against it: the power agrees, the field up to the forward phase modulation on the
way in and out, which is first order in δn and not part of first-order Born.
"""

from __future__ import annotations

import numpy as np
import pytest

from scattuccino.microstructure import WhittleMatern, born_reflection, random_field

N0, LAM = 1.38, 0.8


def test_born_reflection_of_a_homogeneous_layer() -> None:
    # Airy: r = (r₁₂ + r₂₃ e^(2ik₁d))/(1 + r₁₂ r₂₃ e^(2ik₁d)), r₂₃ = −r₁₂, k₁ = n₁k₀
    k0 = 2 * np.pi / (LAM * 1e-3)
    delta, h, n = 1e-4, 1e-4, 8
    for nz in (1, 3, 7, 12, 30):
        d = nz * h
        f = born_reflection(np.full((nz, n, n), delta), N0, LAM, h)
        r_born = f[0, 0, 0] / n**2
        n1 = N0 + delta
        r12 = (N0 - n1) / (N0 + n1)
        ph = np.exp(2j * n1 * k0 * d)
        r_exact = (r12 - r12 * ph) / (1 - r12 * r12 * ph)
        # first order: Born keeps the host phase k₀n₀ in the layer, an error ≈ 2δk₀d relative to the
        # reflection scale δ/n₀ (absolute: r itself vanishes where 2k d ≈ 2π, e.g. nz = 3)
        assert abs(r_born - r_exact) < (1e-4 + 3 * delta * k0 * d) * delta / N0, nz
        assert np.allclose(f[..., 0].ravel()[1:], 0, atol=1e-12 * n**2)  # uniform: specular only
        assert np.allclose(f[..., 1], 0, atol=1e-12 * n**2)


def test_bwpm_back_scattering_against_born() -> None:
    fluct = pytest.importorskip("bwpm.fluctuation")
    if not hasattr(fluct, "born_emission"):
        pytest.skip("BWPM without back scattering at δn (B8)")
    from bwpm import fast
    from bwpm.volume3d import Bidirectional3D, Grid3D

    h, n = 1e-4, 64
    k = 2 * np.pi * N0 / (LAM * 1e-3)
    kx = 2 * np.pi * np.fft.fftfreq(n, h)
    kxx, kyy = np.meshgrid(kx, kx)
    use = kxx**2 + kyy**2 < (0.96 * k) ** 2  # |cos θ| > 0.28: away from grazing cells
    err = {}
    for sigma in (0.0036, 0.0143):
        dn = random_field(
            WhittleMatern(sigma, 0.5e-3, 4.0), (40, n, n), h, np.random.default_rng(1)
        )
        solver = Bidirectional3D(
            Grid3D(n, n, h, h, h, LAM),
            [complex(N0)],
            np.zeros(dn.shape, dtype=np.int64),
            complex(N0),
            None,
            fluctuation=dn,
        )
        e = np.zeros((n, n, 2), dtype=np.complex128)
        e[:, :, 0] = 1.0
        rb = np.fft.fft2(fast.run_3d(solver, e, iterations=1).reflected, axes=(0, 1))
        ref = born_reflection(dn, N0, LAM, h)
        p_ref = np.sum(np.abs(ref[use]) ** 2)
        err[sigma] = float(np.sqrt(np.sum(np.abs(rb[use] - ref[use]) ** 2) / p_ref))
        if sigma == 0.0036:  # Born regime: power measured 1.019 (seeds 1–3: 0.996–1.019)
            assert abs(np.sum(np.abs(rb[use]) ** 2) / p_ref - 1) < 0.05
    # the field differs by the forward phase modulation, first order in σ_n: measured 0.079 and
    # 0.311, ratio 3.94 (×2 in σ_n: 1.99–2.01)
    assert err[0.0036] < 0.15
    assert 3.5 < err[0.0143] / err[0.0036] < 4.5
