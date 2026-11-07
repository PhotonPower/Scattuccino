"""S9b-3a: hybrid OCT chain (wave optics → hand-off → Monte Carlo → mode coupling) (T17a).

Mechanics only: power at the hand-off plane, the mode coupling of the exit packets and the
coherence-gated A-scan exist and are finite. The statistical estimator of the mode-coupled signal
is too noisy for quantitative comparisons (a few packets carry the signal) – a next-event
estimator towards the detection mode follows in S9b-4 (docs/testfaelle.md, T17).
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("rtt_fusion")
pytest.importorskip("bwpm")

from scattuccino.hybrid import (  # noqa: E402
    Beam,
    gaussian_beam_field,
    hybrid_oct,
    mode_coupled_signal,
)
from scattuccino.microstructure import WhittleMatern  # noqa: E402


def test_gaussian_beam_focuses_at_z_f() -> None:
    from rtt_fusion.ports import GriddedField, Port

    lam, n0, w_f, z_f, nx, dx = 0.8, 1.38, 2e-3, 0.01, 128, 2.5e-4
    e = gaussian_beam_field(Beam(w_f, z_f), n0, lam, nx, dx)
    k = 2 * np.pi * n0 / (lam * 1e-3)
    spec = np.fft.fft2(e[:, :, 0])
    kk = 2 * np.pi * np.fft.fftfreq(nx, dx)
    kx, ky = np.meshgrid(kk, kk)
    field_f = np.fft.ifft2(spec * np.exp(1j * np.sqrt(k**2 - kx**2 - ky**2 + 0j) * z_f))
    x = (np.arange(nx) - nx // 2) * dx
    xx, _ = np.meshgrid(x, x)
    inten = np.abs(field_f) ** 2
    assert abs(np.sqrt(np.sum(inten * xx**2) / inten.sum()) / (w_f / 2) - 1) < 0.02  # rms = w/2
    assert GriddedField(Port(0.0, 1, complex(n0), lam), dx, dx, e).field.shape == (nx, nx, 2)


def test_mode_coupling_of_a_ray_through_the_focus() -> None:
    # a ray leaving the surface that traces back exactly through the waist centre along the axis
    # couples with the peak density; a displaced one with exp(−2 d²/w²)·…
    beam, n0, lam = Beam(2e-3, 0.01), 1.38, 0.8
    u = np.array([[0.0, 0.0, -1.0], [0.0, 0.0, -1.0]])
    w = np.array([1.0, 1.0])
    c = mode_coupled_signal(np.array([0.0, 1e-3]), np.array([0.0, 0.0]), u, w, beam, n0, lam)
    assert (
        abs(c[0] - (2 * np.pi) ** 2 / ((2 * np.pi) ** 2 * (beam.w_f**2 / 4) * (1 / beam.w_f**2)))
        < 1e-9
    )
    assert abs(c[1] / c[0] - np.exp(-((1e-3) ** 2) / (2 * beam.w_f**2 / 4))) < 1e-12


def test_hybrid_chain_runs_and_conserves_power() -> None:
    med = WhittleMatern(0.0572, 0.5e-3)
    z = np.linspace(0.0, 0.03, 31).astype(np.float64)
    r0 = hybrid_oct(
        med,
        1.38,
        1.0,
        0.8,
        Beam(2e-3, 0.004),
        0.0,
        z,
        0.004,
        20_000,
        seed=3,
        nx=96,
        dx=1.5e-4,
        max_path=0.1,
        sigma_husimi=0.6e-3,
    )
    r1 = hybrid_oct(
        med,
        1.38,
        1.0,
        0.8,
        Beam(2e-3, 0.004),
        0.003,
        z,
        0.004,
        20_000,
        seed=3,
        nx=96,
        dx=1.5e-4,
        max_path=0.1,
        sigma_husimi=0.6e-3,
    )
    assert 0.9 < r0.power_handoff <= 1.0 + 1e-12  # a few per cent evanescent in the Husimi window
    assert 0.95 < r1.power_handoff < 1.02  # forward FPWPM through 3 µm of the random medium
    for r in (r0, r1):
        assert np.all(np.isfinite(r.signal)) and np.all(r.signal >= 0) and r.signal.sum() > 0
