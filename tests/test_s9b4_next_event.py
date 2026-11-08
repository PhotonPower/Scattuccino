"""S9b-4: next-event estimator for the coupling into a detection mode (T17b)."""

from __future__ import annotations

import numpy as np
import pytest

from scattuccino._types import FloatArray
from scattuccino.hybrid import Beam, mode_coupled_signal, next_event_signal
from scattuccino.layers import Layer, Stack
from scattuccino.microstructure import WhittleMatern, born
from scattuccino.transport import simulate

BEAM = Beam(2e-3, 0.004)
ACCEPT = ((10e-3) ** 2, 1500.0**2)  # multimode acceptance: 10 µm, k-spread 1500/mm (NA ≈ 0.14)


@pytest.fixture(scope="module")
def runs() -> list[tuple[float, float]]:
    b = born(WhittleMatern(0.0572, 0.5e-3), 1.38, 0.8)
    stack = Stack((Layer(np.inf, 1.0, b.mus, 0.0, 1.38, phase=b.phase),), 1.38, 1.38)
    out = []
    for seed in (1, 2, 3):
        sim = simulate(stack, 60_000, seed=seed, max_path=0.1, record_scattering=True)
        p = sim.protocol
        top = p.side == -1
        u = np.column_stack([p.ux, p.uy, p.uz])[top]
        direct = mode_coupled_signal(
            p.x[top], p.y[top], u, p.weight[top], BEAM, 1.38, 0.8, ACCEPT
        ).sum()
        assert sim.scatter_events is not None
        nee = next_event_signal(
            sim.scatter_events, 1.0 + b.mus, b.phase.pdf, BEAM, 1.38, 0.8, ACCEPT
        )[0].sum()
        out.append((float(direct), float(nee)))
    return out


def test_next_event_equals_direct_estimator(runs: list[tuple[float, float]]) -> None:
    direct = np.mean([r[0] for r in runs])
    nee = np.mean([r[1] for r in runs])
    assert (
        abs(nee / direct - 1) < 0.08
    )  # paraxial limit of the acceptance (NA 0.14): measured 3–4 %


def test_callback_receives_the_same_events() -> None:
    stack = Stack((Layer(1.0, 0.2, 8.0, 0.8, 1.4),))
    stored = simulate(stack, 3_000, seed=4, record_scattering=True).scatter_events
    batches: list[dict[str, FloatArray]] = []
    simulate(stack, 3_000, seed=4, record_scattering=batches.append)
    assert stored is not None
    for key in ("x", "z", "w", "uz", "opl"):
        assert np.array_equal(np.concatenate([b[key] for b in batches]), stored[key])


def test_coupling_of_a_point_scatterer_follows_the_mode_intensity() -> None:
    # one event on the axis at the waist vs. 1 µm off axis: ratio exp(−ρ²/(2 w²/4))
    zero: FloatArray = np.zeros(2)
    ev: dict[str, FloatArray] = {
        "y": zero,
        "ux": zero,
        "uy": zero,
        "opl": zero,
        "layer": zero,
        "x": np.asarray([0.0, 1e-3]),
        "z": np.full(2, BEAM.z_f),
        "w": np.ones(2),
        "uz": np.ones(2),
    }
    val, _ = next_event_signal(ev, 0.0, lambda mu: np.ones_like(mu), BEAM, 1.38, 0.8)
    assert abs(val[1] / val[0] - np.exp(-((1e-3) ** 2) / (2 * BEAM.w_f**2 / 4))) < 1e-12
