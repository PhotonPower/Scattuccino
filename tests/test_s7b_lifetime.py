"""S7b: fluorescence lifetime in the time and frequency domain (T8b).

The emission is delayed by an exponential waiting time with lifetime τ, independent of the paths
(single-exponential decay, lakowicz2006): the detected time profile is the one for τ = 0
convolved with e^(−t/τ)/τ, i.e. in the frequency domain multiplied by 1/(1 + iωτ) (phase
tan φ = ωτ, modulation 1/√(1 + ω²τ²)). Excitation timing: absorption at depth z in a
non-scattering slab happens after the optical path n z.
"""

from __future__ import annotations

import numpy as np

from scattuccino.fluorescence import simulate_fluorescence
from scattuccino.layers import Layer, Stack
from scattuccino.transport import simulate

N = 1.4


def test_absorption_events_carry_the_optical_path() -> None:
    # matched index outside: no reflections, every packet absorbed on its first pass
    stack = Stack((Layer(1.0, 0.5, 0.0, 0.0, N),), N, N)
    sim = simulate(stack, 2000, seed=1, record_absorption=True)
    assert sim.events is not None
    # μs = 0: the first interaction absorbs all weight; later zero-weight events carry nothing
    live = sim.events["dw"] > 0
    assert live.sum() > 500
    assert np.allclose(sim.events["opl"][live], N * sim.events["z"][live], rtol=1e-12, atol=1e-12)


def _pair(lifetime_ps: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:  # type: ignore[type-arg]
    stack_x = Stack((Layer(np.inf, 0.05, 10.0, 0.8, N),))
    stack_m = Stack((Layer(np.inf, 0.02, 9.0, 0.8, N),))
    res = simulate_fluorescence(
        stack_x, stack_m, np.array([0.03]), 0.5, 20_000, 40_000, seed=4, lifetime_ps=lifetime_ps
    )
    p = res.emission.protocol
    top = p.side == -1
    return p.time_ps[top], p.weight[top], p.rho[top]


def test_lifetime_delays_each_exit_by_an_exponential_time() -> None:
    tau = 2000.0  # ps
    t0, w0, r0 = _pair(0.0)
    t1, w1, r1 = _pair(tau)
    # same excitation, sources and emission transport: only the delays differ
    assert np.array_equal(w0, w1) and np.array_equal(r0, r1)
    d = t1 - t0
    assert np.all(d >= 0)
    n = d.size
    assert abs(d.mean() / tau - 1) < 4 / np.sqrt(n)  # exponential: sd = mean
    assert abs(d.std() / tau - 1) < 8 / np.sqrt(n)
    # without lifetime the exits are later than the shortest straight path (excitation + emission)
    assert np.all(t0 > 0)


def test_frequency_domain_phase_and_modulation() -> None:
    # F_τ(ω)/F_0(ω) = 1/(1 + iωτ) with F(ω) = Σ w e^(−iωt); statistical error from the delays only
    tau = 2000.0
    t0, w, _ = _pair(0.0)
    t1, _, _ = _pair(tau)
    for omega_tau in (0.3, 1.0, 3.0):
        om = omega_tau / tau  # rad/ps
        f0 = np.sum(w * np.exp(-1j * om * t0))
        f1 = np.sum(w * np.exp(-1j * om * t1))
        ref = 1 / (1 + 1j * omega_tau)
        var = 1 - abs(ref) ** 2  # Var e^(−iωd) for an exponential d
        se = np.sqrt(np.sum(w * w) * var) / abs(f0)
        assert abs(f1 / f0 - ref) < 4 * se + 1e-12, omega_tau
        # phase and modulation
        ratio = f1 / f0
        assert abs(np.tan(-np.angle(ratio)) - omega_tau) < 4 * se * (1 + omega_tau**2) * 2
        assert abs(abs(ratio) - 1 / np.sqrt(1 + omega_tau**2)) < 4 * se
