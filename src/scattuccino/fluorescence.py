"""Fluorescence and Raman (stage S7a, Python reference): inelastic emission in layered media.

Two-stage Monte Carlo:

1. **Excitation** at λ_x: the stack carries the total absorption μa = μa,tissue + μa,f per layer;
   every absorption event (implicit capture) gives the fluorophore its share μa,f/μa of the
   absorbed weight.
2. **Emission** at λ_m: isotropic packets from these events, weight × quantum yield Q (Raman:
   the scattering efficiency into the band), transported with the optical properties at λ_m (the
   fluorophore's absorption at λ_m included, i.e. re-absorption). Re-absorbed emission is not
   re-emitted (single generation; documented limitation).
3. **Time** (S7b): an emission packet starts at the time of its excitation event (optical path of
   the excitation packet / c) plus an exponential delay with the fluorescence lifetime τ
   (single-exponential decay, lakowicz2006); `Protocol.time_ps` of the emission includes both.
   The delays come from their own random stream, so the transport is the same for every τ.
4. **Emission spectrum** (S7c): the emission is split into wavelength bins with the normalised
   spectrum f_b (Σ f_b = 1); bin b carries the weight Q f_b × fluorophore absorption and is
   transported with its own optical properties (re-absorption per wavelength: inner filter). One
   excitation serves all bins; the packets are shared out ∝ f_b.

Source sampling: N_m emission packets are drawn from the excitation events with probability
∝ fluorophore share, each with weight (total fluorophore absorption × Q)/N_m – unbiased, and the
number of emission packets is independent of the number of events.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from scattuccino._types import FloatArray, IntArray
from scattuccino.layers import Stack
from scattuccino.transport import C_MM_PER_PS, Simulation, Sources, simulate


@dataclass
class FluorescenceResult:
    excitation: Simulation
    emission: Simulation
    absorbed_fluorophore: float  # per incident power
    emitted: float  # Q × absorbed_fluorophore


def _excite(
    stack_x: Stack, mua_f_x: FloatArray, packets_x: int, seed: int, cos_incident: float
) -> tuple[Simulation, IntArray, FloatArray, float]:
    """Excitation run; returns it, the layer of each absorption event, the fluorophore's share
    of the absorbed weight per event and its total (per incident power)."""
    mua_f = np.asarray(mua_f_x, dtype=np.float64)
    mua_x = np.array([layer.mua for layer in stack_x.layers])
    if np.any(mua_f > mua_x + 1e-15):
        raise ValueError("fluorophore absorption must be part of the layer absorption at λ_x")
    exc = simulate(stack_x, packets_x, seed=seed, cos_incident=cos_incident, record_absorption=True)
    assert exc.events is not None
    lay = exc.events["layer"].astype(np.int64)
    with np.errstate(invalid="ignore", divide="ignore"):
        share = np.where(mua_x[lay] > 0, mua_f[lay] / mua_x[lay], 0.0)
    w_f = exc.events["dw"] * share
    return exc, lay, w_f, float(w_f.sum())


def _emit(
    stack_m: Stack,
    exc: Simulation,
    lay: IntArray,
    w_f: FloatArray,
    emitted: float,
    packets_m: int,
    seed_pick: int,
    seed_transport: int,
    seed_delay: int,
    lifetime_ps: float,
) -> Simulation:
    """Emission run: packets_m isotropic sources drawn from the excitation events ∝ w_f, each with
    weight emitted/packets_m, started at the excitation time plus the lifetime delay."""
    assert exc.events is not None
    total_f = float(w_f.sum())
    rng = np.random.default_rng(seed_pick)
    pick: IntArray
    if total_f <= 0:
        pick = np.zeros(packets_m, dtype=np.int64)
        weight = np.zeros(packets_m)
    else:
        pick = rng.choice(w_f.size, size=packets_m, p=w_f / total_f)
        weight = np.full(packets_m, emitted / packets_m)
    start = exc.events["opl"][pick] / C_MM_PER_PS
    if lifetime_ps > 0:
        start = start + np.random.default_rng(seed_delay).exponential(lifetime_ps, packets_m)
    src = Sources(
        exc.events["x"][pick],
        exc.events["y"][pick],
        exc.events["z"][pick],
        lay[pick],
        weight,
        time_ps=start,
    )
    return simulate(stack_m, packets_m, seed=seed_transport, sources=src)


def simulate_fluorescence(
    stack_x: Stack,
    stack_m: Stack,
    mua_f_x: FloatArray,
    quantum_yield: float,
    packets_x: int,
    packets_m: int,
    seed: int = 1,
    cos_incident: float = 1.0,
    lifetime_ps: float = 0.0,
) -> FluorescenceResult:
    """stack_x: optical properties at λ_x including the fluorophore (μa,f,x per layer in mua_f_x);
    stack_m: properties at λ_m (including re-absorption by the fluorophore); lifetime_ps:
    fluorescence lifetime τ (ps, 0 = instantaneous; Raman: 0)."""
    exc, lay, w_f, total_f = _excite(stack_x, mua_f_x, packets_x, seed, cos_incident)
    emitted = total_f * quantum_yield
    emi = _emit(
        stack_m,
        exc,
        lay,
        w_f,
        emitted,
        packets_m,
        seed + 1_000_003,
        seed + 2,
        seed + 2_000_003,
        lifetime_ps,
    )
    return FluorescenceResult(exc, emi, total_f, emitted)


@dataclass
class SpectralFluorescenceResult:
    excitation: Simulation
    emission: list[Simulation]  # one per wavelength bin
    absorbed_fluorophore: float  # per incident power
    emitted: FloatArray  # Q f_b × absorbed_fluorophore per bin
    packets: IntArray  # emission packets per bin

    def detected(self, side: int = -1) -> FloatArray:
        """Escaping fluorescence per bin through `side` (per incident power)."""
        return np.array(
            [float(e.protocol.weight[e.protocol.side == side].sum()) for e in self.emission]
        )


def simulate_fluorescence_spectrum(
    stack_x: Stack,
    stacks_m: Sequence[Stack],
    spectrum: FloatArray,
    mua_f_x: FloatArray,
    quantum_yield: float,
    packets_x: int,
    packets_m: int,
    seed: int = 1,
    cos_incident: float = 1.0,
    lifetime_ps: float = 0.0,
) -> SpectralFluorescenceResult:
    """Fluorescence with an emission spectrum (S7c): stacks_m[b] holds the optical properties at
    the wavelength of bin b (including the fluorophore's re-absorption there), spectrum[b] the
    emission in bin b (normalised here to Σ = 1); packets_m emission packets in total, shared out
    ∝ spectrum (at least one per bin)."""
    f_b = np.asarray(spectrum, dtype=np.float64)
    if f_b.ndim != 1 or f_b.size != len(stacks_m) or np.any(f_b < 0) or f_b.sum() <= 0:
        raise ValueError("spectrum needs one non-negative value per stack, not all zero")
    f_b = f_b / f_b.sum()
    exc, lay, w_f, total_f = _excite(stack_x, mua_f_x, packets_x, seed, cos_incident)
    n_b = np.maximum(np.rint(packets_m * f_b).astype(np.int64), 1)
    emitted = total_f * quantum_yield * f_b
    emission = [
        _emit(
            stacks_m[b],
            exc,
            lay,
            w_f,
            float(emitted[b]),
            int(n_b[b]),
            seed + 1_000_003 + 10_007 * (b + 1),
            seed + 3 + 10_007 * (b + 1),
            seed + 2_000_003 + 10_007 * (b + 1),
            lifetime_ps,
        )
        for b in range(f_b.size)
    ]
    return SpectralFluorescenceResult(exc, emission, total_f, emitted, n_b)
