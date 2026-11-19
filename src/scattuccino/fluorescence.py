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

Source sampling: N_m emission packets are drawn from the excitation events with probability
∝ fluorophore share, each with weight (total fluorophore absorption × Q)/N_m – unbiased, and the
number of emission packets is independent of the number of events.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from scattuccino._types import FloatArray
from scattuccino.layers import Stack
from scattuccino.transport import C_MM_PER_PS, Simulation, Sources, simulate


@dataclass
class FluorescenceResult:
    excitation: Simulation
    emission: Simulation
    absorbed_fluorophore: float  # per incident power
    emitted: float  # Q × absorbed_fluorophore


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
    total_f = float(w_f.sum())
    rng = np.random.default_rng(seed + 1_000_003)
    if total_f <= 0:
        pick = np.zeros(packets_m, dtype=np.int64)
        weight = np.zeros(packets_m)
    else:
        pick = rng.choice(w_f.size, size=packets_m, p=w_f / total_f)
        weight = np.full(packets_m, total_f * quantum_yield / packets_m)
    start = exc.events["opl"][pick] / C_MM_PER_PS
    if lifetime_ps > 0:
        start = start + np.random.default_rng(seed + 2_000_003).exponential(lifetime_ps, packets_m)
    src = Sources(
        exc.events["x"][pick],
        exc.events["y"][pick],
        exc.events["z"][pick],
        lay[pick],
        weight,
        time_ps=start,
    )
    emi = simulate(stack_m, packets_m, seed=seed + 2, sources=src)
    return FluorescenceResult(exc, emi, total_f, total_f * quantum_yield)
