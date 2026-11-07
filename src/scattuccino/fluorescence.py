"""Fluorescence and Raman (stage S7a, Python reference): inelastic emission in layered media.

Two-stage Monte Carlo:

1. **Excitation** at λ_x: the stack carries the total absorption μa = μa,tissue + μa,f per layer;
   every absorption event (implicit capture) gives the fluorophore its share μa,f/μa of the
   absorbed weight.
2. **Emission** at λ_m: isotropic packets from these events, weight × quantum yield Q (Raman:
   the scattering efficiency into the band), transported with the optical properties at λ_m (the
   fluorophore's absorption at λ_m included, i.e. re-absorption). Re-absorbed emission is not
   re-emitted (single generation; documented limitation).

Source sampling: N_m emission packets are drawn from the excitation events with probability
∝ fluorophore share, each with weight (total fluorophore absorption × Q)/N_m – unbiased, and the
number of emission packets is independent of the number of events.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from scattuccino._types import FloatArray
from scattuccino.layers import Stack
from scattuccino.transport import Simulation, Sources, simulate


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
) -> FluorescenceResult:
    """stack_x: optical properties at λ_x including the fluorophore (μa,f,x per layer in mua_f_x);
    stack_m: properties at λ_m (including re-absorption by the fluorophore)."""
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
    src = Sources(
        exc.events["x"][pick], exc.events["y"][pick], exc.events["z"][pick], lay[pick], weight
    )
    emi = simulate(stack_m, packets_m, seed=seed + 2, sources=src)
    return FluorescenceResult(exc, emi, total_f, total_f * quantum_yield)
