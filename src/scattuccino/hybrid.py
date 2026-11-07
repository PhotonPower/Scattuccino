"""Hybrid OCT: wave optics above the hand-off plane, photon transport below (stage S9b-3).

Decisions D8–D10 (`docs/entwurf-s9b-hybrid.md`):

1. A Gaussian beam focused into the tissue (waist w_f at depth z_f, medium index n₀) enters a
   BWPM volume holding one realisation of the microstructure (`microstructure.random_field`)
   from z = 0 to the hand-off depth z_h (forward FPWPM; the reflection of this region is not part
   of the depth-resolved signal below z_h).
2. The field at z_h becomes photon packets through its Husimi distribution
   (`rtt_fusion.phasespace.field_to_packets`); z_h = 0 is the pure Monte-Carlo limit (packets from
   the incident field at the surface).
3. Below z_h – and on the way back up through 0 < z < z_h – the packets are transported by Monte
   Carlo with the statistical parameters of the same tissue (μs, phase function from Born).
4. Packets leaving through the top couple into the detection mode (the same focused beam,
   reciprocity) with the Wigner density of the mode at its waist: the exit ray is traced back to
   z_f, coupled power = (2π)² Σ wᵢ p_m(r_f,i, k⊥,i) (`rtt_fusion.phasespace`).
5. Depth z_OPL = (n₀ z_h + OPL_MC)/2 with coherence gating as S9a.

Index-matched surface (n_above = n₀) in this stage: the mode is defined inside the medium.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from scattuccino._types import FloatArray
from scattuccino.layers import Layer, Stack
from scattuccino.microstructure import WhittleMatern, born, quantise, random_field
from scattuccino.transport import Sources, simulate


@dataclass(frozen=True)
class Beam:
    w_f: float  # waist radius (1/e² intensity), mm
    z_f: float  # focus depth below the surface, mm


def gaussian_beam_field(
    beam: Beam, n0: float, wavelength_um: float, nx: int, dx: float
) -> np.ndarray:  # type: ignore[type-arg]
    """(ny, nx, 2) field at z = 0 of a Gaussian beam (x polarised) with its waist at z_f, in n₀."""
    k = 2 * np.pi * n0 / (wavelength_um * 1e-3)
    z_r = k * beam.w_f**2 / 2
    q0 = -beam.z_f - 1j * z_r  # q(z) = z − z_f − i z_R for e^(i(kz − ωt)): waist at z_f
    x = (np.arange(nx) - nx // 2) * dx
    xx, yy = np.meshgrid(x, x)
    e = np.zeros((nx, nx, 2), dtype=np.complex128)
    e[:, :, 0] = (-1j * z_r / q0) * np.exp(1j * k * (xx**2 + yy**2) / (2 * q0))
    return e


def mode_coupled_signal(
    x: FloatArray,
    y: FloatArray,
    u: FloatArray,
    weight: FloatArray,
    beam: Beam,
    n0: float,
    wavelength_um: float,
) -> FloatArray:
    """Coupled power per packet: exit (x, y) at z = 0 with direction u (u_z < 0, inside the
    medium) traced back to the waist plane z_f; Wigner density of the Gaussian mode there."""
    from rtt_fusion.phasespace import gaussian_mode_density

    k = 2 * np.pi * n0 / (wavelength_um * 1e-3)
    # the time-reversed ray goes down with direction −u: r(z_f) = r₀ + z_f (−u⊥)/(−u_z)
    xf = x + beam.z_f * u[:, 0] / u[:, 2]
    yf = y + beam.z_f * u[:, 1] / u[:, 2]
    kx, ky = -k * u[:, 0], -k * u[:, 1]
    return np.asarray((2 * np.pi) ** 2 * weight * gaussian_mode_density(xf, yf, kx, ky, beam.w_f))


@dataclass
class HybridResult:
    z: FloatArray  # depth grid of the A-scan (mm)
    signal: FloatArray  # mode-coupled, coherence-gated OCT signal per mm
    power_handoff: float  # forward power at z_h (incident power 1)
    packets: int


def hybrid_oct(
    medium: WhittleMatern,
    n0: float,
    mua: float,
    wavelength_um: float,
    beam: Beam,
    z_h: float,
    z: FloatArray,
    l_c: float,
    packets: int,
    seed: int = 1,
    nx: int = 128,
    dx: float = 1e-4,
    sigma_husimi: float | None = None,
    backend: str = "auto",
    max_path: float = np.inf,
) -> HybridResult:
    from rtt_fusion.phasespace import field_to_packets
    from rtt_fusion.ports import GriddedField, Port, power

    rng = np.random.default_rng(seed)
    e0 = gaussian_beam_field(beam, n0, wavelength_um, nx, dx)
    port0 = Port(0.0, 1, complex(n0), wavelength_um)
    p_in = power(GriddedField(port0, dx, dx, e0))
    if z_h > 0:
        from bwpm import fast
        from bwpm.volume3d import Bidirectional3D, Grid3D

        nz = round(z_h / dx)
        dn = random_field(medium, (nz, nx, nx), dx, rng)
        lab, idx = quantise(dn, n0)
        solver = Bidirectional3D(
            Grid3D(nx, nx, dx, dx, dx, wavelength_um), idx, lab, complex(n0), None
        )
        res = fast.run_3d(solver, e0, iterations=0, keep=(nz,), backend=backend)
        e_h = res.planes[nz][0]
    else:
        e_h = e0
    g_h = GriddedField(Port(z_h, 1, complex(n0), wavelength_um), dx, dx, e_h)
    sig = sigma_husimi if sigma_husimi is not None else 0.5 * beam.w_f
    pk = field_to_packets(g_h, packets, sig, rng, periodic=True)
    p_h = float(pk.weight.sum()) / p_in
    b = born(medium, n0, wavelength_um)
    stack = Stack(
        (
            Layer(max(z_h, 1e-12), mua, b.mus, 0.0, n0, phase=b.phase),
            Layer(np.inf, mua, b.mus, 0.0, n0, phase=b.phase),
        ),
        n0,
        n0,
    )
    m = pk.x.size
    src = Sources(
        np.asarray(pk.x).astype(np.float64),
        np.asarray(pk.y).astype(np.float64),
        np.full(m, max(z_h, 1e-12)),
        np.ones(m, dtype=np.int64),
        np.asarray(pk.weight / p_in).astype(np.float64),
        np.column_stack([pk.ux, pk.uy, pk.uz]).astype(np.float64),
    )
    sim = simulate(stack, m, seed=seed + 7, sources=src, max_path=max_path)
    pr = sim.protocol
    top = pr.side == -1
    u = np.column_stack([pr.ux, pr.uy, pr.uz])[top]
    coupled = mode_coupled_signal(pr.x[top], pr.y[top], u, pr.weight[top], beam, n0, wavelength_um)
    zo = (n0 * z_h + (pr.path @ pr.n_layer)[top]) / 2
    sigma = l_c / (2 * np.sqrt(2 * np.log(2)))
    d = (z[None, :] - zo[:, None]) / sigma
    sig_z = np.sum(coupled[:, None] * np.exp(-0.5 * d * d), axis=0) / (sigma * np.sqrt(2 * np.pi))
    return HybridResult(np.asarray(z), np.asarray(sig_z), p_h, m)
