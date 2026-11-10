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
from scattuccino.phase import PhaseFunction
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
    variances: tuple[float, float] | None = None,
) -> FloatArray:
    """Coupled power per packet: exit (x, y) at z = 0 with direction u (u_z < 0, inside the
    medium) traced back to the waist plane z_f; Wigner density of the Gaussian mode there."""
    k = 2 * np.pi * n0 / (wavelength_um * 1e-3)
    # the time-reversed ray goes down with direction −u: r(z_f) = r₀ + z_f (−u⊥)/(−u_z)
    xf = x + beam.z_f * u[:, 0] / u[:, 2]
    yf = y + beam.z_f * u[:, 1] / u[:, 2]
    kx, ky = -k * u[:, 0], -k * u[:, 1]
    # Wigner density of the Gaussian mode at its waist (as rtt_fusion.phasespace), or a general
    # Gaussian acceptance (multimode detector) with the same normalisation
    vx, vk = variances if variances is not None else (beam.w_f**2 / 4, 1 / beam.w_f**2)
    dens = np.exp(-(xf**2 + yf**2) / (2 * vx) - (kx**2 + ky**2) / (2 * vk)) / (
        (2 * np.pi) ** 2 * vx * vk
    )
    return np.asarray((2 * np.pi) ** 2 * weight * dens)


def next_event_signal(
    events: dict[str, FloatArray],
    mut: float,
    phase_pdf: object,
    beam: Beam,
    n0: float,
    wavelength_um: float,
    variances: tuple[float, float] | None = None,
) -> tuple[FloatArray, FloatArray]:
    """Next-event estimator of the mode-coupled power (S9b-4).

    At every scattering event (position r₀, depth z₀, weight w after absorption, incoming
    direction u) the expected power scattered into the detection mode and leaving the surface
    without further interaction, in the paraxial limit of the mode (its angular width ≪ the
    width of the phase function):

        E = w p(u·(−ẑ)) e^(−μt z₀) (2π/k²) exp(−ρ₀²/(2V))/V,  V = w_f²/4 + ((z_f − z₀)/k)²/w_f²

    – the Gaussian-mode integral over the scattering directions done analytically; V is the
    variance of the mode's intensity at depth z₀, so E ∝ |E_mode(r₀)|² (reciprocity). Matched
    index, homogeneous μt above the event. `variances` = (v_x, v_k) replaces the mode by a general
    Gaussian acceptance (V = v_x + ((z_f − z₀)/k)² v_k). Returns (contribution, optical path)."""
    k = 2 * np.pi * n0 / (wavelength_um * 1e-3)
    z0 = events["z"]
    vx, vk = variances if variances is not None else (beam.w_f**2 / 4, 1 / beam.w_f**2)
    v = vx + ((beam.z_f - z0) / k) ** 2 * vk
    rho2 = events["x"] ** 2 + events["y"] ** 2
    p = phase_pdf(-events["uz"])  # type: ignore[operator]
    e = events["w"] * p * np.exp(-mut * z0) * (2 * np.pi / k**2) * np.exp(-rho2 / (2 * v)) / v
    return np.asarray(e), np.asarray(events["opl"] + n0 * z0)


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
    estimator: str = "next_event",
    acceptance: tuple[float, float] | None = None,
    scattering: tuple[float, PhaseFunction] | None = None,
) -> HybridResult:
    """scattering: (μs, phase function) for the Monte Carlo below and on the way back above the
    hand-off; default Born. Use `Calibration.scattering()` from `calibrate_with_bwpm` on the same
    grid so that both regions describe the same medium (D9, S9b-2 finding)."""
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
    mus_mc, phase_mc = scattering if scattering is not None else (b.mus, b.phase)
    stack = Stack(
        (
            Layer(max(z_h, 1e-12), mua, mus_mc, 0.0, n0, phase=phase_mc),
            Layer(np.inf, mua, mus_mc, 0.0, n0, phase=phase_mc),
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
    if estimator not in ("next_event", "direct"):
        raise ValueError("estimator must be 'next_event' or 'direct'")
    nee = estimator == "next_event"
    sigma = l_c / (2 * np.sqrt(2 * np.log(2)))
    sig_z = np.zeros(z.size)

    def gate(values: FloatArray, depth: FloatArray) -> None:
        for i0 in range(0, depth.size, 50_000):
            d = (z[None, :] - depth[i0 : i0 + 50_000, None]) / sigma
            sig_z[:] += np.sum(values[i0 : i0 + 50_000, None] * np.exp(-0.5 * d * d), axis=0)

    def on_events(ev: dict[str, FloatArray]) -> None:
        val, opl = next_event_signal(
            ev, mua + mus_mc, phase_mc.pdf, beam, n0, wavelength_um, acceptance
        )
        keep = val > 0
        gate(val[keep], (n0 * z_h + opl[keep]) / 2)

    sim = simulate(
        stack,
        m,
        seed=seed + 7,
        sources=src,
        max_path=max_path,
        record_scattering=on_events if nee else False,
    )
    if not nee:
        pr = sim.protocol
        top = pr.side == -1
        u = np.column_stack([pr.ux, pr.uy, pr.uz])[top]
        coupled = mode_coupled_signal(
            pr.x[top], pr.y[top], u, pr.weight[top], beam, n0, wavelength_um, acceptance
        )
        gate(coupled, (n0 * z_h + (pr.path @ pr.n_layer)[top]) / 2)
    sig_z /= sigma * np.sqrt(2 * np.pi)
    return HybridResult(np.asarray(z), np.asarray(sig_z), p_h, m)
