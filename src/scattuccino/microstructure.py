"""Tissue microstructure as a random refractive-index field (stage S9b-2, decision D9).

δn(r) is a Gaussian random field with the Whittle–Matérn correlation (Rogers, Backman et al.):
power spectral density (3D, Φ(q) = (2π)⁻³ ∫ C(r) e^(−iq·r) d³r)

    Φ(q) = σ_n² l_c³ Γ(D/2) / (π^(3/2) |Γ((D − 3)/2)| (1 + q² l_c²)^(D/2)),   3 < D < 6,

D = 4 being the exponential correlation C(r) = σ_n² e^(−r/l_c). First-order Born approximation
for a scalar wave in n₀ + δn (k = n₀k₀): differential scattering coefficient per unit volume and
solid angle σ(θ) = (2π k⁴/n₀²) Φ(2k sin(θ/2)), μs = ∫ σ dΩ, p(θ) = σ/μs.

The field is generated on a periodic grid by filtering white noise with √Φ (exact second-order
statistics on the grid) and handed to BWPM as labels of a quantised index.

Units: mm (lengths, l_c, wavenumbers in 1/mm), µm for the wavelength.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import gamma

from scattuccino._types import FloatArray, IntArray
from scattuccino.phase import Tabulated


@dataclass(frozen=True)
class WhittleMatern:
    sigma_n: float  # rms of δn
    l_c: float  # correlation length, mm
    d: float = 4.0  # shape parameter, 3 < D < 6

    def psd(self, q: FloatArray) -> FloatArray:
        c = gamma(self.d / 2) / (np.pi**1.5 * abs(gamma((self.d - 3) / 2)))
        return np.asarray(
            self.sigma_n**2 * self.l_c**3 * c / (1 + (q * self.l_c) ** 2) ** (self.d / 2)
        )


@dataclass
class BornScattering:
    mus: float  # 1/mm
    g: float
    phase: Tabulated


def born(
    medium: WhittleMatern, n0: float, wavelength_um: float, n_mu: int = 20_001
) -> BornScattering:
    k = 2 * np.pi * n0 / (wavelength_um * 1e-3)
    mu = np.unique(np.cos(np.pi * np.linspace(0, 1, n_mu) ** 2)).astype(np.float64)
    q = 2 * k * np.sqrt(np.clip((1 - mu) / 2, 0, None))
    sigma = 2 * np.pi * k**4 / n0**2 * medium.psd(q)
    mus = float(2 * np.pi * np.trapezoid(sigma, mu))
    tab = Tabulated(mu, sigma)
    return BornScattering(mus, tab.g, tab)


def random_field(
    medium: WhittleMatern, shape: tuple[int, int, int], spacing: float, rng: np.random.Generator
) -> FloatArray:
    """δn on a periodic grid (nz, ny, nx) of isotropic spacing (mm)."""
    noise = rng.standard_normal(shape)
    qs = [2 * np.pi * np.fft.fftfreq(n, d=spacing) for n in shape]
    qz, qy, qx = np.meshgrid(*qs, indexing="ij")
    q = np.sqrt(qx**2 + qy**2 + qz**2)
    filt = np.sqrt(medium.psd(q) * (2 * np.pi) ** 3 / spacing**3)
    return np.asarray(np.real(np.fft.ifftn(np.fft.fftn(noise) * filt)))


def quantise(dn: FloatArray, n0: float, levels: int = 64) -> tuple[IntArray, list[complex]]:
    """Labels and indices for BWPM: δn in `levels` equal steps over its range (the step is far
    below σ_n for 64 levels; mean preserved)."""
    lo, hi = float(dn.min()), float(dn.max())
    step = (hi - lo) / (levels - 1) if hi > lo else 1.0
    lab = np.clip(np.rint((dn - lo) / step), 0, levels - 1).astype(np.int64)
    return lab, [complex(n0 + lo + j * step) for j in range(levels)]


@dataclass
class Calibration:
    depths: FloatArray  # mm
    coherent: FloatArray  # |⟨E⟩|² / |E₀|² per depth
    total: FloatArray  # transmitted power / incident (forward)
    mus_fit: float  # from the coherent decay
    g_forward: float  # ⟨cos θ⟩ of the scattered light at the shallowest depth > 0 (forward half)
    born: BornScattering
    g_forward_born: float
    phase_measured: Tabulated | None = None  # forward half from BWPM, backward half Born-shaped

    def scattering(self) -> tuple[float, Tabulated]:
        """(μs, phase function) as measured in BWPM, for the Monte Carlo below the hand-off (D9:
        both regions describe the same medium – the realisation, not its Born limit)."""
        assert self.phase_measured is not None
        return self.mus_fit, self.phase_measured


def calibrate_with_bwpm(
    medium: WhittleMatern,
    n0: float,
    wavelength_um: float,
    shape: tuple[int, int, int],
    spacing: float,
    rng: np.random.Generator,
    keep_every: int = 10,
    backend: str = "auto",
    fit_range: tuple[float, float] = (0.9, 0.05),
    n_bins: int = 90,
) -> Calibration:
    """Propagate a plane wave through one realisation in BWPM (forward FPWPM) and measure the
    extinction of the coherent field (μs) and the angular distribution of the scattered field
    (g, forward half). The coherent power is fitted where it lies between fit_range (fractions);
    a laterally finite box leaves a floor of about 1/(number of speckles)."""
    from bwpm import fast
    from bwpm.volume3d import Bidirectional3D, Grid3D

    nz, ny, nx = shape
    dn = random_field(medium, shape, spacing, rng)
    lab, idx = quantise(dn, n0)
    solver = Bidirectional3D(
        Grid3D(nx, ny, spacing, spacing, spacing, wavelength_um), idx, lab, complex(n0), None
    )
    e = np.zeros((ny, nx, 2), dtype=np.complex128)
    e[:, :, 0] = 1.0
    keep = tuple(range(0, nz + 1, keep_every))
    res = fast.run_3d(solver, e, iterations=0, keep=keep, backend=backend)
    coh, tot = [], []
    for plane in keep:
        f = res.planes[plane][0]
        coh.append(float(np.abs(np.mean(f[:, :, 0])) ** 2 + np.abs(np.mean(f[:, :, 1])) ** 2))
        tot.append(float(np.mean(np.sum(np.abs(f) ** 2, axis=2))))
    depths = np.array(keep, dtype=np.float64) * spacing
    c = np.array(coh)
    sel = (c < fit_range[0]) & (c > fit_range[1])
    sel[0] = False
    mus_fit = (
        float(-np.polyfit(depths[sel], np.log(c[sel]), 1)[0]) if sel.sum() >= 2 else float("nan")
    )
    # angular distribution at the first kept depth > 0
    f = res.planes[keep[1]][0]
    k = 2 * np.pi * n0 / (wavelength_um * 1e-3)
    kx = 2 * np.pi * np.fft.fftfreq(nx, spacing)
    ky = 2 * np.pi * np.fft.fftfreq(ny, spacing)
    kxx, kyy = np.meshgrid(kx, ky)
    kp2 = kxx**2 + kyy**2
    prop = kp2 < k**2
    ct = np.sqrt(np.clip(1 - kp2 / k**2, 0, 1))
    s = np.zeros((ny, nx))
    for comp in range(2):
        s += np.abs(np.fft.fft2(f[:, :, comp] - np.mean(f[:, :, comp]))) ** 2
    # power per solid angle ∝ |F|² cos²θ; a k⊥ cell covers dΩ = dk⊥²/(k² cos θ), so the weight of
    # the cell in an average over directions is |F|² cos θ (found in the S9b-2 server runs)
    w = np.where(prop, s * ct, 0.0)
    g_fwd = float(np.sum(w * ct) / np.sum(w))
    b = born(medium, n0, wavelength_um)
    fw = b.phase.mu > 0
    g_b = float(
        np.trapezoid(b.phase.mu[fw] * b.phase.p[fw], b.phase.mu[fw])
        / np.trapezoid(b.phase.p[fw], b.phase.mu[fw])
    )
    phase = _measured_phase(ct[prop], w[prop], b, n_bins)
    return Calibration(depths, c, np.array(tot), mus_fit, g_fwd, b, g_b, phase)


def _measured_phase(
    cos_t: FloatArray, weight: FloatArray, b: BornScattering, n_bins: int
) -> Tabulated:
    """Phase function from the scattered power per k⊥ cell (weight ∝ σ dΩ): histogram over θ in
    the forward half → density per sr; the backward half (not seen by the forward propagation)
    takes the Born shape, scaled to meet the measured density near 90°; normalised over the
    sphere."""
    theta = np.arccos(np.clip(cos_t, 0.0, 1.0))
    edges = np.linspace(0.0, np.pi / 2, n_bins + 1)
    h, _ = np.histogram(theta, bins=edges, weights=weight)
    d_omega = 2 * np.pi * (np.cos(edges[:-1]) - np.cos(edges[1:]))
    dens = h / d_omega
    mu_f = np.cos(0.5 * (edges[:-1] + edges[1:]))[::-1]  # ascending in μ
    p_f = dens[::-1]
    tail = slice(0, max(n_bins // 20, 3))  # bins next to 90°
    born_90 = b.phase.pdf(mu_f[tail])
    scale = float(np.sum(p_f[tail]) / np.sum(born_90)) if np.sum(born_90) > 0 else 0.0
    mu_b = np.linspace(-1.0, 0.0, n_bins // 2 + 1)[:-1].astype(np.float64)
    p_b = scale * b.phase.pdf(mu_b)
    mu = np.concatenate([mu_b, [0.0], mu_f, [1.0]])
    p = np.concatenate([p_b, [0.5 * (p_b[-1] + p_f[0])], p_f, [p_f[-1]]])
    order = np.argsort(mu)
    mu, p = mu[order], np.maximum(p[order], 0.0)
    keep = np.concatenate([[True], np.diff(mu) > 0])
    return Tabulated(mu[keep].astype(np.float64), p[keep].astype(np.float64))
