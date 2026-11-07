"""Layered media and the Monte-Carlo reference of stage S1 (docs/spezifikation.md §2).

Units: mm, 1/mm. Tissue occupies 0 < z < Σ d_i below an ambient medium (z < 0, index n_above)
and above another (z > Σ d_i, index n_below); light enters from z < 0. The algorithm follows
MCML (Wang, Jacques, Zheng 1995): weighted packets, implicit capture, Fresnel decisions at
interfaces, Russian roulette. Vectorised over packets with NumPy.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from scattuccino._types import FloatArray
from scattuccino.phase import HenyeyGreenstein, PhaseFunction


@dataclass(frozen=True)
class Layer:
    """Homogeneous layer: thickness d (mm, may be np.inf for the last layer), μa, μs (1/mm),
    Henyey-Greenstein anisotropy g, refractive index n and optionally another phase function
    (`scattuccino.phase`; S3, used by `scattuccino.transport` and `adding_doubling`)."""

    d: float
    mua: float
    mus: float
    g: float = 0.0
    n: float = 1.0
    phase: PhaseFunction | None = None

    @property
    def phase_function(self) -> PhaseFunction:
        return self.phase if self.phase is not None else HenyeyGreenstein(self.g)

    @property
    def mut(self) -> float:
        return self.mua + self.mus


@dataclass(frozen=True)
class Stack:
    layers: tuple[Layer, ...]
    n_above: float = 1.0
    n_below: float = 1.0

    def boundaries(self) -> FloatArray:
        """z of the interfaces: 0, d1, d1 + d2, …"""
        return np.concatenate([[0.0], np.cumsum([layer.d for layer in self.layers])])


@dataclass
class Result:
    """Per incident power 1; standard errors from the per-packet contributions."""

    r_specular: float
    r_diffuse: float
    t_diffuse: float
    absorbed: FloatArray  # per layer
    r_diffuse_err: float
    t_diffuse_err: float
    absorbed_err: FloatArray
    packets: int
    seed: int
    extra: dict[str, object] = field(default_factory=dict)

    @property
    def total(self) -> float:
        return self.r_specular + self.r_diffuse + self.t_diffuse + float(np.sum(self.absorbed))


def fresnel_unpolarised(n1: float, n2: float, cos_i: FloatArray) -> tuple[FloatArray, FloatArray]:
    """Unpolarised power reflectance R and cos θ_t for incidence from n1 into n2 (TIR: R = 1)."""
    cos_i = np.abs(cos_i)
    if n1 == n2:
        return np.zeros_like(cos_i), cos_i
    sin_t2 = (n1 / n2) ** 2 * (1 - cos_i**2)
    tir = sin_t2 >= 1.0
    cos_t = np.sqrt(np.clip(1 - sin_t2, 0.0, 1.0))
    rs = (n1 * cos_i - n2 * cos_t) / (n1 * cos_i + n2 * cos_t)
    rp = (n2 * cos_i - n1 * cos_t) / (n2 * cos_i + n1 * cos_t)
    r = np.where(tir, 1.0, 0.5 * (rs**2 + rp**2))
    return r, cos_t


def sample_hg(g: FloatArray, xi: FloatArray) -> FloatArray:
    """cos θ from Henyey-Greenstein (isotropic where g = 0)."""
    g = np.asarray(g, dtype=np.float64)
    iso = np.abs(g) < 1e-12
    gs = np.where(iso, 0.5, g)  # placeholder to avoid 0/0, replaced below
    t = (1 - gs * gs) / (1 - gs + 2 * gs * xi)
    cos_t = (1 + gs * gs - t * t) / (2 * gs)
    return np.clip(np.where(iso, 2 * xi - 1, cos_t), -1.0, 1.0)


def _scatter(
    ux: FloatArray, uy: FloatArray, uz: FloatArray, cos_t: FloatArray, phi: FloatArray
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """New direction after deflection θ, azimuth φ (MCML eq. 3.24-3.26 form)."""
    sin_t = np.sqrt(np.clip(1 - cos_t**2, 0.0, 1.0))
    cp, sp = np.cos(phi), np.sin(phi)
    near_pole = np.abs(uz) > 1 - 1e-12
    den = np.sqrt(np.clip(1 - uz * uz, 1e-300, None))
    nx = sin_t * (ux * uz * cp - uy * sp) / den + ux * cos_t
    ny = sin_t * (uy * uz * cp + ux * sp) / den + uy * cos_t
    nz = -sin_t * cp * den + uz * cos_t
    px = sin_t * cp
    py = sin_t * sp
    pz = np.sign(uz) * cos_t
    return (np.where(near_pole, px, nx), np.where(near_pole, py, ny), np.where(near_pole, pz, nz))


def run(
    stack: Stack,
    packets: int,
    seed: int = 1,
    cos_incident: float = 1.0,
    w_threshold: float = 1e-4,
    roulette_m: float = 10.0,
    batch: int = 200_000,
) -> Result:
    """Monte Carlo for a layered medium, collimated incidence from above (cos θ_i = cos_incident,
    azimuth 0). Returns reflectance, transmittance and absorption per layer with standard errors."""
    rng = np.random.Generator(np.random.PCG64(seed))
    z_b = stack.boundaries()
    nl = len(stack.layers)
    n_layer = np.array([layer.n for layer in stack.layers])
    mua = np.array([layer.mua for layer in stack.layers])
    mut = np.array([layer.mut for layer in stack.layers])
    g = np.array([layer.g for layer in stack.layers])
    n_out_top, n_out_bottom = stack.n_above, stack.n_below
    r_sp_arr, cos_t0 = fresnel_unpolarised(n_out_top, n_layer[0], np.array([cos_incident]))
    r_sp = float(r_sp_arr[0])
    sin_i = np.sqrt(1 - cos_incident**2)
    sin_t0 = n_out_top / n_layer[0] * sin_i
    sums = {"r": 0.0, "t": 0.0, "r2": 0.0, "t2": 0.0}
    a_sum = np.zeros(nl)
    a_sq = np.zeros(nl)
    done = 0
    while done < packets:
        m = min(batch, packets - done)
        done += m
        x = np.zeros(m)
        y = np.zeros(m)
        z = np.zeros(m)
        ux = np.full(m, sin_t0)
        uy = np.zeros(m)
        uz = np.full(m, float(cos_t0[0]))
        w = np.full(m, 1.0 - r_sp)
        layer = np.zeros(m, dtype=np.int64)
        alive = np.ones(m, dtype=bool)
        rd = np.zeros(m)
        td = np.zeros(m)
        ab = np.zeros((m, nl))
        s_left = np.zeros(m)  # dimensionless remaining step (in units of 1/μt), 0 = sample new
        while np.any(alive):
            idx = np.nonzero(alive)[0]
            new = s_left[idx] <= 0
            s_left[idx[new]] = -np.log(rng.random(int(new.sum())))
            li = layer[idx]
            mt = mut[li]
            step = np.where(mt > 0, s_left[idx] / np.where(mt > 0, mt, 1.0), np.inf)
            # distance to the layer boundary along uz
            zt, zbot = z_b[li], z_b[li + 1]
            uzi = uz[idx]
            with np.errstate(divide="ignore", invalid="ignore"):
                d_b = np.where(
                    uzi < 0, (zt - z[idx]) / uzi, np.where(uzi > 0, (zbot - z[idx]) / uzi, np.inf)
                )
            hit = d_b <= step
            # move
            dist = np.where(hit, d_b, step)
            x[idx] += ux[idx] * dist
            y[idx] += uy[idx] * dist
            z[idx] += uzi * dist
            s_left[idx] = np.where(hit, s_left[idx] - d_b * mt, 0.0)
            # --- boundary crossings
            hb = idx[hit]
            if hb.size:
                up = uz[hb] < 0
                lb = layer[hb]
                n1 = n_layer[lb]
                nb = np.where(
                    up,
                    np.where(lb == 0, n_out_top, n_layer[np.maximum(lb - 1, 0)]),
                    np.where(lb == nl - 1, n_out_bottom, n_layer[np.minimum(lb + 1, nl - 1)]),
                )
                r = np.empty(hb.size)
                cos_t = np.empty(hb.size)
                for k in np.unique(np.column_stack([n1, nb]), axis=0):
                    sel = (n1 == k[0]) & (nb == k[1])
                    r[sel], cos_t[sel] = fresnel_unpolarised(float(k[0]), float(k[1]), uz[hb[sel]])
                refl = rng.random(hb.size) < r
                # reflect
                rr = hb[refl]
                uz[rr] = -uz[rr]
                # transmit (refract)
                tt = hb[~refl]
                if tt.size:
                    upt = up[~refl]
                    ratio = n1[~refl] / nb[~refl]
                    ux[tt] *= ratio
                    uy[tt] *= ratio
                    uz[tt] = np.where(upt, -cos_t[~refl], cos_t[~refl])
                    lt = layer[tt]
                    leave_top = upt & (lt == 0)
                    leave_bot = (~upt) & (lt == nl - 1)
                    rd[tt[leave_top]] += w[tt[leave_top]]
                    td[tt[leave_bot]] += w[tt[leave_bot]]
                    alive[tt[leave_top | leave_bot]] = False
                    inside = tt[~(leave_top | leave_bot)]
                    layer[inside] += np.where(upt[~(leave_top | leave_bot)], -1, 1)
                    # remaining step continues with the new layer's μt (dimensionless s_left kept)
            # --- interactions
            ib = idx[~hit]
            if ib.size:
                lb = layer[ib]
                dw = w[ib] * np.where(
                    mut[lb] > 0, mua[lb] / np.where(mut[lb] > 0, mut[lb], 1.0), 0.0
                )
                np.add.at(ab, (ib, lb), dw)
                w[ib] -= dw
                ct = sample_hg(g[lb], rng.random(ib.size))
                phi = 2 * np.pi * rng.random(ib.size)
                ux[ib], uy[ib], uz[ib] = _scatter(ux[ib], uy[ib], uz[ib], ct, phi)
                low = ib[w[ib] < w_threshold]
                if low.size:
                    survive = rng.random(low.size) < 1.0 / roulette_m
                    w[low[survive]] *= roulette_m
                    w[low[~survive]] = 0.0
                    alive[low[~survive]] = False
        sums["r"] += rd.sum()
        sums["t"] += td.sum()
        sums["r2"] += (rd**2).sum()
        sums["t2"] += (td**2).sum()
        a_sum += ab.sum(axis=0)
        a_sq += (ab**2).sum(axis=0)

    def err(s: float, s2: float) -> float:
        mean = s / packets
        return float(np.sqrt(max(s2 / packets - mean**2, 0.0) / packets))

    return Result(
        r_sp,
        sums["r"] / packets,
        sums["t"] / packets,
        a_sum / packets,
        err(sums["r"], sums["r2"]),
        err(sums["t"], sums["t2"]),
        np.sqrt(np.maximum(a_sq / packets - (a_sum / packets) ** 2, 0) / packets),
        packets,
        seed,
    )
