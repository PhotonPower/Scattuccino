"""Photon transport with a photon protocol and volume scoring (stage S2, docs/spezifikation.md §3).

Same physics as the S1 reference (`scattuccino.layers.run`) plus

- a **protocol** of every packet that leaves the medium: exit position (x, y), direction, weight,
  total and per-layer geometric path length (mm), number of scattering events, maximum depth;
- **volume scoring** of the absorbed power on a cylindrical grid (r, z), from which the fluence is
  Φ = A/(μa V);
- two absorption schemes:
  - "implicit" (as S1, MCML): steps with μt, absorbed fraction μa/μt at each interaction;
  - "mbl" (microscopic Beer–Lambert): steps with μs only, the weight is multiplied by
    exp(−Σ μa,i L_i) along the path. The paths then do not depend on μa, so a protocol computed once
    can be re-weighted exactly to any other absorption (`reweight`), the basis of fast fitting
    ("white Monte Carlo") and of the partial path lengths.

Units: mm, 1/mm; time of flight t = Σ n_i L_i / c with c in mm/ps.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from scattuccino._types import FloatArray, IntArray
from scattuccino.layers import Stack, _scatter, fresnel_unpolarised

C_MM_PER_PS = 0.299792458


@dataclass
class Protocol:
    """Packets that left the medium. side: −1 top (reflection), +1 bottom (transmission)."""

    side: np.ndarray  # type: ignore[type-arg]
    x: FloatArray
    y: FloatArray
    ux: FloatArray
    uy: FloatArray
    uz: FloatArray
    weight: FloatArray
    path: FloatArray  # (N, layers) geometric path length per layer, mm
    scatterings: np.ndarray  # type: ignore[type-arg]
    max_depth: FloatArray
    n_layer: FloatArray  # refractive index per layer (for times of flight)
    packets: int
    mua_used: FloatArray  # μa per layer applied in the weights ("mbl": inside `weight`)
    packet: IntArray | None = None  # packet number of each exit (voxel transport, S4/S5)
    momentum: FloatArray | None = None  # (N, layers) Σ (1 − cos θ) per layer (S8, dynamics)

    @property
    def total_path(self) -> FloatArray:
        return np.asarray(self.path.sum(axis=1))

    @property
    def time_ps(self) -> FloatArray:
        return np.asarray(self.path @ self.n_layer / C_MM_PER_PS)

    @property
    def rho(self) -> FloatArray:
        return np.asarray(np.hypot(self.x, self.y))


@dataclass
class Grid:
    """Cylindrical scoring grid: radial step dr (mm), nr bins, axial step dz, nz bins."""

    dr: float
    nr: int
    dz: float
    nz: int


@dataclass
class Sources:
    """Isotropic point sources inside the medium (S7: fluorescence, Raman): positions (mm),
    layer index and weight (power) per packet. Replaces the incident beam."""

    x: FloatArray
    y: FloatArray
    z: FloatArray
    layer: IntArray
    weight: FloatArray


@dataclass
class Simulation:
    stack: Stack
    r_specular: float
    protocol: Protocol
    absorbed: FloatArray  # per layer
    absorbed_grid: FloatArray | None  # (nr, nz), per incident power
    grid: Grid | None
    seed: int
    events: dict[str, FloatArray] | None = None  # absorption events (S7): x, y, z, layer, dw

    def fluence(self) -> FloatArray:
        """Φ(r, z) per incident power (1/mm²) from the absorbed grid: A/(μa · V)."""
        assert self.grid is not None and self.absorbed_grid is not None
        g = self.grid
        r_edges = np.arange(g.nr + 1) * g.dr
        vol = np.pi * (r_edges[1:] ** 2 - r_edges[:-1] ** 2)[:, None] * g.dz
        z_c = (np.arange(g.nz) + 0.5) * g.dz
        bounds = self.stack.boundaries()
        layer = np.clip(
            np.searchsorted(bounds, z_c, side="right") - 1, 0, len(self.stack.layers) - 1
        )
        mua = np.array([self.stack.layers[i].mua for i in layer])
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.asarray(
                np.where(mua[None, :] > 0, self.absorbed_grid / (mua[None, :] * vol), np.nan)
            )


def simulate(
    stack: Stack,
    packets: int,
    seed: int = 1,
    cos_incident: float = 1.0,
    absorption: str = "implicit",
    grid: Grid | None = None,
    max_path: float = np.inf,
    w_threshold: float = 1e-4,
    roulette_m: float = 10.0,
    batch: int = 100_000,
    sources: Sources | None = None,
    record_absorption: bool = False,
) -> Simulation:
    """Monte Carlo with protocol (see module docstring). max_path ends packets whose path exceeds it
    (their weight is lost; use only where its share is negligible).

    sources (S7): isotropic packets from points inside the medium instead of the incident beam
    (`packets` must equal their number; weights are taken as given, not divided). With
    record_absorption the absorption events (position, layer, absorbed weight per incident power)
    of the implicit scheme are returned in `Simulation.events`."""
    if absorption not in ("implicit", "mbl"):
        raise ValueError("absorption must be 'implicit' or 'mbl'")
    mbl = absorption == "mbl"
    rng = np.random.Generator(np.random.PCG64(seed))
    z_b = stack.boundaries()
    nl = len(stack.layers)
    n_layer = np.array([layer.n for layer in stack.layers])
    mua = np.array([layer.mua for layer in stack.layers])
    mus = np.array([layer.mus for layer in stack.layers])
    mut = mus if mbl else mua + mus
    phases = [layer.phase_function for layer in stack.layers]
    r_sp_arr, cos_t0 = fresnel_unpolarised(stack.n_above, n_layer[0], np.array([cos_incident]))
    r_sp = float(r_sp_arr[0])
    sin_t0 = stack.n_above / n_layer[0] * np.sqrt(1 - cos_incident**2)
    out: dict[str, list[FloatArray]] = {
        k: [] for k in ("side", "x", "y", "ux", "uy", "uz", "w", "path", "nsc", "zmax", "mom")
    }
    absorbed = np.zeros(nl)
    agrid = np.zeros((grid.nr, grid.nz)) if grid else None
    ev: dict[str, list[FloatArray]] = {k: [] for k in ("x", "y", "z", "layer", "dw")}
    if sources is not None:
        if len(sources.x) != packets:
            raise ValueError("packets must equal the number of sources")
        if mbl:
            raise ValueError("internal sources are implemented for the implicit scheme")
        r_sp = 0.0
    done = 0
    while done < packets:
        m = min(batch, packets - done)
        done += m
        x: FloatArray = np.zeros(m)
        y: FloatArray = np.zeros(m)
        z: FloatArray = np.zeros(m)
        ux: FloatArray = np.full(m, sin_t0)
        uy: FloatArray = np.zeros(m)
        uz: FloatArray = np.full(m, float(cos_t0[0]))
        w: FloatArray = np.full(m, 1.0 - r_sp)
        layer: IntArray = np.zeros(m, dtype=np.int64)
        if sources is not None:
            sl = slice(done - m, done)
            x, y, z = (
                sources.x[sl].astype(np.float64),
                sources.y[sl].astype(np.float64),
                sources.z[sl].astype(np.float64),
            )
            layer = sources.layer[sl].astype(np.int64)
            w = sources.weight[sl].astype(np.float64) * packets  # undone by the final division
            cz = 2 * rng.random(m) - 1
            ph = 2 * np.pi * rng.random(m)
            sz = np.sqrt(1 - cz * cz)
            ux, uy, uz = sz * np.cos(ph), sz * np.sin(ph), cz
        alive = np.ones(m, dtype=bool)
        path = np.zeros((m, nl))
        mom = np.zeros((m, nl))
        nsc = np.zeros(m, dtype=np.int64)
        zmax = np.zeros(m)
        s_left = np.zeros(m)

        def deposit(
            dw: FloatArray,
            zpos: FloatArray,
            lay: IntArray,
            xpos: FloatArray,
            ypos: FloatArray,
        ) -> None:
            np.add.at(absorbed, lay, dw)
            if agrid is not None and grid is not None:
                ir = np.minimum((np.hypot(xpos, ypos) / grid.dr).astype(np.int64), grid.nr - 1)
                iz = np.minimum((zpos / grid.dz).astype(np.int64), grid.nz - 1)
                ok = zpos >= 0
                np.add.at(agrid, (ir[ok], iz[ok]), dw[ok])

        def leave(
            ids: IntArray,
            side: int,
            w: FloatArray = w,
            path: FloatArray = path,
            x: FloatArray = x,
            y: FloatArray = y,
            ux: FloatArray = ux,
            uy: FloatArray = uy,
            uz: FloatArray = uz,
            nsc: IntArray = nsc,
            zmax: FloatArray = zmax,
            mom: FloatArray = mom,
        ) -> None:
            wout = w[ids]
            if mbl:
                wout = wout * np.exp(-(path[ids] @ mua))
            for key, val in (
                ("side", np.full(ids.size, side)),
                ("x", x[ids]),
                ("y", y[ids]),
                ("ux", ux[ids]),
                ("uy", uy[ids]),
                ("uz", uz[ids]),
                ("w", wout),
                ("path", path[ids]),
                ("nsc", nsc[ids]),
                ("zmax", zmax[ids]),
                ("mom", mom[ids]),
            ):
                out[key].append(np.asarray(val))

        while np.any(alive):
            idx = np.nonzero(alive)[0]
            new = s_left[idx] <= 0
            s_left[idx[new]] = -np.log(rng.random(int(new.sum())))
            li = layer[idx]
            mt = mut[li]
            step = np.where(mt > 0, s_left[idx] / np.where(mt > 0, mt, 1.0), np.inf)
            zt, zbot = z_b[li], z_b[li + 1]
            uzi = uz[idx]
            with np.errstate(divide="ignore", invalid="ignore"):
                d_b = np.where(
                    uzi < 0, (zt - z[idx]) / uzi, np.where(uzi > 0, (zbot - z[idx]) / uzi, np.inf)
                )
            hit = d_b <= step
            dist = np.where(hit, d_b, step)
            x0, y0, z0 = x[idx].copy(), y[idx].copy(), z[idx].copy()
            x[idx] += ux[idx] * dist
            y[idx] += uy[idx] * dist
            z[idx] += uzi * dist
            path[idx, li] += dist
            zmax[idx] = np.maximum(zmax[idx], z[idx])
            if mbl:  # continuous absorption along the segment, scored at its midpoint
                dw = (
                    w[idx]
                    * np.exp(-(path[idx] @ mua - mua[li] * dist))
                    * (1 - np.exp(-mua[li] * dist))
                )
                deposit(dw, 0.5 * (z0 + z[idx]), li, 0.5 * (x0 + x[idx]), 0.5 * (y0 + y[idx]))
            s_left[idx] = np.where(hit, s_left[idx] - d_b * mt, 0.0)
            hb = idx[hit]
            if hb.size:
                up = uz[hb] < 0
                lb = layer[hb]
                n1 = n_layer[lb]
                nb = np.where(
                    up,
                    np.where(lb == 0, stack.n_above, n_layer[np.maximum(lb - 1, 0)]),
                    np.where(lb == nl - 1, stack.n_below, n_layer[np.minimum(lb + 1, nl - 1)]),
                )
                r = np.empty(hb.size)
                cos_t = np.empty(hb.size)
                for k in np.unique(np.column_stack([n1, nb]), axis=0):
                    sel = (n1 == k[0]) & (nb == k[1])
                    r[sel], cos_t[sel] = fresnel_unpolarised(float(k[0]), float(k[1]), uz[hb[sel]])
                refl = rng.random(hb.size) < r
                uz[hb[refl]] = -uz[hb[refl]]
                tt = hb[~refl]
                if tt.size:
                    upt = up[~refl]
                    ratio = n1[~refl] / nb[~refl]
                    ux[tt] *= ratio
                    uy[tt] *= ratio
                    uz[tt] = np.where(upt, -cos_t[~refl], cos_t[~refl])
                    lt = layer[tt]
                    top = upt & (lt == 0)
                    bot = (~upt) & (lt == nl - 1)
                    leave(tt[top], -1)
                    leave(tt[bot], +1)
                    alive[tt[top | bot]] = False
                    inside = ~(top | bot)
                    layer[tt[inside]] += np.where(upt[inside], -1, 1)
            ib = idx[~hit]
            if ib.size:
                lb = layer[ib]
                if not mbl:
                    dw = w[ib] * np.where(
                        mut[lb] > 0, mua[lb] / np.where(mut[lb] > 0, mut[lb], 1.0), 0.0
                    )
                    deposit(dw, z[ib], lb, x[ib], y[ib])
                    if record_absorption:
                        for key, val in (
                            ("x", x[ib]),
                            ("y", y[ib]),
                            ("z", z[ib]),
                            ("layer", lb.astype(np.float64)),
                            ("dw", dw),
                        ):
                            ev[key].append(np.asarray(val, dtype=np.float64))
                    w[ib] -= dw
                xi_s = rng.random(ib.size)
                ct = np.empty(ib.size)
                for li_ in np.unique(lb):
                    sel_ = lb == li_
                    ct[sel_] = phases[int(li_)].sample(xi_s[sel_])
                phi = 2 * np.pi * rng.random(ib.size)
                ux[ib], uy[ib], uz[ib] = _scatter(ux[ib], uy[ib], uz[ib], ct, phi)
                nsc[ib] += 1
                mom[ib, lb] += 1.0 - ct
                w_eff = w[ib] * (np.exp(-(path[ib] @ mua)) if mbl else 1.0)
                low = ib[w_eff < w_threshold]
                if low.size:
                    survive = rng.random(low.size) < 1.0 / roulette_m
                    w[low[survive]] *= roulette_m
                    if mbl:  # the lost weight is absorbed in expectation; account it as absorbed
                        pass
                    w[low[~survive]] = 0.0
                    alive[low[~survive]] = False
            if np.isfinite(max_path):
                long = idx[alive[idx] & (path[idx].sum(axis=1) > max_path)]
                alive[long] = False
    cat = {
        k: (np.concatenate(v) if v else np.zeros((0, nl)) if k in ("path", "mom") else np.zeros(0))
        for k, v in out.items()
    }
    proto = Protocol(
        cat["side"].astype(np.int64),
        cat["x"],
        cat["y"],
        cat["ux"],
        cat["uy"],
        cat["uz"],
        cat["w"] / packets,
        cat["path"].reshape(-1, nl),
        cat["nsc"].astype(np.int64),
        cat["zmax"],
        n_layer,
        packets,
        mua.copy(),
        None,
        cat["mom"].reshape(-1, nl),
    )
    events = None
    if record_absorption:
        events = {k: (np.concatenate(v) if v else np.zeros(0)) for k, v in ev.items()}
        events["dw"] = events["dw"] / packets
    return Simulation(
        stack,
        r_sp,
        proto,
        absorbed / packets,
        agrid / packets if agrid is not None else None,
        grid,
        seed,
        events,
    )


def reweight(sim: Simulation, mua_new: FloatArray) -> FloatArray:
    """Weights of the protocol for other absorption coefficients per layer, exact for "mbl":
    w′ = w · exp(−Σ_i (μa,i′ − μa,i) L_i)."""
    p = sim.protocol
    delta = np.asarray(mua_new, dtype=np.float64) - p.mua_used
    return np.asarray(p.weight * np.exp(-(p.path @ delta)))
