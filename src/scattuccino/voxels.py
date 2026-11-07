"""Photon transport in voxel geometry (stage S4a).

A rectangular grid of voxels (nx, ny, nz) of size (dx, dy, dz) in mm, each voxel holding a
material label; materials are `scattuccino.layers.Layer` objects (thickness unused). The grid spans
x ∈ [−nx dx/2, nx dx/2], y likewise, z ∈ [0, nz dz]; light enters at (0, 0, 0) from z < 0
(index n_above), below z = nz dz lies n_below. Laterally the grid is either periodic (laterally
infinite media, comparable to layers) or open (packets leaving the side are recorded with
side = 0).

Algorithm as `scattuccino.transport` (implicit capture, Russian roulette) with counter-based
random numbers (one stream per packet, `scattuccino.random`) and a voxel walk:
each step ends at the next interaction or the next voxel face. At a face between different
refractive indices the packet is reflected or refracted (Fresnel, normal along the axis of the
face); between equal indices it passes, and the dimensionless remaining step s (in units of
1/μt) carries over, so voxels of equal material can be refined without changing a single random
number.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from scattuccino._types import FloatArray, IntArray
from scattuccino.layers import Layer, _scatter, fresnel_unpolarised
from scattuccino.random import CounterRandom
from scattuccino.transport import Protocol


@dataclass(frozen=True)
class VoxelMedium:
    labels: IntArray  # (nx, ny, nz)
    voxel: tuple[float, float, float]  # dx, dy, dz in mm
    materials: tuple[Layer, ...]
    n_above: float = 1.0
    n_below: float = 1.0
    lateral: str = "periodic"  # or "open"

    def __post_init__(self) -> None:
        lab = np.asarray(self.labels)
        if lab.ndim != 3 or lab.min() < 0 or lab.max() >= len(self.materials):
            raise ValueError("labels must be (nx, ny, nz) indices into materials")
        if self.lateral not in ("periodic", "open"):
            raise ValueError("lateral must be 'periodic' or 'open'")


@dataclass
class VoxelSimulation:
    medium: VoxelMedium
    r_specular: float
    protocol: Protocol  # path lengths per material
    absorbed: FloatArray  # per material
    absorbed_voxels: FloatArray  # (nx, ny, nz)
    seed: int


def simulate_voxels(
    medium: VoxelMedium,
    packets: int,
    seed: int = 1,
    cos_incident: float = 1.0,
    w_threshold: float = 1e-4,
    roulette_m: float = 10.0,
    batch: int = 100_000,
) -> VoxelSimulation:
    rng = CounterRandom(seed, packets)
    lab = np.asarray(medium.labels, dtype=np.int64)
    nx, ny, nz = lab.shape
    dx, dy, dz = medium.voxel
    size = np.array([dx, dy, dz])
    lo = np.array([-nx * dx / 2, -ny * dy / 2, 0.0])
    dims = np.array([nx, ny, nz])
    mats = medium.materials
    nm = len(mats)
    n_mat = np.array([m.n for m in mats])
    mua = np.array([m.mua for m in mats])
    mut = np.array([m.mua + m.mus for m in mats])
    phases = [m.phase_function for m in mats]
    periodic = medium.lateral == "periodic"
    # entry voxel: column through (0, 0)
    i0 = np.minimum((-lo[:2] / size[:2]).astype(np.int64), dims[:2] - 1)
    entry_mat = int(lab[i0[0], i0[1], 0])
    r_sp_arr, cos_t0 = fresnel_unpolarised(
        medium.n_above, n_mat[entry_mat], np.array([cos_incident])
    )
    r_sp = float(r_sp_arr[0])
    sin_t0 = medium.n_above / n_mat[entry_mat] * np.sqrt(1 - cos_incident**2)
    absorbed = np.zeros(nm)
    avox = np.zeros(lab.shape)
    out: dict[str, list[FloatArray]] = {
        k: [] for k in ("side", "x", "y", "ux", "uy", "uz", "w", "path", "nsc", "zmax")
    }
    done = 0
    while done < packets:
        m = min(batch, packets - done)
        gid = np.arange(done, done + m)  # global packet ids (streams)
        done += m
        pos = np.zeros((m, 3))
        u = np.zeros((m, 3))
        u[:, 0], u[:, 2] = sin_t0, float(cos_t0[0])
        idx = np.tile(np.array([i0[0], i0[1], 0]), (m, 1))
        w = np.full(m, 1.0 - r_sp)
        alive = np.ones(m, dtype=bool)
        path = np.zeros((m, nm))
        nsc = np.zeros(m, dtype=np.int64)
        zmax = np.zeros(m)
        s_left = np.zeros(m)
        while np.any(alive):
            a = np.nonzero(alive)[0]
            new = s_left[a] <= 0
            s_left[a[new]] = -np.log(rng.uniform(gid[a[new]]))
            mat = lab[idx[a, 0], idx[a, 1], idx[a, 2]]
            mt = mut[mat]
            step = np.where(mt > 0, s_left[a] / np.where(mt > 0, mt, 1.0), np.inf)
            # distance to the faces of the current voxel
            face_lo = lo + idx[a] * size
            ua = u[a]
            with np.errstate(divide="ignore", invalid="ignore"):
                d_ax = np.where(
                    ua > 0,
                    (face_lo + size - pos[a]) / ua,
                    np.where(ua < 0, (face_lo - pos[a]) / ua, np.inf),
                )
            d_ax = np.maximum(d_ax, 0.0)
            axis = np.argmin(d_ax, axis=1)
            d_face = d_ax[np.arange(a.size), axis]
            hit = d_face <= step
            dist = np.where(hit, d_face, step)
            pos[a] += ua * dist[:, None]
            path[a, mat] += dist
            zmax[a] = np.maximum(zmax[a], pos[a, 2])
            s_left[a] = np.where(hit, s_left[a] - d_face * mt, 0.0)
            # ---- faces
            hb, ax = a[hit], axis[hit]
            if hb.size:
                direction = np.sign(u[hb, ax]).astype(np.int64)
                nb_idx = idx[hb].copy()
                nb_idx[np.arange(hb.size), ax] += direction
                n1 = n_mat[lab[idx[hb, 0], idx[hb, 1], idx[hb, 2]]]
                out_z_top = (ax == 2) & (nb_idx[:, 2] < 0)
                out_z_bot = (ax == 2) & (nb_idx[:, 2] >= nz)
                out_lat = (ax < 2) & (
                    (nb_idx[np.arange(hb.size), ax] < 0)
                    | (nb_idx[np.arange(hb.size), ax] >= dims[ax])
                )
                if periodic and np.any(out_lat):
                    k = np.nonzero(out_lat)[0]
                    axk = ax[k]
                    nb_idx[k, axk] %= dims[axk]
                    pos[hb[k], axk] -= direction[k] * dims[axk] * size[axk]
                    out_lat[:] = False
                inside = ~(out_z_top | out_z_bot | out_lat)
                n2 = np.empty(hb.size)
                n2[out_z_top] = medium.n_above
                n2[out_z_bot] = medium.n_below
                n2[out_lat] = n1[out_lat]
                n2[inside] = n_mat[lab[nb_idx[inside, 0], nb_idx[inside, 1], nb_idx[inside, 2]]]
                cos_i = np.abs(u[hb, ax])
                r = np.zeros(hb.size)
                cos_t = cos_i.copy()
                diff = n1 != n2
                for pair in np.unique(np.column_stack([n1[diff], n2[diff]]), axis=0):
                    sel = diff & (n1 == pair[0]) & (n2 == pair[1])
                    r[sel], cos_t[sel] = fresnel_unpolarised(
                        float(pair[0]), float(pair[1]), cos_i[sel]
                    )
                # random numbers only where the index changes: equal-index faces (and periodic
                # wraps) consume none, so refining voxels of one material keeps the random stream
                refl = np.zeros(hb.size, dtype=bool)
                refl[diff] = rng.uniform(gid[hb[diff]]) < r[diff]
                # reflection: flip the normal component, stay in the voxel
                rr = np.nonzero(refl)[0]
                u[hb[rr], ax[rr]] *= -1
                # transmission: refract, move to the neighbour or leave
                tt = np.nonzero(~refl)[0]
                if tt.size:
                    ratio = n1[tt] / n2[tt]
                    ut = u[hb[tt]] * ratio[:, None]
                    ut[np.arange(tt.size), ax[tt]] = direction[tt] * cos_t[tt]
                    u[hb[tt]] = ut
                    leaving = (out_z_top | out_z_bot | out_lat)[tt]
                    ids = hb[tt[leaving]]
                    if ids.size:
                        side = np.where(
                            out_z_top[tt[leaving]], -1, np.where(out_z_bot[tt[leaving]], 1, 0)
                        )
                        for key, val in (
                            ("side", side.astype(np.float64)),
                            ("x", pos[ids, 0]),
                            ("y", pos[ids, 1]),
                            ("ux", u[ids, 0]),
                            ("uy", u[ids, 1]),
                            ("uz", u[ids, 2]),
                            ("w", w[ids]),
                            ("path", path[ids]),
                            ("nsc", nsc[ids].astype(np.float64)),
                            ("zmax", zmax[ids]),
                        ):
                            out[key].append(np.asarray(val))
                        alive[ids] = False
                    stay = tt[~leaving]
                    idx[hb[stay]] = nb_idx[stay]
            # ---- interactions
            ib = a[~hit]
            if ib.size:
                mb = lab[idx[ib, 0], idx[ib, 1], idx[ib, 2]]
                dw = w[ib] * np.where(
                    mut[mb] > 0, mua[mb] / np.where(mut[mb] > 0, mut[mb], 1.0), 0.0
                )
                np.add.at(absorbed, mb, dw)
                np.add.at(avox, (idx[ib, 0], idx[ib, 1], idx[ib, 2]), dw)
                w[ib] -= dw
                xi = rng.uniform(gid[ib])
                ct = np.empty(ib.size)
                for mm in np.unique(mb):
                    sel = mb == mm
                    ct[sel] = phases[int(mm)].sample(xi[sel])
                phi = 2 * np.pi * rng.uniform(gid[ib])
                u[ib, 0], u[ib, 1], u[ib, 2] = _scatter(u[ib, 0], u[ib, 1], u[ib, 2], ct, phi)
                nsc[ib] += 1
                low = ib[w[ib] < w_threshold]
                if low.size:
                    survive = rng.uniform(gid[low]) < 1.0 / roulette_m
                    w[low[survive]] *= roulette_m
                    w[low[~survive]] = 0.0
                    alive[low[~survive]] = False
    cat = {
        k: (np.concatenate(v) if v else (np.zeros((0, nm)) if k == "path" else np.zeros(0)))
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
        cat["path"].reshape(-1, nm),
        cat["nsc"].astype(np.int64),
        cat["zmax"],
        n_mat,
        packets,
        mua.copy(),
    )
    return VoxelSimulation(medium, r_sp, proto, absorbed / packets, avox / packets, seed)
