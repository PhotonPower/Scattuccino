"""Back scattering of a random index slab: BWPM against Born and Monte Carlo (stage S9c-4b).

A plane wave (x polarised, λ = 0.8 µm, host n₀ = 1.38 everywhere) falls on a slab of thickness L
holding one realisation of a Whittle-Matérn δn (l_c = 0.5 µm, D = 4) given to BWPM as the weak
fluctuation (B7/B8, `model="increment"`, iterations = 1). Per realisation the reflected field is
binned into dR/dΩ over |cos θ| (flux per plane wave with E_z, `plane_wave_power`, per solid angle
dΩ = dk⊥²/(k² cos θ)); the same is done for the first-order Born field of the same realisation
(`microstructure.born_reflection`). The Monte-Carlo run bins the reflected packets of the same
slab with the vector Born phase function.

    python tools/backscatter_bwpm.py run --set A --h 0.05 --seed 1 --backend cuda
    python tools/backscatter_bwpm.py mc --set B --packets 20000000 --seed 1
    python tools/backscatter_bwpm.py table

Sets: A single scattering (L = 2 µm), B multiple scattering (L = 40 µm, h = 0.05 µm only);
σ_n = 0.0143 (vector Born μs ≈ 11.9/mm, stable increment). Every run appends one JSON line to
tools/out/backscatter_runs.jsonl. Units: mm in the code, µm in the output.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

import scattuccino.microstructure as ms
from scattuccino.layers import Layer, Stack
from scattuccino.transport import simulate

OUT = Path(__file__).parent / "out"
RUNS = OUT / "backscatter_runs.jsonl"
N0, LAM, SIGMA, L_C = 1.38, 0.8, 0.0143, 0.5e-3
SETS = {"A": 0.002, "B": 0.04}  # slab thickness, mm
BOX = 0.0128  # lateral box, mm (n = BOX/h)
EDGES = np.linspace(0.3, 0.9, 7)  # |cos θ| bins, away from grazing and from the slab-edge peak


def medium() -> ms.WhittleMatern:
    return ms.WhittleMatern(SIGMA, L_C, 4.0)


def binned(field: np.ndarray, h: float) -> np.ndarray:  # type: ignore[type-arg]
    """dR/dΩ per |cos θ| bin of a reflected field (ny, nx, 2) in real space (incident flux 1)."""
    n = field.shape[0]
    k = 2 * np.pi * N0 / (LAM * 1e-3)
    kx = 2 * np.pi * np.fft.fftfreq(n, h)
    kxx, kyy = np.meshgrid(kx, kx)
    ct = np.sqrt(np.clip(1 - (kxx**2 + kyy**2) / k**2, 0, 1))
    use = (ct >= EDGES[0]) & (ct < EDGES[-1])
    flux = ms.plane_wave_power(field, kxx, kyy, k) / n**4 * ct
    d_omega = (kx[1] - kx[0]) ** 2 / (k**2 * np.where(use, ct, 1.0))
    i = np.digitize(ct[use], EDGES) - 1
    return np.asarray(
        np.bincount(i, weights=flux[use], minlength=6)
        / np.bincount(i, weights=d_omega[use], minlength=6)
    )


def born_continuous(thickness: float) -> np.ndarray:  # type: ignore[type-arg]
    """Single scattering with the continuous vector Born phase function, attenuated:
    dR/dΩ = μs p(−c) (1 − e^(−x))/(x/L), x = μs L (1 + 1/c), at the bin centres c."""
    b = ms.born(medium(), N0, LAM)
    c = 0.5 * (EDGES[1:] + EDGES[:-1])
    x = b.mus * thickness * (1 + 1 / c)
    return np.asarray(b.mus * b.phase.pdf(-c) * (1 - np.exp(-x)) / (x / thickness))


def run_bwpm(name: str, h: float, seed: int, backend: str) -> dict[str, object]:
    from bwpm import fast
    from bwpm.volume3d import Bidirectional3D, Grid3D

    thickness = SETS[name]
    n = round(BOX / h)
    t0 = time.time()
    dn = ms.random_field(medium(), (round(thickness / h), n, n), h, np.random.default_rng(seed))
    lab, idx, fl = ms.bwpm_medium(dn, N0, "increment")
    solver = Bidirectional3D(
        Grid3D(n, n, h, h, h, LAM), idx, lab, complex(N0), None, fluctuation=fl
    )
    e = np.zeros((n, n, 2), dtype=np.complex128)
    e[:, :, 0] = 1.0
    res = fast.run_3d(solver, e, iterations=1, backend=backend)
    t_bwpm = time.time() - t0
    born_field = np.fft.ifft2(ms.born_reflection(dn, N0, LAM, h), axes=(0, 1))
    return {
        "kind": "bwpm",
        "set": name,
        "h_um": h * 1e3,
        "seed": seed,
        "n": n,
        "L_um": thickness * 1e3,
        "bwpm": binned(np.asarray(res.reflected), h).tolist(),
        "born_discrete": binned(born_field, h).tolist(),
        "transmitted_power": float(np.mean(np.sum(np.abs(res.transmitted) ** 2, axis=-1))),
        "time_s": t_bwpm,
    }


def run_mc(name: str, packets: int, seed: int) -> dict[str, object]:
    thickness = SETS[name]
    b = ms.born(medium(), N0, LAM)
    stack = Stack((Layer(thickness, 0.0, b.mus, 0.0, N0, phase=b.phase),), N0, N0)
    t0 = time.time()
    p = simulate(stack, packets, seed=seed).protocol
    top = p.side == -1
    c = -p.uz[top]
    hist, _ = np.histogram(c, bins=EDGES, weights=p.weight[top])
    hist2, _ = np.histogram(c, bins=EDGES, weights=p.weight[top] ** 2)
    d_omega = 2 * np.pi * np.diff(EDGES)
    return {
        "kind": "mc",
        "set": name,
        "seed": seed,
        "packets": packets,
        "L_um": thickness * 1e3,
        "mc": (hist / d_omega).tolist(),
        "mc_se": (np.sqrt(hist2) / d_omega).tolist(),
        "time_s": time.time() - t0,
    }


def table() -> None:
    rows = [json.loads(line) for line in RUNS.read_text().splitlines() if line.strip()]
    c = 0.5 * (EDGES[1:] + EDGES[:-1])
    print("|cos θ| bins:", " ".join(f"{x:.2f}" for x in c))
    for name in SETS:
        ref = born_continuous(SETS[name])
        mc = [r for r in rows if r["kind"] == "mc" and r["set"] == name]
        if mc:
            m = np.mean([r["mc"] for r in mc], axis=0)
            print(f"set {name} MC / continuous Born:", " ".join(f"{x:.3f}" for x in m / ref))
        for h in sorted({r["h_um"] for r in rows if r["kind"] == "bwpm" and r["set"] == name}):
            sel = [r for r in rows if r["kind"] == "bwpm" and r["set"] == name and r["h_um"] == h]
            bw = np.array([r["bwpm"] for r in sel])
            bd = np.array([r["born_discrete"] for r in sel])
            se = bw.std(axis=0) / np.sqrt(len(sel)) / ref
            print(f"set {name} h {h:.3f} µm, {len(sel)} realisations")
            print("  BWPM / continuous Born:   ", " ".join(f"{x:.3f}" for x in bw.mean(0) / ref))
            print("  ± (standard error)        ", " ".join(f"{x:.3f}" for x in se))
            print("  discrete Born / continuous:", " ".join(f"{x:.3f}" for x in bd.mean(0) / ref))
            print(
                "  BWPM / discrete Born:     ",
                " ".join(f"{x:.3f}" for x in bw.mean(0) / bd.mean(0)),
            )
            if mc:
                print("  BWPM / MC:                ", " ".join(f"{x:.3f}" for x in bw.mean(0) / m))


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--set", choices=tuple(SETS), required=True)
    r.add_argument("--h", type=float, required=True, help="grid step, µm")
    r.add_argument("--seed", type=int, required=True)
    r.add_argument("--backend", default="auto")
    m = sub.add_parser("mc")
    m.add_argument("--set", choices=tuple(SETS), required=True)
    m.add_argument("--packets", type=int, required=True)
    m.add_argument("--seed", type=int, required=True)
    sub.add_parser("table")
    a = ap.parse_args()
    if a.cmd == "table":
        table()
        return
    if a.cmd == "run":
        row = run_bwpm(a.set, a.h * 1e-3, a.seed, a.backend)
    else:
        row = run_mc(a.set, a.packets, a.seed)
    OUT.mkdir(exist_ok=True)
    with RUNS.open("a") as fh:
        fh.write(json.dumps(row) + "\n")
    print({k: v for k, v in row.items() if not isinstance(v, list)})


if __name__ == "__main__":
    main()
