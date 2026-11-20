"""Acceptance runs of the hybrid OCT with the coherent hand-off, tight focus (issue 34, S9b-6).

    python tools/abnahme_s9b6.py cal                                      # step 0 (set B only)
    python tools/abnahme_s9b6.py run --set A --zh 0.005 --seed 1
    python tools/abnahme_s9b6.py run --set A --zh 0.005 --seed 1 --packets 1000000 --pilot
    python tools/abnahme_s9b6.py eval                                     # tables and verdict

Units: mm for lengths, µm for the wavelength. Outputs go to tools/out/abnahme_s9b6/.
Set A: vector Born (default), σ_H = 0.4 µm. Set C: as A with σ_H = 0.8 µm (z_h > 0).
Set B: scattering measured in BWPM (from `cal`), σ_H = 0.4 µm.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from scattuccino.hybrid import Beam, hybrid_oct
from scattuccino.microstructure import WhittleMatern, born, calibrate_with_bwpm
from scattuccino.phase import Tabulated

OUT = Path(__file__).parent / "out" / "abnahme_s9b6"
MED = WhittleMatern(0.0572, 0.5e-3)
N0, WAVELENGTH_UM = 1.38, 0.8
Z = np.linspace(0.0, 0.06, 31)
ZH_A = (0.0, 0.0025, 0.005, 0.010)
SIGMA = {"A": 0.4e-3, "B": 0.4e-3, "C": 0.8e-3}
SETS = {
    "A": (ZH_A, (1, 2, 3, 4, 5)),
    "C": ((0.005, 0.010), (1, 2, 3, 4, 5)),
    "B": ((0.0, 0.005), (1, 2, 3)),
}
Z_MIN = 0.015
WIN = (0.016 - 1e-9, 0.022 + 1e-9)  # window of the near-field mean ratio (mm)


def calibrate() -> None:
    cal = calibrate_with_bwpm(
        MED, N0, WAVELENGTH_UM, (134, 256, 256), 1.5e-4, np.random.default_rng(11),
        keep_every=10, backend="cuda",
    )  # fmt: skip
    mus, phase = cal.scattering()
    OUT.mkdir(parents=True, exist_ok=True)
    row = {
        "mus_fit": cal.mus_fit, "mus_born": cal.born.mus, "g_forward": cal.g_forward,
        "g_forward_born": cal.g_forward_born, "g_measured_phase": phase.g,
        "max_power_dev": float(np.max(np.abs(cal.total - 1.0))),
        "coherent": cal.coherent.tolist(), "depths": cal.depths.tolist(),
    }  # fmt: skip
    (OUT / "kalibrierung.json").write_text(json.dumps(row, indent=1))
    np.savez(OUT / "kalibrierung_phase.npz", mus=mus, mu=phase.mu, p=phase.p)
    print({k: v for k, v in row.items() if not isinstance(v, list)})


def load_scattering() -> tuple[float, Tabulated]:
    d = np.load(OUT / "kalibrierung_phase.npz")
    return float(d["mus"]), Tabulated(d["mu"], d["p"])


def run(satz: str, z_h: float, seed: int, packets: int, pilot: bool) -> None:
    sc = load_scattering() if satz == "B" else None
    t0 = time.time()
    res = hybrid_oct(
        MED, N0, 1.0, WAVELENGTH_UM, Beam(2e-3, 0.012), z_h, Z, 0.004, packets, seed=seed,
        nx=256, dx=1.5e-4, max_path=0.15, sigma_husimi=SIGMA[satz], backend="cuda",
        scattering=sc, coherent="wigner",
    )  # fmt: skip
    dt = time.time() - t0
    name = f"lauf_{satz}_{z_h:g}_{seed}" + ("_pilot" if pilot else "")
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez(
        OUT / f"{name}.npz", z=res.z, signal=res.signal, power_handoff=res.power_handoff,
        power_coherent=res.power_coherent, runtime_s=dt, packets=packets,
    )  # fmt: skip
    print(
        f"{name}: {dt:.0f} s, power_handoff {res.power_handoff:.4f}, "
        f"power_coherent {res.power_coherent:.4f}"
    )


def load(satz: str, z_h: float, key: str = "signal") -> np.ndarray:
    seeds = SETS[satz][1]
    return np.array([np.load(OUT / f"lauf_{satz}_{z_h:g}_{s}.npz")[key] for s in seeds])


def stats(a: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return a.mean(0), a.std(0, ddof=1) / np.sqrt(len(a))


def ratio(
    m1: np.ndarray, e1: np.ndarray, m0: np.ndarray, e0: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    r = m1 / m0
    return r, np.abs(r) * np.hypot(e1 / m1, e0 / m0)


def window_mean(r: np.ndarray, e: np.ndarray) -> tuple[float, float]:
    w = (Z >= WIN[0]) & (Z <= WIN[1])
    return float(r[w].mean()), float(np.sqrt(np.sum(e[w] ** 2)) / w.sum())


def evaluate() -> None:
    sel = Z >= Z_MIN - 1e-12
    mean, sem = {}, {}
    for zh in ZH_A:
        mean[zh], sem[zh] = stats(load("A", zh))
    head = " | ".join(f"z_h={zh * 1e3:g} µm" for zh in ZH_A)
    print(f"Signal (mean ± standard error, 5 seeds, set A)\n\n| z µm | {head} |")
    print("|" + " --- |" * (len(ZH_A) + 1))
    for i in np.where(sel)[0]:
        line = " | ".join(f"{mean[zh][i]:.4g} ± {sem[zh][i]:.2g}" for zh in ZH_A)
        print(f"| {Z[i] * 1e3:.0f} | {line} |")
    print("\nRatio to z_h = 0 (✗: criterion missed)\n")
    print("| z µm | " + " | ".join(f"z_h={zh * 1e3:g} µm" for zh in ZH_A[1:]) + " |")
    print("|" + " --- |" * len(ZH_A))
    rat = {zh: ratio(mean[zh], sem[zh], mean[0.0], sem[0.0]) for zh in ZH_A[1:]}
    bad = []
    for i in np.where(sel)[0]:
        cells = []
        for zh in ZH_A[1:]:
            r, e = rat[zh][0][i], rat[zh][1][i]
            ok = abs(r - 1) <= 3 * e or abs(r - 1) <= 0.15
            if not ok:
                bad.append((round(Z[i] * 1e3), zh * 1e3, round(float(r), 3), round(float(e), 3)))
            cells.append(f"{r:.3f} ± {e:.3f}{'' if ok else ' ✗'}")
        print(f"| {Z[i] * 1e3:.0f} | " + " | ".join(cells) + " |")
    print(f"\nDeviating (z µm, z_h µm, ratio, error): {bad}")
    print("PASSED" if not bad else "NOT PASSED")
    print("\n1. Mean ratio over 16-22 µm:")
    for zh in ZH_A[1:]:
        m, e = window_mean(*rat[zh])
        print(f"   z_h = {zh * 1e3:g} µm: {m:.3f} ± {e:.3f}")
    mus = born(MED, N0, WAVELENGTH_UM).mus
    print(f"\n2. power_coherent against exp(-mu_s z_h), mu_s = {mus:.2f} /mm:")
    for zh in ZH_A:
        pc = load("A", zh, "power_coherent")
        m, e = pc.mean(), pc.std(ddof=1) / np.sqrt(len(pc))
        ex = np.exp(-mus * zh)
        print(f"   z_h = {zh * 1e3:g} µm: {m:.4f} ± {e:.4f}, exp {ex:.4f}, ratio {m / ex:.3f}")
    print("\n3. Set C against set A (ratio C/A):")
    for zh in SETS["C"][0]:
        mc, ec = stats(load("C", zh))
        r, e = ratio(mc, ec, mean[zh], sem[zh])
        print(f"   z_h = {zh * 1e3:g} µm, window mean {window_mean(r, e)}, per depth:")
        print("   " + " ".join(f"{x:.3f}" for x in r[sel]))
    print("\n4. Set B against set A (ratio B/A):")
    for zh in SETS["B"][0]:
        mb, eb = stats(load("B", zh))
        r, e = ratio(mb, eb, mean[zh], sem[zh])
        m = float(r[sel].mean())
        em = float(np.sqrt(np.sum(e[sel] ** 2)) / sel.sum())
        print(f"   z_h = {zh * 1e3:g} µm, mean over z >= 15 µm {m:.3f} ± {em:.3f}, per depth:")
        print("   " + " ".join(f"{x:.3f}" for x in r[sel]))
    print("\nRuntimes (s):")
    for f in sorted(OUT.glob("lauf_*_[0-9].npz")):
        print(f"   {f.stem} {float(np.load(f)['runtime_s']):.0f}")


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("cal")
    r = sub.add_parser("run")
    r.add_argument("--set", dest="satz", choices=["A", "B", "C"], required=True)
    r.add_argument("--zh", type=float, required=True, help="hand-off depth, mm")
    r.add_argument("--seed", type=int, required=True)
    r.add_argument("--packets", type=int, default=10_000_000)
    r.add_argument("--pilot", action="store_true", help="suffix _pilot, ignored by eval")
    sub.add_parser("eval")
    a = ap.parse_args()
    if a.cmd == "cal":
        calibrate()
    elif a.cmd == "run":
        run(a.satz, a.zh, a.seed, a.packets, a.pilot)
    else:
        evaluate()


if __name__ == "__main__":
    main()
