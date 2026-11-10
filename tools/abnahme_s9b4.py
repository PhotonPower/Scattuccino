"""Acceptance runs of the hybrid OCT: independence of the hand-off depth z_h (issue 29, S9b-4).

    python tools/abnahme_s9b4.py cal                                  # step 0: BWPM calibration
    python tools/abnahme_s9b4.py run --set A --zh 0.005 --seed 1      # one hybrid run
    python tools/abnahme_s9b4.py run --set A --zh 0.005 --seed 1 --packets 1000000 --pilot
    python tools/abnahme_s9b4.py eval                                 # tables and verdict

Units: mm for lengths, µm for the wavelength. Outputs go to tools/out/abnahme_s9b4/.
Set A uses the scattering measured in BWPM (from `cal`), set B the Born default.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from scattuccino.hybrid import Beam, hybrid_oct
from scattuccino.microstructure import WhittleMatern, calibrate_with_bwpm
from scattuccino.phase import Tabulated

OUT = Path(__file__).parent / "out" / "abnahme_s9b4"
MED = WhittleMatern(0.0572, 0.5e-3)  # Born: μs = 200 /mm, g = 0.935, 1/μs = 5 µm
N0, WAVELENGTH_UM = 1.38, 0.8
Z = np.linspace(0.0, 0.06, 31)
ZH_A = (0.0, 0.0025, 0.005, 0.010)
ZH_B = (0.0, 0.005)
SEEDS_A = (1, 2, 3)
Z_MIN = 0.015  # evaluate depths from here (behind all hand-offs)


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
    sc = load_scattering() if satz == "A" else None
    t0 = time.time()
    res = hybrid_oct(
        MED, N0, 1.0, WAVELENGTH_UM, Beam(2e-3, 0.012), z_h, Z, 0.004, packets, seed=seed,
        nx=256, dx=1.5e-4, max_path=0.15, sigma_husimi=0.6e-3, backend="cuda", scattering=sc,
    )  # fmt: skip
    dt = time.time() - t0
    name = f"lauf_{satz}_{z_h:g}_{seed}" + ("_pilot" if pilot else "")
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez(
        OUT / f"{name}.npz", z=res.z, signal=res.signal, power_handoff=res.power_handoff,
        runtime_s=dt, packets=packets,
    )  # fmt: skip
    print(f"{name}: {dt:.0f} s, power_handoff {res.power_handoff:.4f}")


def evaluate() -> None:
    def runs(satz: str, z_h: float) -> np.ndarray:
        files = sorted(OUT.glob(f"lauf_{satz}_{z_h:g}_[0-9].npz"))
        return np.array([np.load(f)["signal"] for f in files])

    sel = Z >= Z_MIN - 1e-12
    base = runs("A", 0.0)
    rows = {zh: runs("A", zh) for zh in ZH_A}
    print("| z µm | " + " | ".join(f"z_h={zh * 1e3:g} µm" for zh in ZH_A) + " |")
    print("|" + " --- |" * (len(ZH_A) + 1))
    mean = {zh: r.mean(0) for zh, r in rows.items()}
    sem = {zh: r.std(0, ddof=1) / np.sqrt(len(r)) for zh, r in rows.items()}
    for i in np.where(sel)[0]:
        cells = [f"{mean[zh][i]:.4g} ± {sem[zh][i]:.2g}" for zh in ZH_A]
        print(f"| {Z[i] * 1e3:.0f} | " + " | ".join(cells) + " |")
    print("\nRatio to z_h = 0 (error from the seeds):\n")
    print("| z µm | " + " | ".join(f"z_h={zh * 1e3:g} µm" for zh in ZH_A[1:]) + " |")
    print("|" + " --- |" * len(ZH_A))
    bad = []
    for i in np.where(sel)[0]:
        cells = []
        for zh in ZH_A[1:]:
            r = mean[zh][i] / mean[0.0][i]
            e = r * np.hypot(sem[zh][i] / mean[zh][i], sem[0.0][i] / mean[0.0][i])
            ok = abs(r - 1) <= 3 * e or abs(r - 1) <= 0.15
            if not ok:
                bad.append((Z[i] * 1e3, zh * 1e3, r, e))
            cells.append(f"{r:.3f} ± {e:.3f}{'' if ok else ' ✗'}")
        print(f"| {Z[i] * 1e3:.0f} | " + " | ".join(cells) + " |")
    print(f"\nDeviating (z µm, z_h µm, ratio, err): {bad}")
    print("PASSED" if not bad else "NOT PASSED")
    del base
    print("\nSet B against set A (seed 1), ratio B/A:\n")
    for zh in ZH_B:
        b = runs("B", zh)
        a = np.load(OUT / f"lauf_A_{zh:g}_1.npz")["signal"]
        if len(b):
            print(f"z_h = {zh * 1e3:g} µm:", np.round(b[0][sel] / a[sel], 3).tolist())
    print("\nRuntimes (s):")
    for f in sorted(OUT.glob("lauf_*_[0-9].npz")):
        print(f.stem, f"{float(np.load(f)['runtime_s']):.0f}")


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("cal")
    r = sub.add_parser("run")
    r.add_argument("--set", dest="satz", choices=["A", "B"], required=True)
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
