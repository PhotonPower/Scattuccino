"""Calibration runs of the random medium in BWPM (stage S9b-2, docs/agent-server-s9b2.md).

One run propagates a plane wave through one realisation of a Whittle-Matérn index field
(`scattuccino.microstructure`) in the forward FPWPM of BWPM (CUDA core) and measures

* μs from the decay of the coherent field |⟨E⟩|² (fit as in `calibrate_with_bwpm`),
* the angular distribution of the scattered field (forward half) and its g at every kept depth,
* the power balance (mean of |E|² over the plane) at every kept depth,

against the first-order Born approximation. The analysis repeats `calibrate_with_bwpm` and adds
the series over depth (the library reports only the first kept depth; for a plane wave the
scattered spectrum builds up over about 2/(k(1 − cos θ)) and is multiply scattered deeper down).

Units: lengths in mm in the code (grid spacing, l_c), µm only for the wavelength and in the output.
Every run appends one JSON line to tools/out/calibration_runs.jsonl.

    python tools/calibrate_bwpm.py run --mus 20 --lc 0.5 --d 4 --seed 1 --n 512 --levels 16
    python tools/calibrate_bwpm.py table            # markdown table of all runs
"""

from __future__ import annotations

import argparse
import json
import resource
import subprocess
import threading
import time
from pathlib import Path

import numpy as np

import scattuccino.microstructure as ms

OUT = Path(__file__).parent / "out"
RUNS = OUT / "calibration_runs.jsonl"
N0, WAVELENGTH_UM = 1.38, 0.8


def sigma_for_mus(mus: float, l_c: float, d: float) -> float:
    """σ_n with Born μs = `mus` (1/mm); μs ∝ σ_n² exactly."""
    unit = ms.born(ms.WhittleMatern(1.0, l_c, d), N0, WAVELENGTH_UM).mus
    return float(np.sqrt(mus / unit))


class GpuPeak:
    """Peak GPU memory (MiB) of this process, polled with nvidia-smi."""

    def __init__(self) -> None:
        self.peak = 0
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._poll, daemon=True)

    def _poll(self) -> None:
        import os

        pid = str(os.getpid())
        while not self._stop.is_set():
            out = subprocess.run(
                [
                    "nvidia-smi",
                    "--query-compute-apps=pid,used_memory",
                    "--format=csv,noheader,nounits",
                ],
                capture_output=True,
                text=True,
            ).stdout
            for line in out.splitlines():
                p, _, mem = line.partition(",")
                if p.strip() == pid:
                    self.peak = max(self.peak, int(mem))
            self._stop.wait(1.0)

    def __enter__(self) -> GpuPeak:
        self._thread.start()
        return self

    def __exit__(self, *exc: object) -> None:
        self._stop.set()
        self._thread.join()


def angular_distribution(
    f: np.ndarray,
    spacing: float,
    n0: float,
    edges: np.ndarray,
) -> tuple[float, float, np.ndarray]:
    """Forward scattered field f (ny, nx, 2) on a plane → (g as in the library, g with the
    solid-angle Jacobian, p(θ) per solid angle on the θ bins `edges`, normalised over the forward
    hemisphere).

    Far field amplitude ∝ cos θ · F(k⊥), so the power per solid angle is |F|² cos²θ and a cell
    of the k⊥ grid covers dΩ = dk⊥²/(k² cos θ). The library's g weights the cells with
    |F|² cos²θ (no Jacobian); the second value uses the cell weight |F|² cos²θ · dΩ.
    """
    ny, nx, _ = f.shape
    k = 2 * np.pi * n0 / (WAVELENGTH_UM * 1e-3)
    kx = 2 * np.pi * np.fft.fftfreq(nx, spacing)
    ky = 2 * np.pi * np.fft.fftfreq(ny, spacing)
    kxx, kyy = np.meshgrid(kx, ky)
    kp2 = kxx**2 + kyy**2
    prop = kp2 < k**2
    ct = np.sqrt(np.clip(1 - kp2 / k**2, 1e-12, 1))
    s = np.zeros((ny, nx))
    for comp in range(2):
        s += np.abs(np.fft.fft2(f[:, :, comp] - np.mean(f[:, :, comp]))) ** 2
    w = s * ct**2 * prop
    g_lib = float(np.sum(w * ct) / np.sum(w))
    dk = (kx[1] - kx[0]) * (ky[1] - ky[0])
    d_omega = dk / (k**2 * ct)
    pw = w * d_omega
    g_jac = float(np.sum(pw * ct) / np.sum(pw))
    theta = np.arccos(ct)
    p_bin, _ = np.histogram(theta[prop], bins=edges, weights=pw[prop])
    omega = 2 * np.pi * (np.cos(edges[:-1]) - np.cos(edges[1:]))
    p = p_bin / omega
    p /= np.sum(p * omega)
    return g_lib, g_jac, p


def run_one(
    mus: float,
    l_c_um: float,
    d: float,
    seed: int,
    n: int,
    levels: int,
    depth_mus: float,
    phase_csv: Path | None = None,
) -> dict[str, object]:
    from bwpm import fast
    from bwpm.volume3d import Bidirectional3D, Grid3D

    l_c = l_c_um * 1e-3
    spacing = WAVELENGTH_UM / (5 * N0) * 1e-3  # mm
    sigma = sigma_for_mus(mus, l_c, d)
    med = ms.WhittleMatern(sigma, l_c, d)
    nz = int(np.ceil(depth_mus / mus / spacing))
    keep_every = max(nz // 20, 1)
    keep = tuple(range(0, nz + 1, keep_every))
    t0 = time.time()
    dn = ms.random_field(med, (nz, n, n), spacing, np.random.default_rng(seed))
    lab, idx = ms.quantise(dn, N0, levels)
    var_ratio = float(np.var(dn) / sigma**2)
    q_err = float(np.std(np.asarray(idx)[lab] - N0 - dn))  # rms quantisation error
    del dn
    t_field = time.time() - t0
    solver = Bidirectional3D(
        Grid3D(n, n, spacing, spacing, spacing, WAVELENGTH_UM), idx, lab, complex(N0), None
    )
    del lab
    e = np.zeros((n, n, 2), dtype=np.complex128)
    e[:, :, 0] = 1.0
    with GpuPeak() as gpu:
        t1 = time.time()
        res = fast.run_3d(solver, e, iterations=0, keep=keep, backend="cuda")
        t_run = time.time() - t1
    depths = np.array(keep, dtype=np.float64) * spacing
    coh, tot = [], []
    for plane in keep:
        f = res.planes[plane][0]
        coh.append(float(np.abs(np.mean(f[:, :, 0])) ** 2 + np.abs(np.mean(f[:, :, 1])) ** 2))
        tot.append(float(np.mean(np.sum(np.abs(f) ** 2, axis=2))))
    c = np.array(coh)
    sel = (c < 0.9) & (c > 0.05)
    sel[0] = False
    mus_fit = (
        float(-np.polyfit(depths[sel], np.log(c[sel]), 1)[0]) if sel.sum() >= 2 else float("nan")
    )
    edges = np.linspace(0, np.pi / 2, 61)
    g_lib, g_jac, p_all = [], [], []
    for plane in keep[1:]:
        gl, gj, p = angular_distribution(res.planes[plane][0], spacing, N0, edges)
        g_lib.append(gl)
        g_jac.append(gj)
        p_all.append(p)
    b = ms.born(med, N0, WAVELENGTH_UM)
    fw = b.phase.mu > 0
    mu, pp = b.phase.mu[fw], b.phase.p[fw]
    g_born = float(np.trapezoid(mu * pp, mu) / np.trapezoid(pp, mu))
    # reference depth for g: nearest kept depth to 0.25/μs_born (15-20 % of the field scattered)
    ref = int(np.argmin(np.abs(depths[1:] - 0.25 / b.mus)))
    if phase_csv is not None:
        th = 0.5 * (edges[:-1] + edges[1:])
        mu_th = np.cos(th)
        order = np.argsort(b.phase.mu)
        p_born = np.interp(mu_th, b.phase.mu[order], b.phase.p[order])
        p_born /= np.sum(p_born * 2 * np.pi * (np.cos(edges[:-1]) - np.cos(edges[1:])))
        phase_csv.parent.mkdir(parents=True, exist_ok=True)
        with open(phase_csv, "w") as fh:
            fh.write(
                "# theta_rad, p_bwpm, p_born: power per solid angle, forward hemisphere, "
                f"normalised to 1 over the hemisphere; BWPM at depth "
                f"{depths[1 + ref] * 1e3:.1f} um\n"
            )
            fh.write("theta,p_bwpm,p_born\n")
            for a, bw, bo in zip(th, p_all[ref], p_born, strict=True):
                fh.write(f"{a:.6f},{bw:.8e},{bo:.8e}\n")
    row: dict[str, object] = {
        "mus_born": b.mus,
        "sigma_n": sigma,
        "l_c_um": l_c_um,
        "D": d,
        "seed": seed,
        "n_lateral": n,
        "nz": nz,
        "levels": levels,
        "spacing_um": spacing * 1e3,
        "depth_um": nz * spacing * 1e3,
        "keep_every": keep_every,
        "mus_fit": mus_fit,
        "mus_ratio": mus_fit / b.mus,
        "g_born": g_born,
        "g_bwpm_ref": g_lib[ref],
        "g_bwpm_ref_jacobian": g_jac[ref],
        "g_ref_depth_um": depths[1 + ref] * 1e3,
        "g_bwpm_first_library": g_lib[0],
        "g_first_depth_um": depths[1] * 1e3,
        "depths_um": (depths * 1e3).tolist(),
        "coherent": c.tolist(),
        "total": tot,
        "g_series": g_lib,
        "g_series_jacobian": g_jac,
        "power_err": float(np.max(np.abs(np.array(tot) - 1.0))),
        "var_ratio_field": var_ratio,
        "quantisation_rms": q_err,
        "time_field_s": t_field,
        "time_run_s": t_run,
        "time_total_s": time.time() - t0,
        "gpu_peak_mib": gpu.peak,
        "host_maxrss_gb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6,
        "n_fit_planes": int(sel.sum()),
    }
    OUT.mkdir(exist_ok=True)
    with open(RUNS, "a") as fh:
        fh.write(json.dumps(row) + "\n")
    return row


def table() -> None:
    rows = [json.loads(line) for line in open(RUNS)]
    print(
        "| σ_n | l_c µm | D | seed | n | L | Born μs | fit μs | ratio | Born g | g (ref) | "
        "g first | max|tot−1| | s |"
    )
    print("|" + " --- |" * 14)
    for r in rows:
        print(
            f"| {r['sigma_n']:.4f} | {r['l_c_um']} | {r['D']} | {r['seed']} | {r['n_lateral']} | "
            f"{r['levels']} | {r['mus_born']:.2f} | {r['mus_fit']:.2f} | {r['mus_ratio']:.3f} | "
            f"{r['g_born']:.4f} | {r['g_bwpm_ref']:.4f} | {r['g_bwpm_first_library']:.4f} | "
            f"{r['power_err']:.1e} | {r['time_total_s']:.0f} |"
        )


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--mus", type=float, required=True, help="Born μs, 1/mm")
    r.add_argument("--lc", type=float, default=0.5, help="correlation length, µm")
    r.add_argument("--d", type=float, default=4.0)
    r.add_argument("--seed", type=int, default=1)
    r.add_argument("--n", type=int, default=512, help="lateral grid points")
    r.add_argument("--levels", type=int, default=16, help="index quantisation levels")
    r.add_argument("--depth-mus", type=float, default=4.0, help="depth in units of 1/μs")
    r.add_argument("--phase-csv", type=Path, default=None)
    sub.add_parser("table")
    a = ap.parse_args()
    if a.cmd == "table":
        table()
        return
    row = run_one(a.mus, a.lc, a.d, a.seed, a.n, a.levels, a.depth_mus, a.phase_csv)
    print({k: v for k, v in row.items() if not isinstance(v, list)})


if __name__ == "__main__":
    main()
