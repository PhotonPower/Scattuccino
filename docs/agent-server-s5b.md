# Auftrag für einen Agenten: CUDA-Kern (S5b) auf dem Server bauen, testen, messen

> Für einen Programmieragenten auf dem Rechner des Maintainers (RTX 6000 Ada). Es gelten
> `AGENTS.md`, `docs/konventionen.md`, `docs/abhaengigkeiten.md`. Erfahrungen aus BWPM B4d:
> nvcc 12.9 mit gcc-toolset-15 als Host-Compiler, statische libstdc++ wegen conda-Python.

## Ziel

1. `scattuccino._core` mit CUDA bauen (`from scattuccino import _core` lädt ohne Umgebungsvariablen).
2. Tests grün: alle, dazu `tests/test_s5_cuda.py` mit `SCATTUCCINO_REQUIRE_CUDA=1` – ohne Tests oder
   Toleranzen zu ändern.
3. Laufzeiten messen (CPU-Kern auf allen Kernen gegen CUDA) und berichten.

## Regeln

- Branch `s5b-cuda`; nicht mergen, nicht auf `main` pushen.
- Erlaubt: Build-Korrekturen (`CMakeLists.txt`, `pyproject.toml`) und Fehlerbehebungen in
  `cpp/src/core_cuda.cu`, jeweils eigener Commit mit Begründung. Nicht ändern: Tests, Toleranzen,
  `cpp/include/scattuccino/photon.hpp` (gemeinsamer Paketcode für CPU und GPU – Änderungen dort
  verändern auch den CPU-Kern; bei Bedarf anhalten und berichten).
- Keine Pfade, Hostnamen oder Tokens in Dateien oder Commits.

## Schritte

```bash
cd ~/Scattuccino && git fetch && git checkout s5b-cuda && git reset --hard origin/s5b-cuda
nvidia-smi --query-gpu=name,driver_version --format=csv && which nvcc && nvcc --version | tail -1
rm -rf build
CXX=<g++-15> CC=<gcc-15> pip install -e . -C cmake.define.SCATTUCCINO_ENABLE_CUDA=ON \
  -C cmake.define.SCATTUCCINO_STATIC_LIBSTDCXX=ON -C cmake.define.CMAKE_CUDA_HOST_COMPILER=<g++-15>
python -c "from scattuccino import _core, fast; print(_core.HAS_CUDA, _core.cuda_device_count(), fast.cuda_available())"
pytest -q
SCATTUCCINO_REQUIRE_CUDA=1 pytest tests/test_s5_cuda.py tests/test_s5_core.py -v
```

Erwartet: `True 1 True`, alles grün. Häufige Fehlerbilder (aus BWPM): nvcc lehnt den
Host-Compiler ab → gcc-toolset-15 oder `-allow-unsupported-compiler`; `GLIBCXX_…` beim Import →
`SCATTUCCINO_STATIC_LIBSTDCXX=ON`; Übersetzungsfehler in `core_cuda.cu` (nur mit clang geprüft) →
beheben, eigener Commit.

**Schlägt `test_cuda_equals_cpu_packet_by_packet` fehl:** den Anteil gleicher Pakete und die
Abweichung der Summen ausgeben (siehe Test) und berichten; die GPU-Mathebibliothek rundet
log/cos/sin/pow/cbrt in den letzten Bits anders, gefordert sind ≥ 99,9 % gleiche Pakete.

## Messung

```bash
python - <<'PY'
import time, numpy as np
from scattuccino import fast
from scattuccino.layers import Layer
from scattuccino.voxels import VoxelMedium
cases = {
    "Haut-Wuerfel 101^3, 1e7": (VoxelMedium(np.zeros((101, 101, 101), dtype=np.int64), (0.05, 0.05, 0.05),
                                (Layer(0.0, 0.05, 10.0, 0.9, 1.4),), lateral="open"), 10_000_000),
    "schwach absorbierend, 1e7": (VoxelMedium(np.zeros((1, 1, 200), dtype=np.int64), (100.0, 100.0, 0.1),
                                  (Layer(0.0, 0.01, 10.0, 0.9, 1.4),)), 10_000_000),
}
for name, (med, n) in cases.items():
    for backend in ("cpu", "cuda"):
        t = time.time()
        s = fast.simulate_voxels(med, n, seed=1, backend=backend)
        print(name, backend, f"{time.time() - t:.1f} s", "R =", s.protocol.weight[s.protocol.side == -1].sum())
PY
```

## Bericht

Als Kommentar in PR (S5b) oder als Datei: Umgebung, Build-Befehl, Testergebnis (Zahlen),
Messtabelle (Fall, Backend, Laufzeit, R), Commits mit Begründung.
