# Auftrag für einen Agenten: Laufzeit des CUDA-Kerns nach S5c messen

> Vorgänger: `docs/agent-server-s5b.md` (Build-Befehl dort, Bericht in PR #11). Es gelten
> `AGENTS.md` und die Regeln aus dem S5b-Auftrag.

**Änderungen S5c:** persistente Threads (jeder Thread holt sich das nächste Paket aus einem
globalen Zähler), Absorption je Material in Registern statt atomarer Additionen bei jeder
Wechselwirkung, Voxelwertung und Protokoll abschaltbar (`score_voxels=False`, `record=False`).

## Schritte

```bash
cd ~/Scattuccino && git fetch && git checkout s5b-cuda && git reset --hard origin/s5b-cuda && rm -rf build
T=/opt/rh/gcc-toolset-15/root/usr/bin
CXX=$T/g++ CC=$T/gcc pip install -e . -C cmake.define.SCATTUCCINO_ENABLE_CUDA=ON \
  -C cmake.define.SCATTUCCINO_STATIC_LIBSTDCXX=ON -C cmake.define.CMAKE_CUDA_HOST_COMPILER=/usr/bin/g++
SCATTUCCINO_REQUIRE_CUDA=1 pytest tests/test_s5_cuda.py tests/test_s5_core.py -v
python - <<'PY'
import time, numpy as np
from scattuccino import fast
from scattuccino.layers import Layer
from scattuccino.voxels import VoxelMedium
cases = {
    "Haut-Wuerfel 101^3": VoxelMedium(np.zeros((101, 101, 101), dtype=np.int64), (0.05, 0.05, 0.05),
                                      (Layer(0.0, 0.05, 10.0, 0.9, 1.4),), lateral="open"),
    "schwach absorbierend": VoxelMedium(np.zeros((1, 1, 200), dtype=np.int64), (100.0, 100.0, 0.1),
                                        (Layer(0.0, 0.01, 10.0, 0.9, 1.4),)),
}
for name, med in cases.items():
    for backend in ("cpu", "cuda"):
        for record, voxels in ((True, True), (False, False)):
            t = time.time()
            s = fast.simulate_voxels(med, 10_000_000, seed=1, backend=backend, record=record,
                                     score_voxels=voxels)
            print(name, backend, "Protokoll" if record else "nur Summen", f"{time.time() - t:.1f} s",
                  "R =", s.totals["r_diffuse"])
PY
```

## Bericht

Kommentar in PR #11: Testergebnis, Tabelle (Fall, Backend, Variante, Laufzeit, R).
