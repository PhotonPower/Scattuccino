# Auftrag für einen Agenten: FP32-Modus für den Voxeltransport (S5d)

> Für einen Programmieragenten auf dem GPU-Rechner des Maintainers (RTX 6000 Ada, 32 Kerne).
> Es gelten `AGENTS.md`, `docs/konventionen.md`, `docs/abhaengigkeiten.md`. Build wie in
> `docs/agent-server-s5c.md` (gcc-toolset-15 für C++ und Linken, System-gcc 14 als
> nvcc-Host-Compiler, statische libstdc++). Entscheidung D7 (Issue #14): FP64 bleibt Standard und
> Referenz, FP32 kommt als Option `precision="single"`, Abnahme statistisch.

## Ausgangslage

S5c-Messung (PR #11): CUDA 13,3 s / 62,9 s gegen 25,4 s / 169,6 s auf 32 CPU-Kernen (10⁷ Pakete,
Würfel 101³ / schwach absorbierend). Der verbleibende Engpass ist FP64: Die RTX 6000 Ada rechnet
FP64 mit 1/64 des FP32-Durchsatzes, und der Paketcode (`cpp/include/scattuccino/photon.hpp`)
rechnet durchgehend in double – auch `log`, `sin`, `cos`, `pow`, `cbrt`, `sqrt` und alle
Divisionen.

## Ziel

1. `precision="single"` in `scattuccino.fast.simulate_voxels` für **CPU und CUDA** (CPU in
   single dient als Referenz für die GPU in single).
2. Mit `precision="double"` (Standard) ändert sich **nichts**: Der erzeugte Code für double und
   damit alle bestehenden Ergebnisse bleiben bitgleich (T13, T14 grün ohne Änderung).
3. Neue Tests T15 (unten) grün; Laufzeiten gemessen und berichtet.

## Regeln

- Branch `s5d-fp32` (enthält diesen Auftrag); nicht mergen, nicht auf `main` pushen; Commits mit
  Begründung, ein Thema je Commit.
- **Nicht ändern:** bestehende Tests und Toleranzen, `scattuccino/voxels.py` (Python-Referenz),
  `photon_pol.hpp`, `photon_tet.hpp` (Polarisation und Tetraeder kommen später).
- **Erlaubt:** `photon.hpp`, `core.cpp`, `core_cuda.cu`, `core.hpp`, `bindings.cpp`,
  `scattuccino/fast.py`, `CMakeLists.txt`, neue Testdatei `tests/test_s5d_single.py`,
  `docs/testfaelle.md` (Abschnitt T15), `CHANGELOG.md`.
- Keine `--use_fast_math` und keine `__logf`/`__sinf`-Intrinsics im Standardweg (Genauigkeit vor
  Tempo); wenn du sie als zusätzliche Variante misst, nur im Bericht, nicht im Code lassen.
- Keine Pfade, Hostnamen oder Tokens in Dateien oder Commits.

## Umsetzung (Vorgabe)

1. **Präzision als Template-Parameter** `template <class Real>` für `trace`, `sample`, `hg`,
   `fresnel`, `scatter`, `Geometry` (Gitterweiten, Ursprung) und `MaterialView`/`PhaseView`
   (Koeffizienten, tabellierte Arrays). Der bisherige Code wird zu `Real = double`; Aufrufer im
   double-Weg ändern sich nur in der Template-Angabe. Prüfe nach dem Umbau zuerst, dass T13 und
   T14 unverändert grün sind, bevor du float einführst.
2. **Was in welcher Genauigkeit:**
   - `Real` (float im Modus single): Richtung, Schrittlänge, Abstände zu den Voxelflächen,
     Phasenfunktion (Ziehen von cos θ inklusive tabellierter Arrays als float-Kopie), Fresnel,
     Streurichtung, `log`/`sin`/`cos`/`pow`/`cbrt`/`sqrt` (die float-Überladungen aus `<cmath>`
     bzw. die CUDA-Funktionen `logf` usw.).
   - **immer double:** Gewicht des Pakets, absorbierte Energie, alle Summen, Protokollwerte im
     Ergebnis, und die Zufallszahl selbst (53 Bit), die erst bei der Verwendung in `Real`
     umgewandelt wird – so bleibt der Zufallsstrom je Paket derselbe.
   - **Position:** zwei Varianten messen – (A) Position in `Real`, (B) Position in double, Rest
     in float. Bei (A) auf Verlust an Voxelflächen achten (Pakete, die eine Fläche wegen Rundung
     nicht erreichen, endlos laufen oder aus dem Gitter fallen): Der Abstand wird schon jetzt mit
     `dmax(d, 0)` begrenzt; prüfe zusätzlich, dass der Voxelindex nach jedem Flächenschritt zur
     Position passt. Behalte im Code die Variante, die T15 sicher besteht und schneller ist;
     begründe die Wahl im Bericht.
3. **Python:** `fast.simulate_voxels(..., precision="double" | "single")`; `"single"` mit
   `backend="python"` → `ValueError` (die Python-Referenz bleibt double).
4. **CUDA:** dieselben Templates; der Kernel wird für `float` und `double` instanziiert; die
   Compiler-Optionen `--fmad=false --prec-div=true --prec-sqrt=true` bleiben für double; für
   float dürfen sie bleiben (Standard) – die Wirkung von `--fmad=true` nur für float als
   Messvariante im Bericht.

## Tests T15 (neu, `tests/test_s5d_single.py`)

Medien wie in `tests/test_s5_core.py` (`_inclusion`, `_layers`) plus eine schwach absorbierende
Schicht (μa = 0,01, μs = 10 /mm, g = 0,9, n = 1,4, 5 mm, 1 Voxel lateral periodisch).

1. **double unverändert:** `precision="double"` liefert für `_inclusion` und `_layers` dieselben
   Gewichte wie die Python-Referenz bitgleich (wie T13).
2. **single gegen double, statistisch** (CPU, je 2·10⁵ Pakete, verschiedene Seeds): R_d, T_d,
   seitliche Verluste und Absorption je Material innerhalb 4σ des kombinierten Standardfehlers;
   R(ρ) auf 8 Ringen (0–4 mm) für `_inclusion` innerhalb 4σ je Ring.
3. **Energie in single:** r_specular + Σ Austritte + Absorption = 1 auf 10⁻⁴ (Roulette).
4. **GPU single gegen CPU single, je Paket** (nur mit CUDA, sonst überspringen; mit
   `SCATTUCCINO_REQUIRE_CUDA=1` Pflicht): ≥ 99 % der Pakete mit gleicher Seite und gleicher Zahl
   der Streuungen; Summen auf 10⁻³ relativ.
5. **Kein Paket verloren:** Zahl der Austritte + Roulette-Abbrüche = Pakete (über das Protokoll
   und einen Zähler, falls nötig ergänzen).

## Messung

Wie im S5c-Auftrag (Würfel 101³ und schwach absorbierend, 10⁷ Pakete, `record=False`,
`score_voxels=False`), jeweils CPU double, CPU single, CUDA double, CUDA single (Variante A und B),
dazu R mit Standardfehler.

## Bericht

Kommentar im PR zu `s5d-fp32` (lege ihn an, Basis `s5b-cuda`, Titel
„S5d: single-precision option (GPU)“): Umgebung, Build, Testergebnisse mit Zahlen
(Abweichungen in σ), Messtabelle, gewählte Positionsvariante mit Begründung, Liste der Commits.
