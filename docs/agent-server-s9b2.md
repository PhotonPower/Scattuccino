# Auftrag für einen Agenten: Kalibrierung des Zufallsmediums in BWPM (S9b-2)

> Für einen Programmieragenten auf dem GPU-Rechner. Es gelten `AGENTS.md` und die
> Abhängigkeitsregeln. Hintergrund: `docs/entwurf-s9b-hybrid.md` (D9), Testfall T16.

## Ziel

Für realistische Gewebeparameter zeigen, ob ein Zufallsfeld δn (Whittle-Matérn) in der BWPM die
Streuparameter liefert, die der Monte Carlo darunter verwendet: μs aus der Extinktion des
kohärenten Feldes, g aus der Winkelverteilung des gestreuten Lichts – jeweils gegen die
Born-Näherung (`scattuccino.microstructure.born`). Ergebnis ist eine Tabelle, aus der für S9b-3
die Parameter (σ_n, l_c, D) zu gegebenem (μs, g) abgelesen werden können, und die Aussage, wo Born
nicht mehr reicht.

## Regeln

- Branch `s9b2-medium`; nichts im Code ändern außer neuen Skripten unter `tools/` (z. B.
  `tools/calibrate_bwpm.py`) und dem Bericht; keine Tests oder Toleranzen ändern.
- BWPM mit CUDA-Kern (wie in BWPM B4d gebaut) installiert: `python -c "from bwpm import fast;
  print(fast.CORE_AVAILABLE, fast.cuda_available())"` → `True True`.
- Pro Lauf Laufzeit und Speicher notieren; Läufe über 30 min vorher abschätzen.

## Läufe

`scattuccino.microstructure.calibrate_with_bwpm(medium, n0=1.38, wavelength_um=0.8,
shape=(nz, n, n), spacing, rng, keep_every, backend="cuda")`, Gitter isotrop mit
`spacing = λ/(5 n0)` (≈ 0,116 µm), lateral `n = 1024` (≈ 119 µm), Tiefe bis μs·z ≈ 4.

| Gruppe | Parameter | Zweck |
| --- | --- | --- |
| A | D = 4, l_c = 0,5 µm; σ_n so, dass Born μs = 10, 20, 50 /mm | Extinktion, g bei kleinem σ_n |
| B | D = 4, μs (Born) = 20 /mm; l_c = 0,2 / 0,5 / 1,0 / 2,0 µm | g über l_c (Born: g steigt mit l_c) |
| C | l_c = 0,5 µm, μs (Born) = 20 /mm; D = 3,5 / 4 / 5 | Form der Phasenfunktion (γ) |
| D | wie A mit 3 Realisierungen (Seeds 1–3) | Streuung zwischen Realisierungen |

Für Gruppe A bei 10 /mm ist die Tiefe 400 µm ≈ 3 450 Ebenen; schätze vorher Zeit und Speicher
(1024² × Ebenen in FP64 auf der GPU) und nutze `keep_every` so, dass etwa 20 Ebenen gespeichert
werden. Reicht der Speicher nicht, lateral 768 oder 512 nehmen und das im Bericht vermerken.

Je Lauf berichten: σ_n, l_c, D, Born-μs, gefittetes μs (Verhältnis), Born-g (vorwärts),
gemessenes g (vorwärts), Leistungserhaltung (max |total − 1|), Laufzeit.

## Zusätzlich

Für einen Lauf aus Gruppe A die Winkelverteilung p(θ) aus der BWPM (vorwärts, normiert) gegen
Born als Datei `tools/out/phase_bwpm_vs_born.csv` (θ, p_BWPM, p_Born) ablegen.

## Bericht

PR „S9b-2: calibration runs of the random medium in BWPM“ (Basis `main`, Branch `s9b2-medium`)
mit Tabelle, Umgebung, Skript und Befund: Wo weicht μs oder g um mehr als 10 % von Born ab?
