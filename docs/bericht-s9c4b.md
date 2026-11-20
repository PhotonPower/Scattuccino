# Bericht S9c-4b: Rückstreuung der BWPM an einem Zufallsmedium (Ensemble)

Auftrag: Issue 49. Werkzeug `tools/backscatter_bwpm.py` (unverändert), Rohdaten
`tools/out/backscatter_runs.jsonl` (90 Zeilen).

## Umgebung

RTX 6000 Ada, nvcc 12.9, Python 3.13.5. Scattuccino `main` `da1fd09`; BWPM `main` `3153738`
(enthält B8 `9b9ca00`). Abweichung vom Auftrag: `pip install -e .` in BWPM baut ohne CUDA
(danach `CORE_AVAILABLE = False`); BWPM wurde deshalb mit dem Befehl aus
`docs/agent-server-b4d-confirm.md` gebaut (gcc-toolset-15, `BWPM_ENABLE_CUDA=ON`,
`BWPM_STATIC_LIBSTDCXX=ON`, `CMAKE_CUDA_HOST_COMPILER=/usr/bin/g++`), danach
`HAS_CUDA, CORE_AVAILABLE, cuda_available = True True True`. Tests `test_s9c4_backscatter.py` und
`test_s9b2_microstructure.py`: 9 bestanden (darunter `test_bwpm_back_scattering_against_born`, nicht
übersprungen), 2 erwartete xfail (`test_bwpm_realisation_reproduces_born_g`,
`test_bwpm_born_increment_conserves_power`).

## Laufzeiten

Pilotläufe: Satz A (h = 0,025 µm, Seed 1) 101 s, Satz B (Seed 1) 189 s (jeweils mit Start); daraus
höchstens 2,8 h für die 80 BWPM-Läufe geschätzt (Grenze 10 h). Die Pilotzeilen wurden aus
`backscatter_runs.jsonl` entfernt, da die Schleife dieselben Läufe (deterministisch) wiederholt.
Gemessene Rechenzeit je Lauf (`time_s`):

| Lauf | Anzahl | je Lauf | Summe |
| --- | --- | --- | --- |
| BWPM A, h = 0,1 µm | 20 | 1–2 s | 25 s |
| BWPM A, h = 0,05 µm | 20 | 4–8 s | 106 s |
| BWPM A, h = 0,025 µm | 20 | 40–59 s | 947 s |
| BWPM B, h = 0,05 µm | 20 | 105–127 s | 2294 s |
| Monte Carlo A (2·10⁷ Pakete) | 5 | 175–179 s | 879 s |
| Monte Carlo B (4·10⁶ Pakete) | 5 | 36 s | 180 s |

GPU-Speicher habe ich nicht gemessen; Satz B (256² × 800 Scheiben) lief ohne Fehler.

## Ausgabe von `table`

```
|cos θ| bins: 0.35 0.45 0.55 0.65 0.75 0.85
set A MC / continuous Born: 1.081 1.071 1.062 1.025 1.036 1.070
set A h 0.025 µm, 20 realisations
  BWPM / continuous Born:    1.167 1.129 1.337 1.218 1.083 2.119
  ± (standard error)         0.028 0.023 0.022 0.022 0.019 0.032
  discrete Born / continuous: 1.140 1.081 1.315 1.203 1.065 2.119
  BWPM / discrete Born:      1.024 1.045 1.017 1.013 1.017 1.000
  BWPM / MC:                 1.079 1.055 1.259 1.189 1.045 1.980
set A h 0.050 µm, 20 realisations
  BWPM / continuous Born:    1.138 1.021 1.293 1.203 1.054 2.108
  ± (standard error)         0.036 0.019 0.023 0.013 0.013 0.034
  discrete Born / continuous: 1.108 0.984 1.284 1.184 1.036 2.114
  BWPM / discrete Born:      1.027 1.038 1.007 1.016 1.018 0.997
  BWPM / MC:                 1.052 0.954 1.218 1.174 1.018 1.970
set A h 0.100 µm, 20 realisations
  BWPM / continuous Born:    1.080 0.877 1.149 1.006 0.882 1.834
  ± (standard error)         0.032 0.022 0.015 0.019 0.013 0.024
  discrete Born / continuous: 1.029 0.841 1.137 1.000 0.864 1.862
  BWPM / discrete Born:      1.050 1.043 1.011 1.006 1.021 0.985
  BWPM / MC:                 0.999 0.818 1.082 0.982 0.851 1.713
set B MC / continuous Born: 3.145 2.699 2.333 2.023 1.893 1.761
set B h 0.050 µm, 20 realisations
  BWPM / continuous Born:    3.757 2.917 2.462 2.124 1.903 1.770
  ± (standard error)         0.118 0.091 0.046 0.035 0.023 0.026
  discrete Born / continuous: 2.073 1.796 1.741 1.614 1.574 1.542
  BWPM / discrete Born:      1.813 1.624 1.414 1.316 1.209 1.148
  BWPM / MC:                 1.195 1.081 1.055 1.050 1.006 1.005
```

## Beobachtung

**Satz A (Einfachstreuung, L = 2 µm).** BWPM gegen diskretes Born stimmt bei h = 0,025 und 0,05 µm
auf 0–4,5 % überein (Verhältnisse 1,000–1,045, Standardfehler der BWPM/kontinuierlich-Verhältnisse
2–3 %) und bei h = 0,1 µm auf 0,6–5 %. Gegen das kontinuierliche Born liegt die BWPM bei
h = 0,025 µm zwischen 1,08 und 1,34 und im letzten Bin (|cos θ| = 0,85) bei 2,12; dieser Wert und
die Abhängigkeit von h folgen dem diskreten Born (2,12 gegen 2,11 bei h = 0,05 µm), die Verhältnisse
zum kontinuierlichen Born ändern sich von h = 0,1 zu 0,05 µm stark (0,88–1,83 → 1,02–2,11) und von
0,05 zu 0,025 µm um höchstens 0,11 (größte Änderung im zweiten Bin: 1,021 → 1,129). Der
Monte Carlo liegt 2,5–8 % über dem kontinuierlichen Born.

**Satz B (Mehrfachstreuung, L = 40 µm, h = 0,05 µm).** BWPM gegen kontinuierliches Born ist 3,76,
2,92, 2,46, 2,12, 1,90, 1,77 in den Bins von |cos θ| = 0,35 bis 0,85, gegen das diskrete Born
(Einfachstreuung) 1,81 bis 1,15; gegen den Monte Carlo liegt sie bei 1,20, 1,08, 1,06, 1,05, 1,01,
1,01, also nur im ersten Bin (|cos θ| = 0,35) um mehr als 10 % darüber. Der Monte Carlo selbst liegt
bei 3,15 bis 1,76 über dem kontinuierlichen Born.
