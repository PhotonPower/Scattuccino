# Bericht: Abnahme S9b-6, Hybrid-OCT mit kohärenter Übergabe (enger Fokus)

Auftrag: Issue 34. Skript `tools/abnahme_s9b6.py`, Daten `tools/out/abnahme_s9b6/`.

## Ergebnis: bestanden

Alle 69 Verhältnisse zu z_h = 0 (23 Tiefen von 16 bis 60 µm × z_h = 2,5 / 5 / 10 µm, 5 Seeds)
liegen innerhalb 3σ bei 1 oder höchstens 15 % daneben; kein Wert verfehlt das Kriterium. In Lauf 1
(Husimi, Issue 29, PR #30) waren es 4 von 69 verfehlt.

Hinweis zur Aussagekraft: Ab etwa 46 µm liegen die Fehler bei 10–25 %, einzelne Verhältnisse bis
1,35 (50 µm bei z_h = 5 und 10 µm: 1,34 ± 0,17 und 1,33 ± 0,17) bestehen nur über das 3σ-Kriterium.
Ein Anstieg über 1 zeigt sich bei großen Tiefen für alle drei z_h (Mittel über 50–60 µm etwa 1,2),
bei nur 5 Seeds ist er nicht von Rauschen zu trennen; Ursache nicht untersucht.

## Umgebung

RTX 6000 Ada, nvcc 12.9, Python 3.13.5. Scattuccino `main` `a7e7610`, RTT-Fusion `9646346`
(0.3.1), BWPM `daf7134` (Branch `b5a-pilot-3d-analysis`, so vorgefunden und nicht gewechselt),
mit CUDA-Kern (`True True`). Tests `test_s9b2_microstructure`, `test_s9b3_hybrid`,
`test_s9b4_next_event`, `test_s9b6_coherent`: 17 bestanden, 1 erwarteter xfail
(`test_bwpm_realisation_reproduces_born_g`, BWPM E9/K17: g vorwärts 0,913 gegen 0,955), keine
Skips. Strahl w_f = 2 µm, Fokus 12 µm, n₀ = 1,38, μa = 1 /mm, 10⁷ Pakete je Lauf,
`max_path` = 0,15 mm, `coherent="wigner"`, 256² bei dx = 0,15 µm, Satz A/B σ_H = 0,4 µm,
Satz C σ_H = 0,8 µm.

## Kalibrierung (Schritt 0, vektoriell, 256² × 134 Ebenen, Seed 11, Laufzeit 109 s)

| Größe | Wert |
| --- | --- |
| μs Born (vektoriell) | 190,31 /mm |
| μs gefittet (BWPM) | 190,73 /mm (Verhältnis 1,002) |
| g vorwärts Born | 0,9553 |
| g vorwärts BWPM (`cal.g_forward`) | 0,9325 |
| g der gemessenen Phasenfunktion (`measured[1].g`, ganze Kugel) | 0,9263 |
| max \|total − 1\| | 0,0032 |

`mus_fit` = 190,7267900540535, identisch zu Lauf 1. Rohdaten in `kalibrierung.json` und
`kalibrierung_phase.npz`.

## Laufzeiten

Pilot: 160 s (10⁶ Pakete) und 400 s (3·10⁶ Pakete), also etwa 40 s fest plus 120 s je 10⁶ Pakete;
daraus 12,4 h für 36 Läufe geschätzt (lineare Hochrechnung des 10⁶-Laufs allein hätte 16 h
ergeben, über der Grenze von 14 h, deshalb die zusätzliche Messung). Tatsächlich 10⁷ Pakete: z_h = 0
1107–1307 s, 2,5 µm 1150–1379 s, 5 µm 1261–1445 s, 10 µm 1380–1609 s (Satz A); Satz C (σ_H = 0,8 µm)
1680–1828 s bei 5 µm und 2288–2356 s bei 10 µm; Satz B 1152–1344 s. Einzelwerte in den
`npz`-Dateien (`runtime_s`), Summe etwa 13 h.

## Signal je Tiefe (Mittel ± Standardfehler, 5 Seeds, Satz A)

| z µm | z_h=0 µm | z_h=2.5 µm | z_h=5 µm | z_h=10 µm |
| --- | --- | --- | --- | --- |
| 16 | 0.0008372 ± 1.8e-05 | 0.0008261 ± 2.3e-05 | 0.0008565 ± 4.1e-05 | 0.0007637 ± 4.2e-05 |
| 18 | 0.0006594 ± 1.3e-05 | 0.0006368 ± 1.8e-05 | 0.0006743 ± 3.5e-05 | 0.0006815 ± 3.3e-05 |
| 20 | 0.0005064 ± 5.4e-06 | 0.000496 ± 7.5e-06 | 0.0005206 ± 2.5e-05 | 0.0005203 ± 2.7e-05 |
| 22 | 0.0003989 ± 6.3e-06 | 0.000392 ± 1e-05 | 0.0004015 ± 2.4e-05 | 0.0003969 ± 2.3e-05 |
| 24 | 0.0003302 ± 7.2e-06 | 0.000318 ± 1.4e-05 | 0.0003177 ± 2.6e-05 | 0.0003178 ± 1.8e-05 |
| 26 | 0.0002597 ± 5.6e-06 | 0.0002549 ± 1.4e-05 | 0.0002648 ± 2.1e-05 | 0.000264 ± 1.6e-05 |
| 28 | 0.0001988 ± 4.7e-06 | 0.0002022 ± 8.3e-06 | 0.0002199 ± 1.4e-05 | 0.0002046 ± 1.8e-05 |
| 30 | 0.0001616 ± 2.9e-06 | 0.000159 ± 5.7e-06 | 0.0001731 ± 1.1e-05 | 0.0001521 ± 1.5e-05 |
| 32 | 0.0001335 ± 6.2e-06 | 0.0001271 ± 6.7e-06 | 0.0001312 ± 9.6e-06 | 0.0001221 ± 1.2e-05 |
| 34 | 0.0001139 ± 8.5e-06 | 0.0001099 ± 4.4e-06 | 0.0001057 ± 5e-06 | 0.0001008 ± 9.8e-06 |
| 36 | 9.099e-05 ± 6e-06 | 9.231e-05 ± 3.8e-06 | 8.718e-05 ± 3.9e-06 | 8.45e-05 ± 9.7e-06 |
| 38 | 6.972e-05 ± 4e-06 | 7.339e-05 ± 3.4e-06 | 7.508e-05 ± 2.8e-06 | 7.197e-05 ± 8.2e-06 |
| 40 | 5.927e-05 ± 3.2e-06 | 6.003e-05 ± 8.7e-07 | 6.275e-05 ± 2.9e-06 | 6.051e-05 ± 4.7e-06 |
| 42 | 5.132e-05 ± 1.9e-06 | 5.366e-05 ± 2.3e-06 | 4.843e-05 ± 1.9e-06 | 4.904e-05 ± 3.4e-06 |
| 44 | 4.188e-05 ± 2e-06 | 4.499e-05 ± 2.1e-06 | 3.885e-05 ± 2.1e-06 | 4.258e-05 ± 1.6e-06 |
| 46 | 3.542e-05 ± 1.6e-06 | 3.815e-05 ± 2.7e-06 | 3.524e-05 ± 3.1e-06 | 3.98e-05 ± 3.3e-06 |
| 48 | 3.029e-05 ± 1.7e-06 | 3.236e-05 ± 2.5e-06 | 3.566e-05 ± 3.3e-06 | 3.52e-05 ± 3.6e-06 |
| 50 | 2.421e-05 ± 1.3e-06 | 2.703e-05 ± 1.4e-06 | 3.249e-05 ± 3.6e-06 | 3.218e-05 ± 3.6e-06 |
| 52 | 2.127e-05 ± 2.1e-06 | 2.599e-05 ± 1.7e-06 | 2.575e-05 ± 2.5e-06 | 2.866e-05 ± 3.2e-06 |
| 54 | 2.023e-05 ± 1.9e-06 | 2.518e-05 ± 2.4e-06 | 2.077e-05 ± 2.6e-06 | 2.407e-05 ± 2.7e-06 |
| 56 | 1.765e-05 ± 1.7e-06 | 2.132e-05 ± 1.8e-06 | 1.943e-05 ± 1.7e-06 | 1.981e-05 ± 1.2e-06 |
| 58 | 1.446e-05 ± 1e-06 | 1.865e-05 ± 2e-06 | 1.912e-05 ± 1.4e-06 | 1.805e-05 ± 1.4e-06 |
| 60 | 1.383e-05 ± 1.1e-06 | 1.862e-05 ± 3.2e-06 | 1.641e-05 ± 1.7e-06 | 1.863e-05 ± 3e-06 |

## Verhältnis zu z_h = 0 (Fehler aus beiden Standardfehlern; ✗ = Kriterium verfehlt)

| z µm | z_h=2.5 µm | z_h=5 µm | z_h=10 µm |
| --- | --- | --- | --- |
| 16 | 0.987 ± 0.035 | 1.023 ± 0.054 | 0.912 ± 0.054 |
| 18 | 0.966 ± 0.033 | 1.023 ± 0.056 | 1.034 ± 0.053 |
| 20 | 0.979 ± 0.018 | 1.028 ± 0.050 | 1.027 ± 0.055 |
| 22 | 0.983 ± 0.030 | 1.007 ± 0.061 | 0.995 ± 0.060 |
| 24 | 0.963 ± 0.048 | 0.962 ± 0.080 | 0.962 ± 0.059 |
| 26 | 0.982 ± 0.057 | 1.020 ± 0.085 | 1.017 ± 0.067 |
| 28 | 1.017 ± 0.048 | 1.106 ± 0.075 | 1.029 ± 0.092 |
| 30 | 0.983 ± 0.040 | 1.071 ± 0.071 | 0.941 ± 0.097 |
| 32 | 0.952 ± 0.067 | 0.982 ± 0.085 | 0.915 ± 0.097 |
| 34 | 0.965 ± 0.081 | 0.928 ± 0.082 | 0.885 ± 0.108 |
| 36 | 1.014 ± 0.079 | 0.958 ± 0.076 | 0.929 ± 0.123 |
| 38 | 1.053 ± 0.077 | 1.077 ± 0.074 | 1.032 ± 0.132 |
| 40 | 1.013 ± 0.057 | 1.059 ± 0.075 | 1.021 ± 0.097 |
| 42 | 1.046 ± 0.059 | 0.944 ± 0.050 | 0.955 ± 0.074 |
| 44 | 1.074 ± 0.072 | 0.928 ± 0.066 | 1.017 ± 0.061 |
| 46 | 1.077 ± 0.092 | 0.995 ± 0.098 | 1.124 ± 0.106 |
| 48 | 1.068 ± 0.103 | 1.177 ± 0.129 | 1.162 ± 0.135 |
| 50 | 1.116 ± 0.084 | 1.342 ± 0.165 | 1.329 ± 0.165 |
| 52 | 1.222 ± 0.146 | 1.211 ± 0.166 | 1.347 ± 0.200 |
| 54 | 1.245 ± 0.168 | 1.027 ± 0.161 | 1.190 ± 0.176 |
| 56 | 1.208 ± 0.154 | 1.101 ± 0.144 | 1.122 ± 0.128 |
| 58 | 1.291 ± 0.166 | 1.323 ± 0.137 | 1.249 ± 0.132 |
| 60 | 1.347 ± 0.258 | 1.187 ± 0.156 | 1.347 ± 0.242 |

## Punkt 1: mittleres Verhältnis über 16–22 µm

| z_h µm | Lauf 1 (Husimi) | jetzt (Wigner) |
| --- | --- | --- |
| 2,5 | 1,10 | 0,979 ± 0,015 |
| 5 | 1,20 | 1,020 ± 0,028 |
| 10 | 1,25 | 0,992 ± 0,028 |

## Punkt 2: kohärente Leistung an der Übergabe gegen exp(−μs z_h), μs = 190,31 /mm

| z_h µm | `power_coherent` (Mittel ± Fehler) | exp(−μs z_h) | Verhältnis |
| --- | --- | --- | --- |
| 0 | 1,0000 | 1,0000 | 1,000 |
| 2,5 | 0,6430 ± 0,0293 | 0,6214 | 1,035 |
| 5 | 0,4405 ± 0,0662 | 0,3861 | 1,141 |
| 10 | 0,1971 ± 0,0624 | 0,1491 | 1,322 |

Das Verhältnis steigt mit z_h über die erwarteten wenigen Prozent hinaus (14 % und 32 %), bei
Fehlern von 15 % und 32 %, also bei 5 und 10 µm nur etwa 1σ vom Erwartungswert. Die Seed-Streuung
ist wegen weniger Speckles in den Realisierungen groß.

## Punkt 3: Satz C (σ_H = 0,8 µm) gegen Satz A (σ_H = 0,4 µm)

Mittleres Verhältnis C/A über 16–22 µm: 1,007 ± 0,036 (z_h = 5 µm) und 1,007 ± 0,038
(z_h = 10 µm). Je Tiefe (16…60 µm, gleiche Reihenfolge wie oben):

- z_h = 5 µm: 0.994 1.010 1.007 1.017 1.017 0.978 0.979 0.988 1.013 1.050 1.119 1.163 1.076 1.055 1.148 1.114 0.988 1.004 1.154 1.204 1.058 0.959 0.967
- z_h = 10 µm: 0.934 0.980 1.037 1.078 1.059 1.007 1.023 1.059 1.064 1.119 1.108 1.034 0.988 1.024 1.092 1.087 1.044 0.908 0.839 0.974 1.124 1.101 1.052

Bis 34 µm liegt C/A innerhalb etwa 12 % (16–28 µm innerhalb 8 %), tiefer schwankt es zwischen
0,84 und 1,20 ohne einheitliche Richtung; im Nahbereich (16–22 µm) hängt der inkohärente Rest
nicht erkennbar von σ_H ab.

## Punkt 4: Satz B (gemessene Streuung) gegen Satz A (vektorielles Born)

Mittleres Verhältnis B/A über z ≥ 15 µm (3 Seeds gegen 5): 0,794 ± 0,017 bei z_h = 0 und
0,744 ± 0,019 bei z_h = 5 µm. Je Tiefe:

- z_h = 0: 0.846 0.824 0.798 0.785 0.758 0.749 0.727 0.668 0.660 0.636 0.665 0.742 0.783 0.781 0.749 0.724 0.768 0.889 0.956 0.992 1.031 0.964 0.771
- z_h = 5 µm: 0.847 0.829 0.826 0.826 0.809 0.724 0.656 0.652 0.680 0.673 0.648 0.675 0.751 0.788 0.703 0.606 0.531 0.605 0.859 0.993 0.840 0.795 0.788

Mit der in der BWPM gemessenen Streuung ist das Signal im Nahbereich 15–25 % niedriger als mit
vektoriellem Born, bei gleichem μs (1,002 des Born-Werts) aber kleinerem g (Phasenfunktion
g = 0,926 gegen 0,955). Ab etwa 50 µm sind die Verhältnisse rauschdominiert.
