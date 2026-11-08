# Bericht S9b-2: Kalibrierläufe des Zufallsmediums in BWPM

Auftrag: `docs/agent-server-s9b2.md`. Skript: `tools/calibrate_bwpm.py`, Rohdaten
`tools/out/calibration_runs.jsonl`, Winkelverteilung `tools/out/phase_bwpm_vs_born.csv`
(Gruppe A, μs = 20 /mm, Seed 1, Tiefe 10 µm).

## Umgebung

NVIDIA RTX 6000 Ada (48 GB), nvcc 12.9, BWPM mit CUDA-Kern (`fast.run_3d(..., iterations=0,
backend="cuda")`), FP64. n₀ = 1,38, λ = 0,8 µm, Gitter isotrop 0,116 µm, ebene Welle,
Polarisation x. Tiefe 4/μs (Born), 20 gespeicherte Ebenen, Anpassung von μs auf dem Bereich
0,9 > |⟨E⟩|² > 0,05.

## Abweichungen vom Auftrag (bitte prüfen)

1. **Lateral 512² statt 1024²**, 16 statt 64 Quantisierungsstufen. Die Kosten der BWPM wachsen mit
   der Zahl verschiedener Indexpaare benachbarter Ebenen (je Ebene eine Matrixanwendung plus
   inverse FFT je Paar), also grob mit Stufen² × n². 1024² mit 64 Stufen wären für die 15 Läufe
   etwa 100 h; 1024² mit 16 Stufen etwa 13 h; 512² mit 16 Stufen 1,5 h gesamt (Tabelle).
2. **Prüfung der Stufenzahl** (256², `tools/out/calibration_validation_levels.jsonl`): 16, 32 und
   64 Stufen geben g auf 10⁻⁴ gleich, das angepasste μs weicht um 0,7 % ab. Quantisierungsfehler
   (rms) 0,003–0,006 gegen σ_n ≥ 0,009.
3. **Prüfung der lateralen Größe**: μs = 50 /mm, Seed 1, 1024² gegen 512²: Verhältnis 0,978 gegen
   0,959, g gleich (0,9623 / 0,9614). Der Unterschied liegt innerhalb der Streuung zwischen
   Realisierungen (unten).
4. Das Skript berechnet g mit und ohne Raumwinkel-Jacobi (siehe Befund 3); die Bibliotheksfunktion
   `calibrate_with_bwpm` wurde nicht verändert.

## Tabelle (512², 16 Stufen, Tiefe bei Born-μs·z = 4)

g (Bib.) = Auswertung wie `calibrate_with_bwpm` (erste Tiefe), g (Ω) = mit Raumwinkel-Jacobi
dΩ = dk⊥²/(k² cos θ) bei ≈ 0,25/μs. Leistung: max |total − 1| über alle Ebenen.

| σ_n | l_c µm | D | Seed | Born μs | fit μs | Verh. | Born g | g (Bib.) | g (Ω) | max\|tot−1\| | Laufzeit s | GPU MiB |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0,0128 | 0,5 | 4 | 1 | 10 | 9,88 | 0,988 | 0,9462 | 0,9634 | 0,9578 | 0,004 | 1015 | 2208 |
| 0,0181 | 0,5 | 4 | 1 | 20 | 19,56 | 0,978 | 0,9462 | 0,9629 | 0,9574 | 0,006 | 521 | 1346 |
| 0,0286 | 0,5 | 4 | 1 | 50 | 47,96 | 0,959 | 0,9462 | 0,9614 | 0,9554 | 0,012 | 202 | 828 |
| 0,0292 | 0,2 | 4 | 1 | 20 | 18,05 | **0,902** | 0,8307 | 0,9056 | 0,8885 | 0,24 | 552 | 1634 |
| 0,0127 | 1,0 | 4 | 1 | 20 | 19,44 | 0,972 | 0,9809 | 0,9847 | 0,9828 | 0,002 | 463 | 1174 |
| 0,0090 | 2,0 | 4 | 1 | 20 | 17,84 | **0,892** | 0,9938 | 0,9940 | 0,9934 | 0,002 | 451 | 976 |
| 0,0236 | 0,5 | 3,5 | 1 | 20 | 19,22 | 0,961 | 0,9217 | 0,9523 | 0,9438 | 0,066 | 538 | 1634 |
| 0,0144 | 0,5 | 5 | 1 | 20 | 19,87 | 0,993 | 0,9725 | 0,9766 | 0,9743 | 0,007 | 465 | 1058 |
| 0,0128 | 0,5 | 4 | 2 | 10 | 9,32 | 0,932 | 0,9462 | 0,9641 | 0,9586 | 0,004 | 1003 | 2208 |
| 0,0181 | 0,5 | 4 | 2 | 20 | 19,23 | 0,962 | 0,9462 | 0,9635 | 0,9580 | 0,007 | 502 | 1346 |
| 0,0286 | 0,5 | 4 | 2 | 50 | 48,53 | 0,971 | 0,9462 | 0,9616 | 0,9558 | 0,013 | 201 | 828 |
| 0,0128 | 0,5 | 4 | 3 | 10 | 9,12 | 0,912 | 0,9462 | 0,9613 | 0,9555 | 0,002 | 1029 | 2208 |
| 0,0181 | 0,5 | 4 | 3 | 20 | 18,81 | 0,941 | 0,9462 | 0,9623 | 0,9566 | 0,004 | 519 | 1634 |
| 0,0286 | 0,5 | 4 | 3 | 50 | 49,41 | 0,988 | 0,9462 | 0,9617 | 0,9559 | 0,010 | 210 | 944 |
| 0,0286 | 0,5 | 4 | 1 | 50 | 48,91 | 0,978 | 0,9462 | 0,9623 | 0,9564 | 0,012 | 925 | 2014 |

Die letzte Zeile ist der 1024²-Lauf. Host-Speicher (Feld in FP64, 3450 × 512²): bis 85 GB bei
μs = 10 /mm; die Felderzeugung (FFT) ist der Speicherengpass, nicht die GPU.

## Befund

**Wo weicht μs um mehr als 10 % von Born ab?** Fast nirgends systematisch: im Mittel liegt das
angepasste μs 1–6 % unter Born (Gruppe A: Mittel über Seeds 0,970 / 0,960 / 0,944 für 50 / 20 /
10 /mm). Die Streuung zwischen Realisierungen beträgt ±3 % (Spanne 0,912–0,988), also größer als
der Unterschied zwischen den Gruppen. Über 10 % liegen die einzelnen Läufe l_c = 2,0 µm (0,892:
vermutlich wenige Korrelationszellen in 512², nicht untersucht) und l_c = 0,2 µm (0,902; hier sind außerdem die Leistungserhaltung (24 %) und
D = 3,5 (7 %) verletzt, vermutlich weil Streuung zu großen Winkeln und Rückstreuung in der
Vorwärts-FPWPM fehlen (nicht untersucht)). Mögliche Ursache des Defizits: die Feldvarianz erreicht nur
92 % (l_c = 0,5, D = 4) bzw. 72–81 % (D = 3,5; l_c = 0,2) von σ_n², weil die hohen Ortsfrequenzen
oberhalb der Gitterfrequenz fehlen (Zusammenhang nicht geprüft).

**Wo weicht g um mehr als 10 % von Born ab?** Das vorwärts gemessene g selbst liegt für alle
Läufe innerhalb 1–3,5 % von Born, aber das ist wegen g → 1 unempfindlich. Auf 1 − g (der Größe,
die in Transportgleichung und μs′ eingeht) gemessen (Raumwinkel-Jacobi berücksichtigt):

| Gruppe | 1 − g Born | 1 − g BWPM | Abweichung |
| --- | --- | --- | --- |
| l_c = 0,5, D = 4 | 0,0538 | 0,043 | −21 % |
| l_c = 0,2, D = 4 | 0,169 | 0,112 | −34 % |
| l_c = 1,0, D = 4 | 0,0191 | 0,0172 | −10 % |
| l_c = 2,0, D = 4 | 0,0062 | 0,0066 | +6 % |
| l_c = 0,5, D = 3,5 | 0,0783 | 0,0562 | −28 % |
| l_c = 0,5, D = 5 | 0,0275 | 0,0257 | −7 % |

Das BWPM streut also weniger zu großen Winkeln als Born (Verhältnis der Dichten bei 40°:
etwa 0,8, `phase_bwpm_vs_born.csv`); der Effekt wächst mit kleinem l_c und kleinem D. Für l_c ≥ 1 µm und D ≥ 5 liegt g innerhalb 10 % von Born. Die Streuung von g (Ω) zwischen Seeds beträgt 0,001.

**Weitere Beobachtungen.**

1. g hängt von der Auswertungstiefe ab: g (Ω) fällt für l_c = 0,5 µm von 0,957 bei 0,25/μs auf
   0,902 bei 4/μs durch Mehrfachstreuung. Die Tabelle gibt das erste Tiefenniveau; die Serie über
   Tiefe steht im JSONL (`g_series`, `g_series_jacobian`).
2. Die Bibliotheksauswertung des g (`calibrate_with_bwpm`) lässt den Raumwinkel-Jacobi
   dΩ = dk⊥²/(k² cos θ) weg und überschätzt g damit um 0,005 (D = 4) bis 0,03 (D = 3,5). Das ist
   ein Befund an der Bibliothek, keine Änderung (Auftrag verbietet Codeänderungen): Issue mit
   Label `frage`.
3. Leistungserhaltung: max |total − 1| zwischen 0,2 % (l_c ≥ 1 µm) und 1,3 % (μs = 50 /mm) für
   D ≥ 4, l_c ≥ 0,5; 6,6 % (D = 3,5) und 24 % (l_c = 0,2 µm), vermutlich wegen Streuung außerhalb der
   Vorwärtsausbreitung.

**Für S9b-3:** μs = Born ·(0,96 ± 0,03) für l_c ≥ 0,5 µm und D ≥ 4 (σ_n entsprechend um 2 %
anheben); g aus Born nur für l_c ≥ 1 µm oder D ≥ 5 auf 10 % von 1 − g, sonst ist die
Kalibrierung über (σ_n, l_c, D) an der gemessenen Winkelverteilung nötig. Offen: 1024² mit
64 Stufen laut Auftrag (Kosten oben, Rückfrage an den Maintainer, ob nötig).
