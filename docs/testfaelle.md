# Testfälle

Status: ✅ umgesetzt und grün, ◐ teilweise, ○ geplant.

| Nr. | Prüfung | Referenz | Stufe | Status |
| --- | --- | --- | --- | --- |
| T1a | Schicht ohne Streuung mit Brechzahlsprung | Beer-Lambert mit inkohärenten Mehrfachreflexionen | S1 | ✅ |
| T1b | Halbunendlich, isotrope Streuung, gleiche Brechzahl (a = 0,5; 0,9; 0,99) | Chandrasekhars H-Funktion: R = 1 − H(1)√(1 − a) | S1 | ✅ |
| T1c | Schichten mit HG-Streuung und Brechzahlsprüngen | Adding-Doubling (eigene Umsetzung) | S1b | ○ |
| T2 | Energieerhaltung (Reflexion + Transmission + Absorption = 1), zwei Schichten, Brechzahlsprünge, schräger Einfall | exakt im Erwartungswert | S1 | ✅ |
| T3 | Ortsaufgelöste Reflexion, Zeit- und Frequenzbereich | Diffusionslösungen (Farrell, Patterson, Kienle-Patterson) im diffusen Bereich | S2 | ○ |
| T4 | Subdiffusive Reflexion | Literatur (Bevilacqua/Depeursinge, Kanick/Amelink) | S3 | ○ |
| T5 | Phasenfunktion aus Mie | eigene Mie-Reihe, CliffordBem | S3 | ○ |
| T6 | Voxel/Netz = Schichten | T1 | S4 | ○ |
| T7 | Polarisation | Ramella-Roman, Prahl, Jacques (2005) | S6 | ○ |
| T8 | Fluoreszenz | Fluoreszenz-Diffusion | S7 | ○ |
| T9 | DCS/LSCI | Korrelationsdiffusion | S8 | ○ |
| T10 | OCT | Einfachstreuung, erweitertes Huygens-Fresnel | S9 | ○ |
| T11 | Jacobi-Matrizen | finite Differenzen | S10 | ○ |
| T12 | Phantome | Mikrosphären (Mie exakt), Lipidemulsionen (Literatur) | S11 | ○ |

## T1a, T1b, T2: Monte-Carlo-Referenz in Schichten (S1)

`scattuccino.layers.run` (NumPy, vektorisiert über Pakete), Statistik mit Standardfehler.

| Fall | Monte Carlo | Referenz | Abweichung |
| --- | --- | --- | --- |
| ohne Streuung, μa = 0,5/mm, d = 1 mm, n = 1,5: R | 0,05309 ± 0,00025 | 0,05357 | 1,9 σ |
| … T | 0,55962 | 0,55931 | (2·10⁵ Pakete) |
| halbunendlich isotrop, a = 0,5 | 0,11427 ± 0,00057 | 0,11523 | 1,7 σ |
| a = 0,9 (10⁶ Pakete) | 0,41488 ± 0,00035 | 0,41495 | 0,2 σ |
| a = 0,99 | 0,75224 ± 0,00102 | 0,75272 | 0,5 σ |

Energiebilanz: |Σ − 1| < 10⁻⁶ in allen Läufen (Russisches Roulette erhält die Energie nur im
Erwartungswert). H-Funktion: Iteration mit 100–400 Gauß-Punkten gleich auf 10⁻¹⁴.
