# Testfälle

Status: ✅ umgesetzt und grün, ◐ teilweise, ○ geplant.

| Nr. | Prüfung | Referenz | Stufe | Status |
| --- | --- | --- | --- | --- |
| T1a | Schicht ohne Streuung mit Brechzahlsprung | Beer-Lambert mit inkohärenten Mehrfachreflexionen | S1 | ✅ |
| T1b | Halbunendlich, isotrope Streuung, gleiche Brechzahl (a = 0,5; 0,9; 0,99) | Chandrasekhars H-Funktion: R = 1 − H(1)√(1 − a) | S1 | ✅ |
| T1c | Schichten mit HG-Streuung und Brechzahlsprüngen, schräger Einfall | Adding-Doubling (eigene Umsetzung); dieses gegen H-Funktion, ohne Streuung exakt und gegen van de Hulsts Werte | S1 | ✅ |
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

## T1c: Adding-Doubling und Monte Carlo (S1)

`scattuccino.adding_doubling.slab`: Matrixoperator-Form mit dünner Startschicht (erste Ordnung,
Δτ = τ/2²⁴), Gauß-Quadratur geteilt am Grenzwinkel der Totalreflexion, kollimierter Strahl
getrennt, Fresnel an beiden Grenzen; ≈ 0,1 s je Fall.

Prüfung des Adding-Doubling selbst: gegen die H-Funktion 5·10⁻⁷ / 7·10⁻⁶ / 5·10⁻⁵ (a = 0,5 / 0,9 /
0,99), ohne Streuung exakt (10⁻¹²), Literaturfall μa = 1, μs = 9 /mm, g = 0,75, d = 0,2 mm, n = 1:
R_d = 0,097395, T_d = 0,660958 gegen van de Hulsts 0,09734 / 0,66096 (zitiert in wang1995).

Monte Carlo gegen Adding-Doubling (2·10⁵ Pakete):

| d (mm) | μa, μs (1/mm) | g | n (oben/Schicht/unten) | cos θ | R_d MC / AD | T_d MC / AD |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 0,1 / 10 | 0,9 | 1 / 1,4 / 1 | 1 | 0,23386 ± 0,00073 / 0,23264 (1,7 σ) | 0,45976 ± 0,00086 / 0,46123 (−1,7 σ) |
| ∞ | 0,1 / 0,9 | 0,75 | 1 / 1,33 / – | 1 | 0,09289 ± 0,00041 / 0,09246 (1,0 σ) | – |
| 0,5 | 0,3 / 6 | 0,7 | 1 / 1,45 / 1,52 | 0,7 | 0,11518 ± 0,00057 / 0,11495 (0,4 σ) | 0,56220 ± 0,00072 / 0,56293 (−1,0 σ) |
| 0,2 | 0,05 / 20 | 0,95 | 1,33 / 1,37 / 1 | 1 | 0,17234 ± 0,00082 / 0,17158 (0,9 σ) | 0,81078 ± 0,00085 / 0,81159 (−1,0 σ) |
