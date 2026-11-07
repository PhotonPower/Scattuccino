# Testfälle

Status: ✅ umgesetzt und grün, ◐ teilweise, ○ geplant.

| Nr. | Prüfung | Referenz | Stufe | Status |
| --- | --- | --- | --- | --- |
| T1a | Schicht ohne Streuung mit Brechzahlsprung | Beer-Lambert mit inkohärenten Mehrfachreflexionen | S1 | ✅ |
| T1b | Halbunendlich, isotrope Streuung, gleiche Brechzahl (a = 0,5; 0,9; 0,99) | Chandrasekhars H-Funktion: R = 1 − H(1)√(1 − a) | S1 | ✅ |
| T1c | Schichten mit HG-Streuung und Brechzahlsprüngen, schräger Einfall | Adding-Doubling (eigene Umsetzung); dieses gegen H-Funktion, ohne Streuung exakt und gegen van de Hulsts Werte | S1 | ✅ |
| T2 | Energieerhaltung (Reflexion + Transmission + Absorption = 1), zwei Schichten, Brechzahlsprünge, schräger Einfall | exakt im Erwartungswert | S1 | ✅ |
| T3a | Protokoll und Detektoren: Buchhaltung (Volumengitter = Absorption, ∫R(ρ) = R_d, f = 0 = CW, ∫TPSF = CW, Σ Teilwege = Weg), Gesamtreflexion beider Absorptionsverfahren gegen Adding-Doubling | exakt bzw. Adding-Doubling | S2 | ✅ |
| T3b | Mikroskopisches Beer-Lambert: Umgewichtung auf andere μa je Schicht exakt; ⟨L_i⟩ = −∂ ln R/∂μa,i | exakt | S2 | ✅ |
| T3c | R(ρ) in 4–8 Transportlängen (n = 1; 1,4) | Kienle-Patterson (Diffusion, Näherung) | S2 | ✅ |
| T3d | lateral integrierte Fluenz in 3–10 mm Tiefe; mittlere Laufzeit 4–7 mm | Diffusion (Dipol mit extrapolierter Grenze, Patterson) | S2 | ✅ |
| T4 | Subdiffusive Reflexion: gleiches g, verschiedenes γ | Verhalten nach Bevilacqua/Depeursinge (nahe der Quelle stark, weit weg schwach abhängig) | S3 | ✅ |
| T5 | Phasenfunktionen (HG, MHG, Gegenbauer, tabelliert, Mie): Normierung, Momente der Stichproben, geschlossene Formen; Mie-Grenzfälle; Monte Carlo gegen Adding-Doubling je Phasenfunktion | analytisch, Rayleigh, Extinktionsparadoxon, Adding-Doubling | S3 | ✅ (Querprobe CliffordBem offen) |
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

## T3: Protokoll, Detektoren und Diffusion (S2)

`scattuccino.transport.simulate` (Protokoll je austretendem Paket, Volumengitter (r, z), Verfahren
„implicit“ und „mbl“), `scattuccino.detectors` (R(ρ), TPSF, Frequenzbereich, Teilweglängen),
`scattuccino.diffusion` (Referenzen). Halbraum μa = 0,05/mm, μs = 10/mm, g = 0,9 (μs′ = 1/mm).

- **Exakt:** Volumengitter = Absorption (10⁻¹¹), ∫R(ρ)dA = R_d (10⁻¹²), f = 0 = CW, ∫TPSF = CW;
  Umgewichtung „mbl“ auf andere μa je Schicht gleich einem direkten Lauf (gleicher Seed, ohne
  Roulette: identische Wege, Gewichte auf 10⁻¹⁵); ⟨L_i⟩ = −∂ ln R/∂μa,i auf 10⁻⁴.
- **R(ρ) gegen Kienle-Patterson** (n = 1 / 1,4), Verhältnis MC/Diffusion:

  | ρ (mm) | 1,25 | 2,25 | 4,25 | 6,25 | 8,25 |
  | --- | --- | --- | --- | --- | --- |
  | n = 1 | 0,83 | 0,85 | 0,93 | 0,97 | 0,97 |
  | n = 1,4 | 0,86 | 0,91 | 0,91 | 0,95 | 0,95 |

  Farrells reiner Flussansatz liegt bei n = 1,4 um bis zu 29 % daneben (Kienle-Patterson: Fluss und
  Fluenz) – deshalb ist Kienle-Patterson die Referenz.
- **Fluenz:** lateral integriert 3–10 mm Tiefe innerhalb 3,5 % der Diffusion; bei 0,1 mm (unter einer
  Transportlänge) 1,87-fach – Diffusion gilt dort nicht.
- **Zeitbereich:** mittlere Laufzeit 4–7 mm: 54,2 ± 0,3 ps gegen 48,9 ps (+11 %; Diffusion zu früh
  bei wenigen Transportlängen).

## T4, T5: Phasenfunktionen und Mie (S3)

`scattuccino.phase` (HG, modifizierte HG, Gegenbauer/Reynolds-McCormick, tabelliert mit exakter
Umkehrung der stückweise linearen Dichte, Mie), `scattuccino.mie`; `Layer(phase=…)` in Transport
und Adding-Doubling.

| Phasenfunktion | g | γ | Monte Carlo gegen Adding-Doubling (d = 1 mm, μa = 0,1, μs = 10 /mm, n = 1,4) |
| --- | --- | --- | --- |
| HG 0,9 | 0,900 | 1,90 | (S1) |
| MHG (g_HG = 0,9, β = 0,8) | 0,720 | 0,97 | innerhalb 4 σ |
| Gegenbauer (g = 0,8, α = 1) | 0,914 | 2,21 | R_d 0,21935 ± 0,00058 gegen 0,21952 (3·10⁵ Pakete) |
| Mie, 1 µm Polystyrol in Wasser, 633 nm | 0,917 | 2,16 | innerhalb 4 σ |

Mie: Rayleigh-Grenzfall Q_sca auf 1,8·10⁻⁴ (x = 0,05), Winkelverteilung ∝ 1 + cos²θ;
Extinktionsparadoxon Q_ext = 2,044 bei x = 300.

**Subdiffusiver Bereich:** HG (γ = 1,90) und MHG (γ = 0,88), beide g = 0,9, also gleiches μs′:
R(ρ) zwischen 0 und 0,1 mm unterscheidet sich um den Faktor 5,2, zwischen 6 und 8 mm um 15 % – nahe
der Quelle bestimmt γ, weit weg μs′ die Reflexion. Für Verfahren mit kleinem Faserabstand (DRS mit
einer Faser, SFDI bei hohen Ortsfrequenzen) reicht μs′ deshalb nicht.
