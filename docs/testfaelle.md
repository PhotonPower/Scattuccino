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
| T6a | Voxel: Verfeinerung und Stapelgröße ändern nichts (zählerbasierter Zufall je Paket); eine und zwei Schichten mit Brechzahlsprüngen und schrägem Einfall; offene Ränder; Einschluss | exakt; Adding-Doubling; Schichttransport; Energie | S4a | ✅ |
| T6b | Tetraedernetz: Voxelgeometrie als Tetraeder je Paket gleich dem Voxelgang; Glaskugel (gekrümmte Grenze) gegen geometrische Optik | Voxelgang (S4a); Fresnel, Snellius | S4b | ✅ |
| T13 | C++-Kern (S5a) gegen die Python-Referenz, je Paket: Zufallszahlen, Phasenfunktionen, Austrittsseite, Gewicht, Zahl der Streuungen bitgleich; Orte, Richtungen, Wege bis auf gesammelte letzte Bits (≤ 7e−13); Threads ändern nichts | Python-Referenz | S5 | ✅ |
| T14 | CUDA-Kern (S5b) gegen den CPU-Kern je Paket (≥ 99,9 % gleiche Entscheidungen), Summen auf 10⁻³ | CPU-Kern | S5 | ○ (Test auf dem Server) |
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

## T6a: Voxelgeometrie (S4a)

`scattuccino.voxels.simulate_voxels` (Voxelgang mit Fresnel an Flächen zwischen verschiedenen
Brechzahlen, seitlich periodisch oder offen) mit **zählerbasierten Zufallszahlen**
(`scattuccino.random`, SplitMix64: die k-te Zahl eines Pakets hängt nur von Seed, Paketnummer und
k ab). Flächen zwischen gleichen Brechzahlen verbrauchen keine Zufallszahl.

- **Invarianz:** dasselbe Medium als 1 Voxel und als 10 × 10 × 10 Voxel, in Stapeln von 20 000 bzw.
  7 000 Paketen: R und T bitgleich (vorher mit gemeinsamem Generator nicht – die Reihenfolge der
  Ziehungen hing von der Zahl der Flächendurchgänge ab). Das ist auch die Voraussetzung dafür, dass
  der GPU-Kern (S5) dieselben Zahlen liefert.
- **Eine Schicht** (d = 1 mm, μa = 0,1, μs = 10 /mm, g = 0,9, n = 1,4) in 4 × 4 × 10 Voxeln,
  10⁵ Pakete: R_d 0,23214 / T_d 0,46170 gegen Adding-Doubling 0,23264 / 0,46123.
- **Zwei Schichten** (n 1,37 über 1,52, Substrat 1,33, 29° Einfall): R 0,17829 / T 0,48270 gegen den
  Schichttransport 0,17722 / 0,48557 (unabhängige Läufe, je 6·10⁴ Pakete).
- **Offene Ränder und Einschluss:** Würfel 2,1 mm, absorbierende Kugel (r = 0,3 mm, μa = 2/mm,
  „Gefäß“) 0,65 mm unter dem Eintrittspunkt: R sinkt von 0,092 auf 0,057, die Kugel absorbiert 41 %;
  Energiebilanz mit seitlichen Verlusten auf 10⁻⁴.

## T13: C++-Kern (S5a)

`cpp/src/core.cpp` (C++20, OpenMP über Pakete, `-ffp-contract=off`) spiegelt
`scattuccino.voxels` Zeile für Zeile, mit denselben zählerbasierten Zufallszahlen und derselben
Reihenfolge der Ziehungen je Paket. Fassade `scattuccino.fast.simulate_voxels(…, backend=…)`.

- Zwei Medien (offener Würfel mit MHG-Einschluss; zwei Schichten mit Mie- und Gegenbauer-
  Phasenfunktion, Brechzahlsprüngen, schrägem Einfall), je 10⁴ Pakete: **jede** Entscheidung
  gleich – Austrittsseite, Gewicht und Zahl der Streuungen bitgleich; Orte, Richtungen und Wege
  weichen nur in den letzten Bits ab (cbrt/pow/log von libm gegen NumPy, gesammelt bis 7·10⁻¹³).
- Threads (1 gegen 4) ändern nichts.
- Laufzeit (ein Kern): Würfel 21³ mit Einschluss, 2·10⁴ Pakete: Python 0,52 s, C++ 0,06 s
  (Faktor 9); schwach absorbierende Schicht (μa = 0,01/mm, μs = 10/mm, 20 mm), 2·10⁵ Pakete: 13 s
  (NumPy-Schichtcode: 80 s); mit OpenMP auf allen Kernen entsprechend schneller.

## T6b: Tetraedernetze (S4b)

Transport im C++-Kern (`cpp/include/scattuccino/photon_tet.hpp`, für CPU und später GPU),
Nachbarn aus gemeinsamen Flächen, Fresnel und Brechung um die Flächennormale, Außenbrechzahl je
Randfläche; Eintritt an einem Punkt einer Randfläche (beliebig geneigt). Abweichung vom Plan: ohne
eigene Python-Fassung, weil die Prüfung gegen den Python-geprüften Voxelgang exakt ist.

- **Voxel als Tetraeder** (11³ Voxel, je 6 Tetraeder nach Kuhn, 7 986 Tetraeder; Gewebe mit
  kugelförmigem Einschluss anderer Brechzahl und MHG-Streuung; Substrat 1,33): für alle 2·10⁴ Pakete
  dieselben Gewichte und Streuzahlen wie der Voxelgang, Orte und Richtungen auf 10⁻¹⁴, Absorption
  je Material gleich.
- **Glaskugel** (n = 1,5, Ikosaeder Stufe 5, 20 480 Tetraeder), Strahl unter 30,4° zur Facette:
  spiegelnde Reflexion = Fresnel (exakt), Anteil in der Richtung des direkt durchgehenden Strahls
  (Ablenkung 2(θᵢ − θₜ)) 0,9188 gegen (1 − R₁)(1 − R₂) = 0,9185; Energie auf 10⁻¹².
