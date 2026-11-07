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
| T14 | CUDA-Kern (S5b) gegen den CPU-Kern je Paket (≥ 99,9 % gleiche Entscheidungen), Summen auf 10⁻³; ohne Protokoll/Voxelwertung gleich (S5c) | CPU-Kern | S5 | ✅ S5b und S5c auf dem Server (2026-11-07) |
| T15 | FP32-Option (S5d): single gegen double statistisch (4σ), Energie und Paketzahl in single, GPU gegen CPU je Paket in single; double bitgleich wie T13 | double (Python-Referenz), CPU single | S5 | ✅ auf dem Server (2026-11-07) |
| T7a | Polarisation, Streuschritt: Stokes aus Feldern (Mie über S1, S2; Rayleigh über das Dipolfeld) | unabhängige Feldrechnung | S6a | ✅ |
| T7b | Polarisation im Transport: unpolarisiert = skalar; reine Zustände bleiben rein; Polarisationsgedächtnis (kleine/große Kugeln, linear/zirkular) | skalarer Transport; exakt; MacKintosh et al. (1989) | S6a | ✅ |
| T7c | Fresnel-Mueller an Brechzahlsprüngen: Fresnel-Rhombus, Brewster, nicht streuende Schicht gegen inkohärente Summe der Wege, Energie | Optik-Lehrbuch; geschlossene Summe | S6b | ✅ |
| T7e | Polarisation im C++-Kern je Paket gleich der Python-Referenz (angepasst; Brechzahlsprünge, schräg, elliptisch; Rayleigh-Halbraum) | Python-Referenz | S6c | ✅ |
| T7d | Benchmarks von Ramella-Roman, Prahl, Jacques (2005) | Literatur (Volltext nötig) | S6 | ○ |
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

## T14: CUDA-Kern auf dem Server (S5b, 2026-11-07)

RTX 6000 Ada, Treiber 590.48, nvcc 12.9 (HPC SDK 25.5), 32 CPU-Kerne, Python 3.13. Build: gcc 14
als nvcc-Host-Compiler, gcc-toolset-15 für C++ und Linken, statische libstdc++. `core_cuda.cu`
übersetzte ohne Änderung; alle Tests grün (50 + 11 mit CUDA). R in CPU und CUDA identisch.

| Fall (10⁷ Pakete) | CPU, 32 Kerne | CUDA (S5b) | R |
| --- | --- | --- | --- |
| Würfel 101³, μa = 0,05, μs = 10 /mm | 26,4 s | 33,4 s | 0,220912 |
| schwach absorbierend, μa = 0,01 /mm | 171,4 s | 339,6 s | 0,603645 |

**Befund:** Die GPU war langsamer als 32 CPU-Kerne. Ursachen: (1) jede Wechselwirkung addierte
atomar auf die Absorption eines einzigen Materials – alle Threads auf dieselbe Adresse; (2) Pakete
leben sehr verschieden lang, Warps warten auf das langsamste; (3) alles in FP64, das die
RTX 6000 Ada nur mit 1/64 des FP32-Durchsatzes rechnet. S5c behebt (1) und (2) ohne die
Physik zu ändern; (3) wäre ein FP32-Modus und ändert die Gleichheit mit der CPU (Entscheidung D7).

**S5c auf dem Server** (Commit 4deed63, 10⁷ Pakete; Tests grün, R zwischen CPU und CUDA auf 10⁻¹³):

| Fall | CPU, 32 Kerne | CUDA S5b | CUDA S5c, Protokoll | CUDA S5c, nur Summen |
| --- | --- | --- | --- | --- |
| Würfel 101³ | 25,4–26,4 s | 33,4 s | 17,3 s | 13,3 s |
| schwach absorbierend | 169,6–170,5 s | 339,6 s | 92,0 s | 62,9 s |

Lokale Akkumulation und persistente Threads machen die GPU 2,5- bzw. 5,4-mal schneller als in S5b
und 1,9- bzw. 2,7-mal schneller als 32 CPU-Kerne. Der verbleibende Engpass ist FP64 (1/64 des
FP32-Durchsatzes der RTX 6000 Ada) → FP32-Option (D7).

## T15: FP32-Option (S5d, 2026-11-07)

`fast.simulate_voxels(…, precision="single")` für CPU und CUDA (`"double"` bleibt Standard und
Referenz, Entscheidung D7). In float: Richtung, Schrittlänge, Abstände zu den Voxelflächen,
Position, Phasenfunktion (Tabellen als float-Kopie), Fresnel, Streurichtung, `log`/`sin`/`cos`/
`pow`/`cbrt`/`sqrt` (float-Funktionen); ohne Fast-Math und ohne `__logf`-Intrinsics. Immer double:
Gewicht, absorbierte Energie, Summen, Protokollwerte und die Zufallszahl (53 Bit, erst bei der
Verwendung nach float gewandelt – der Zufallsstrom je Paket bleibt derselbe), außerdem der
Absorptionsanteil je Wechselwirkung (μa/(μa+μs) in double). Nur für float übersetzt: stabile Formen
für die tabellierte Phasenfunktion und die Drehung nahe der Polachse, Renormierung der Richtung,
Paket beim Flächenübertritt auf die Fläche gesetzt (Voxelindex und Position stimmen immer),
Zufallszahl nie genau 1.

**double unverändert:** Protokoll (Seite, Orte, Richtungen, Gewichte, Wege, Streuzahlen) beider
Kerne für beide Testmedien (10⁵ Pakete) bitgleich zum Stand vor dem Umbau; die atomaren
CUDA-Summen (Absorption) weichen nur durch die Summationsreihenfolge um 10⁻¹⁵ ab.

Tests (`tests/test_s5d_single.py`, 10 Tests grün; 61 mit der ganzen Suite, CUDA erforderlich):

- single gegen double (CPU): 20 Chargen à 10⁴ Pakete je Genauigkeit (andere Seeds), R_d, T_d,
  seitlicher Verlust, Absorption je Material und R(ρ) auf 8 Ringen 0–4 mm; größte Abweichung in
  σ des kombinierten Standardfehlers: Einschluss 2,14 σ, Schichten 1,23 σ, schwach absorbierend
  2,61 σ (Grenze 4 σ).
- Energie in single (2·10⁵ Pakete): |r_spec + ΣAustritte + Absorption − 1| < 10⁻⁴; Austritte +
  Roulette-Abbrüche = Pakete (Zähler `n_exit`, `n_roulette`), Protokoll = Austritte.
- GPU single gegen CPU single je Paket (Variante A, im Code): 100 % gleiche Seite und Streuzahl
  in allen drei Medien, Summen gleich (Abweichung 0); mit Variante B 99,985 % bei der schwach
  absorbierenden Schicht, Summen auf 1,4·10⁻⁵.

**Positionsvariante:** (A) Position in float, (B) Position in double bei float-Richtung. Beide
bestehen T15; A ist auf der GPU 2,3- bis 2,6-mal schneller → A im Code, B entfernt.

**Laufzeiten** (10⁷ Pakete, `record=False`, `score_voxels=False`, RTX 6000 Ada, 32 CPU-Kerne):

| Fall | CPU double | CPU single A | CPU single B | CUDA double | CUDA single A | CUDA single B |
| --- | --- | --- | --- | --- | --- | --- |
| Würfel 101³ | 25,7 s | 23,5 s | 24,1 s | 13,3 s | **1,2 s** | 2,8 s |
| schwach absorbierend | 172,7 s | 143,4 s | 147,6 s | 61,6 s | **4,8 s** | 12,7 s |

R (Würfel) = 0,220912 ± 0,000114 in double, 0,220911 in single; R (schwach absorbierend) =
0,603645 ± 0,000112 in double, 0,603643–0,603644 in single. `--fmad=true` nur für float
(Messvariante, nicht im Code): 1,3–1,5 s bzw. 4,7–4,8 s, keine Wirkung.
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

## T7a, T7b: Polarisation (S6a)

`scattuccino.polarisation`: Stokes-Vektoren im mitgeführten Bezugssystem (u, v, w = u × v),
Mueller-Matrix der Kugel aus S1, S2 (Bohren & Huffman), θ aus S11, φ durch Verwerfen aus
S11 + S12 (Q cos 2φ + U sin 2φ); Schichten mit angepasster Brechzahl (Fresnel-Mueller folgt in
S6b). Python-Referenz; der C++-Kern folgt.

- **Streuschritt gegen Felder** (200 zufällige Richtungen, Bezugssysteme, Polarisationen und
  Winkel): Mie über die Amplituden S1, S2 und Rayleigh über das Dipolfeld −û′ × (û′ × E) auf 10⁻⁷
  (Interpolation der Tabelle) – Drehung des Bezugssystems, Mueller-Matrix und Vorzeichen von V
  stimmen.
- **Unpolarisiert** (Mie, 1 µm Polystyrol in Wasser): R 0,2006 / T 0,6373 gegen den skalaren
  Transport 0,2018 / 0,6363 (unabhängige Läufe).
- **Reine Zustände bleiben rein** auf 10⁻¹⁴ (nach Wiederherstellen von m₁₂² + m₃₃² + m₃₄² = 1
  in der interpolierten Tabelle; vorher 2,6·10⁻⁴).
- **Polarisationsgedächtnis** in Rückstreuung aus dem Halbraum:

  | Kugeln | g | linear: DOP_lin | zirkular: DOP_zirk |
  | --- | --- | --- | --- |
  | x = 0,5 | 0,04 | 0,31 | −0,20 (Helizität kehrt um) |
  | x = 10 | 0,93 | 0,03 | +0,33 (Helizität bleibt) |

  wie bei MacKintosh et al. (1989): kleine Streuer erhalten lineare, große zirkulare Polarisation.

## T7c: Fresnel an Brechzahlsprüngen, polarisiert (S6b)

`scattuccino.polarisation.interface`: Bezugssystem in die Einfallsebene drehen (p̂ = ŝ × û),
reflektieren mit der Wahrscheinlichkeit R(S) = ½(|r_p|² + |r_s|²) + ½(|r_p|² − |r_s|²) Q/I,
Stokes mit den komplexen Fresnel-Amplituden (auch Totalreflexion mit Phase).

- **Fresnel-Rhombus:** zwei Totalreflexionen bei 54,6° in Glas (n = 1,51) machen aus linearem 45°-Licht
  zirkulares: V = 0,9999999.
- **Brewster:** spiegelnde Reflexion von unpolarisiertem Licht bei θ_B (n = 1,33) ist rein
  s-polarisiert: S = (1, −1, 0, 0) auf 10⁻¹².
- **Nicht streuende Schicht** (n = 1,5, 50°, μa = 0,3/mm, d = 1 mm, linear 45°): Monte Carlo gegen
  die inkohärente Summe über alle Wege (jeder Weg behält seine p-s-Kohärenz), 2·10⁵ Pakete:
  T = (0,6299; 0,0705; 0,6240; 0) gegen (0,6303; 0,0706; 0,6245; 0), R_d = (0,0232; −0,0215;
  −0,0085; 0) gegen (0,0229; −0,0213; −0,0084; 0), alle Komponenten innerhalb 4 σ.
- **Energie** mit Streuung, zwei Brechzahlsprüngen und schrägem Einfall: 1 auf 10⁻¹²; kein Paket
  über dem Polarisationsgrad 1.

## T7e: Polarisation im C++-Kern (S6c)

`cpp/include/scattuccino/photon_pol.hpp` spiegelt `scattuccino.polarisation` Zeile für Zeile
(Streuschritt mit Verwerfen, Drehung des Bezugssystems, Mueller-Verhältnisse wie `numpy.interp`
samt Wiederherstellen der Reinheit, Fresnel-Amplituden mit Totalreflexion, polarisierter Eintritt)
und läuft im Voxelgang; die Python-Referenz nutzt dafür jetzt ebenfalls die zählerbasierten
Zufallszahlen je Paket. Fassade `scattuccino.fast.simulate_polarised_slab`.

Drei Fälle, je 4 000 Pakete (angepasst/normal/linear; Brechzahlsprünge 1,4 und 1,5, 34°,
elliptisch; Rayleigh-Halbraum mit Brechzahlsprung, zirkular): Austrittsseite, Zahl der Streuungen
und Gewicht **bitgleich**; Stokes-Vektoren, Richtungen und Bezugsvektoren gleich bis 2·10⁻¹⁴
(komplexe Division in libstdc++ gegen NumPy), im angepassten Fall bitgleich. Ein Kern: 5–7-mal
schneller als die vektorisierte Python-Referenz, mit OpenMP über alle Kerne.
