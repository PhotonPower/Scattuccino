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
| T8 | Fluoreszenz und Raman: ohne Streuung gegen geschlossene Lösung; Energie der Emission; vorwärts gegen adjungiert (Fluchtfunktion); räumlich versetzte Detektion (SORS) | E₂-Integral; exakt; unabhängige Läufe; Verhalten | S7a | ✅ |
| T8b | Fluoreszenzlebensdauer: Zeitpunkt der Anregung (opl = n z ohne Streuung), Verzögerung je Austritt exponentiell, Frequenzbereich F_τ/F_0 = 1/(1 + iωτ) (Phase, Modulation) | exakt (Faltung mit e^(−t/τ)/τ); lakowicz2006 | S7b | ✅ |
| T8c | Emissionsspektrum: ohne Streuung je Wellenlängen-Bin gegen geschlossene Form (Inner-Filter-Effekt), eine Anregung für alle Bins, Energie je Bin, gleiche Bins = eine Wellenlänge | exakt; Monte Carlo | S7c | ✅ |
| T9 | Dynamik: Impulsübertrag (⟨1 − cos θ⟩ = 1 − g), DCS gegen Korrelationsdiffusion, Speckle-Kontrast gegen geschlossene Form, Tiefenempfindlichkeit in zwei Schichten | exakt; Boas & Yodh; Bandyopadhyay | S8a | ✅ |
| T9b | Laser-Doppler: gekennzeichnete bewegte Streuer, g1 mit unverschobenem Anteil gegen Gl. 22 (P_m aus dem Protokoll, I₁ aus Gl. 2), Referenz bonner1981 gegen eigene Quadraturen, Spektrum und erster Moment | bonner1981 | S8b | ✅ |
| T10a | OCT (Intensitätsmodell): Einfach-Rückstreuung gegen exaktes Integral; Kohärenz-Gating erhält das Signal; Anteil der Mehrfachstreuung und Beer-Gesetz des OCT | exakt; Integral; Verhalten | S9a | ✅ |
| T10b | OCT mit Fokus (Fieldtracing), Übergang BWPM, erweitertes Huygens-Fresnel; S9c-1: geschlossene Form gegen Quadratur von Gl. 21 und A1, Grenzfälle, thrane2000 Abb. 2 und 3; S9c-2: allgemeines EHF mit Gaußscher MTF = geschlossene Form | Thrane, Lutomirski–Yura; Querprobe BWPM | S9b/c | ◐ (S9c-1, S9c-2 ✅) |
| T10c | Rückstreuung der BWPM (B8) an δn: Born erster Ordnung derselben Realisierung (homogene Schicht gegen Airy), Leistung im Born-Bereich, Feldfehler linear in σ_n | exakt (Airy); Born erster Ordnung | S9c-4a | ✅ |
| T17a | Hybrid-OCT, Kette: Fokus des Gauß-Strahls, Modenkopplung eines Strahls, Leistung an der Übergabe | analytisch | S9b-3a | ✅ |
| T17b | Next-Event-Schätzer der Modenkopplung: gleich dem direkten Schätzer (Mehrmoden-Akzeptanz), Ereignisse stapelweise = gespeichert, Kopplung eines Punktstreuers ∝ Modenintensität | direkter Schätzer; exakt | S9b-4 | ✅ |
| T16b | Gemessene Streuparameter aus der BWPM für den Monte Carlo: μs aus der Extinktion; Phasenfunktion vorwärts gemessen, rückwärts und unter der Winkelauflösung des Gitters Born (absolut) – normiert, Rückwärtsteil = Born auf 5 % | exakt / Born | S9b-5 | ✅ |
| T18 | Kohärente Übergabe: Wigner-Pakete mit den Momenten des Gaußstrahls; Projektion teilt die Leistung exakt; Einfachstreuung im engen Fokus (z_h = 0) gegen geschlossene Form; kohärenter Anteil hinter 3 µm Zufallsmedium gegen e^(−μs z_h) | analytisch; Born | S9b-6 | ✅ |
| T17 | Hybrid-OCT: Unabhängigkeit von z_h (Server; #29: w_f = 2 µm mit Husimi nicht bestanden, Ursache Husimi-Glättung → S9b-6; Wiederholung mit engem Fokus nach S9b-6), gegen erweitertes Huygens-Fresnel | S9b-4, S9c | ◐ |
| T16 | Zufallsmedium (Whittle-Matérn): Born gegen geschlossene Form, Feldstatistik, eine BWPM-Realisierung reproduziert Born-μs und g | exakt; Statistik; Born | S9b-2 | ✅ (Kalibrierläufe auf dem Server offen) |
| T16c | Vektorielle Born-Näherung (Dipolfaktor) gegen geschlossene Form; Leistung ebener Wellen mit E_z; BWPM-Realisierung: μs gegen vektorielles Born, g vorwärts als bekannte Abweichung (BWPM E9, `xfail`); S9b-8: δn als Born-Inkrement (BWPM B7), Leistungserhaltung als bekannte Abweichung (`xfail`) | exakt; vektorielles Born | S9b-7 | ◐ (g: Grenze der FPWPM) |
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

## T8: Fluoreszenz und Raman (S7a)

`scattuccino.fluorescence.simulate_fluorescence`: zweistufig – Anregung mit Absorptionsereignissen
(Anteil μa,f/μa geht an den Fluorophor), Emission als isotrope innere Quellen
(`transport.Sources`) mit Quantenausbeute bzw. Raman-Effizienz und den optischen Eigenschaften bei
λ_m (Re-Absorption inklusive, ohne zweite Emission).

- **Ohne Streuung** (μa,x = 0,8 davon μa,f = 0,3, μa,m = 0,4 /mm, Q = 0,2, n = 1): absorbierter Anteil
  exakt μa,f/μa,x = 0,375; austretende Fluoreszenz 0,016812 gegen
  Q ∫ μa,f e^(−μa,x z) ½E₂(μa,m z) dz = 0,016901 (1,3 σ, 2·10⁵ Pakete).
- **Energie:** emittiert = austretend + re-absorbiert auf 10⁻¹².
- **Vorwärts gegen adjungiert** (Schicht 2 mm, μs = 10/9 /mm, g = 0,9, n = 1,4): 0,006974 gegen
  Q Σ A_f(z) E_m(z) = 0,006930 (+0,6 %; Fluchtfunktion aus isotropen Punktquellen in 20 Tiefen).
- **SORS:** Raman-aktiv nur die obere Schicht (0,5 mm) oder nur die Schicht darunter: Verhältnis
  tief/oben 0,54 bei 0–0,5 mm Abstand, 3,16 bei 4–6 mm – der tiefe Anteil wächst mit dem Abstand um
  den Faktor 5,8.

## T8b: Fluoreszenzlebensdauer im Zeit- und Frequenzbereich (S7b)

`simulate_fluorescence(lifetime_ps=τ)`: Ein Emissionspaket startet zur Zeit seines
Anregungsereignisses (optischer Weg / c, neu in den Absorptionsereignissen: `events["opl"]`) plus
einer exponentiellen Verzögerung mit Lebensdauer τ aus einem eigenen Zufallsstrom;
`Sources.time_ps`, `Protocol.time_offset_ps`, `Protocol.time_ps` enthält beides. Test
`tests/test_s7b_lifetime.py`.

- Zeitpunkt der Anregung: Platte ohne Streuung mit angepasster Brechzahl außen, opl = n z für alle
  Ereignisse mit Gewicht auf 10⁻¹². Bei μs = 0 läuft ein Paket nach der vollständigen Absorption
  gewichtslos weiter, bis das Roulette es beendet; diese Ereignisse tragen das Gewicht 0.
- Lebensdauer (Halbraum μs = 10/mm, HG 0,8, n = 1,4; 4·10⁴ Emissionspakete, τ = 2 ns): Transport
  und Gewichte sind für τ = 0 und τ > 0 gleich, die Zeit jedes Austritts ist um die gezogene
  Verzögerung später; Mittel 2017 ps und Streuung 2014 ps bei 32 257 Austritten (exponentiell:
  beide τ).
- Frequenzbereich: F_τ(ω)/F_0(ω) gegen 1/(1 + iωτ) innerhalb 4 σ (σ aus der Varianz von e^(−iωd));
  gemessen tan φ = 0,303 / 1,008 / 3,065 und Modulation 0,957 / 0,705 / 0,318 bei ωτ = 0,3 / 1 / 3
  (erwartet ωτ und 0,958 / 0,707 / 0,316).

## T8c: Emissionsspektren (S7c)

`simulate_fluorescence_spectrum(stack_x, stacks_m, spectrum, …)`: Die Emission wird in
Wellenlängen-Bins mit dem normierten Spektrum f_b aufgeteilt; Bin b trägt Q f_b × Absorption des
Fluorophors und wird mit seinen optischen Eigenschaften transportiert (Re-Absorption je
Wellenlänge). Eine Anregung für alle Bins; Pakete ∝ f_b. Test `tests/test_s7c_spectrum.py`.

- **Ohne Streuung** (μa,x = 0,8 davon μa,f = 0,3 /mm, Q = 0,2; μa,m = 0,8 / 0,4 / 0,1 /mm, f_b =
  0,2 / 0,5 / 0,3; 1,2·10⁵ Emissionspakete): jedes Bin gegen Q f_b ∫ μa,f e^(−μa,x z) ½E₂(μa,m,b z) dz
  innerhalb 4 σ. Detektiertes Spektrum 0,12 / 0,44 / 0,43 statt 0,2 / 0,5 / 0,3 – der
  Inner-Filter-Effekt verschiebt es zu den schwach re-absorbierten Bins.
- **Eine Anregung:** Fluorophor-Absorption gleich dem Lauf mit einer Wellenlänge (gleicher Seed,
  bitgleich); Energie je Bin (emittiert = austretend + re-absorbiert) auf 10⁻⁴ – das Russische
  Roulette erhält die Energie nur im Erwartungswert (gemessen 10⁻⁵).
- **Gleiche Eigenschaften in allen Bins:** Summe der detektierten Fluoreszenz gleich dem Lauf mit
  einer Wellenlänge innerhalb 4 σ (unabhängige Läufe).
- `simulate_fluorescence` (S7a/S7b) ist in Anregung und Emission aufgeteilt; Seeds und Ergebnisse
  unverändert (Tests T8, T8b grün).

## T9: Dynamik – DCS und Speckle-Kontrast (S8a)

Protokoll mit Impulsübertrag je Schicht Y_t = Σ(1 − cos θ) (`Protocol.momentum`),
`scattuccino.dynamics` (g1, g2, Speckle-Kontrast; Brownsche Bewegung und ungerichtete Strömung),
`diffusion.g1_correlation_diffusion`.

- ⟨1 − cos θ⟩ aus dem Protokoll: 0,10006 (HG g = 0,9; exakt 0,1).
- **DCS** (Halbraum μa = 0,01–0,05, μs′ = 1 /mm, n = 1,4, λ = 785 nm, αD_B = 10⁻⁸ cm²/s): g1 bei
  2–3 mm folgt der Korrelationsdiffusion auf 0,02 (z. B. τ = 1 ms: 0,378 gegen 0,400); bei 1–2 mm
  (unter zwei Transportlängen) bis 0,11 daneben – Diffusion gilt dort nicht.
- **Speckle-Kontrast**: numerische Integration von g1 = e^(−τ/τ_c) gegen die geschlossene Form
  K = 0,42427 auf 10⁻⁸.
- **Zwei Schichten**, Bewegung nur unter 1 mm: g1(0,3 ms) bei 3–4 mm deutlich kleiner als bei
  0,5–1 mm – die Tiefenempfindlichkeit von DCS.

## T9b: Laser-Doppler mit gekennzeichneten bewegten Streuern (S8b)

`Layer.moving_fraction`, `Layer.phase_moving`; das Protokoll führt `momentum_moving` und
`n_moving` je Schicht. `dynamics.g1_tagged`, `unshifted_fraction`, `doppler_spectrum` (bonner1981
Gl. 26), `mean_frequency` (Gl. 32); Referenz `reference.bonner_*`. Test `tests/test_s8b_doppler.py`.

- **Referenz gegen ihre Definitionen:** I₁ nach Gl. 18 gegen die Quadratur von Gl. 2 (mit Gl. 12′,
  13′) auf 10⁻⁸, Gl. 20 als Grenzfall; S_j nach Gl. 31 gegen die Kosinus-Quadratur von Gl. 30 auf
  10⁻⁸ – **Gl. 31 druckt S₁ = e^(−W), Gl. 30 gibt ½ e^(−W)** (passend zu S₂, S₃ und Gl. A3;
  Issue mit Label `paper`); ⟨ω⟩ nach Gl. 36/A4 gegen den numerisch berechneten ersten Moment des
  Spektrums von Gl. 24 auf 2·10⁻³ (m̄ = 0,2 / 1 / 3), auch direkt aus G(τ).
- **Monte Carlo** (Halbraum μs = 10/mm, HG 0,8, μa = 0,01/mm, n = 1,4, λ = 785 nm; 2 % der Stöße an
  bewegten Zellen mit dem Strukturfaktor von Gl. 13′, a = 2,75 µm; Gaußsche Geschwindigkeiten,
  v = 1 mm/s; 4·10⁴ Pakete, Abstand 1–3 mm): Anteil der Stöße 0,02 und ⟨1 − cos θ⟩ der bewegten
  Stöße wie die Phasenfunktion der Zellen. g1 gegen Gl. 22 mit P_m aus dem Protokoll und I₁ aus
  Gl. 2 (Quadratur über die Phasenfunktion, unabhängig vom Transport) auf 10⁻³ (0,9902 / 0,8239 /
  0,3680 / 0,2916 gegen 0,9903 / 0,8246 / 0,3684 / 0,2910 bei τ = 0,2 / 1 / 5 / 20 ms).
- **Befund Poisson-Annahme (bonner1981 Gl. 23):** Die Zahl der Stöße mit bewegten Zellen ist für
  die detektierten Photonen nicht Poisson-verteilt (m̄ = 2,24, P₀ = 0,285 statt e^(−m̄) = 0,107),
  weil die Weglängen breit verteilt sind. g1 nach Gl. 23 liegt bei großen τ um den Faktor 2–2,6
  darunter (0,194 statt 0,368 bei 5 ms). Gl. 22 mit dem gemessenen P_m ist die Referenz.
- Ohne bewegte Streuer zieht der Transport keine zusätzlichen Zufallszahlen; die bisherigen Läufe
  bleiben bitgleich (Tests S1–S9 unverändert grün). `Motion.alpha` (S8a, mittlerer Faktor) bleibt
  für DCS; es kennt keinen unverschobenen Anteil.

## T10a: OCT-A-Scans aus dem Protokoll (S9a)

`scattuccino.oct`: Detektor am Eintrittsort (Radius, Akzeptanzwinkel), Tiefe z = OPL/2 mit
OPL = Σ nᵢLᵢ, Kohärenzfunktion (Gauß, FWHM l_c); mittleres, speckle-freies Signal wie im
Monte-Carlo-OCT nach Yao & Wang (1999).

- **Einfach-Rückstreuung** (Halbraum μa = 0,1, μs = 2 /mm, HG g = 0,3, n = 1, Akzeptanz 30°,
  4·10⁵ Pakete): Pakete mit genau einer Streuung gegen das exakte Integral
  ∫∫ μs e^(−μt z) p(−cos θ) e^(−μt z/cos θ) 2π sin θ dθ dz, gebinnt nach z(1 + 1/cos θ)/2 – alle
  10 Tiefenbins innerhalb 1,3 σ.
- **Kohärenz-Gating** (l_c = 10 µm) erhält das detektierte Signal auf 10⁻³.
- **Mehrfachstreuung** (Detektor 1 mm, 30°): Anteil 20 % bei 0,05 mm Tiefe, 99 % bei 1,9 mm; der
  Einfachstreu-Anteil fällt mit 2μt auf 1,5 % genau (Beer-Gesetz des OCT).

## T16: Zufallsmedium und Kalibrierung (S9b-2)

`scattuccino.microstructure`: Whittle-Matérn-Spektrum, Born-Streuung (μs, Phasenfunktion),
Gaußsches Zufallsfeld auf periodischem Gitter, Quantisierung für BWPM,
`calibrate_with_bwpm` (Extinktion des kohärenten Feldes → μs, Winkelspektrum → g, vorwärts).

- Born für D = 4 gegen die geschlossene Form μs = 8k⁴σ_n²l_c³/(n₀²(1 + 4k²l_c²)) und g: 10⁻⁶.
- Feldvarianz über 8 Realisierungen innerhalb 8 % der Gittererwartung (eine Realisierung schwankt
  um ±6–9 % durch die wenigen langwelligen Moden); Korrelation e^(−r/l_c) auf 0,06.
- **Eine BWPM-Realisierung** (96² × 60 bei 0,1 µm, λ = 0,8 µm, n₀ = 1,38, l_c = 0,5 µm,
  σ_n = 0,057 → Born μs = 200 /mm, g = 0,935): μs aus der Extinktion 192 /mm (−4 %), g vorwärts
  0,942 gegen Born 0,946, Leistung auf 0,5 % erhalten. Bei kleinem lateralem Feld bleibt ein
  kohärenter Boden von einigen Prozent (endlich viele Speckles); die Kalibrierläufe auf dem Server
  (`docs/agent-server-s9b2.md`) nutzen 1024² und realistische μs.

## T17a: Hybrid-OCT-Kette (S9b-3a)

`scattuccino.hybrid`: Gauß-Strahl mit Taille in der Tiefe z_f → BWPM durch eine Realisierung der
Mikrostruktur bis z_h → Husimi-Pakete (`rtt_fusion.phasespace`, periodische Fenster) → Monte Carlo
mit den Born-Parametern (auch auf dem Rückweg durch 0 < z < z_h) → Kopplung der austretenden
Pakete in die Detektionsmode (Rückverfolgung zur Taille, Wigner-Dichte) → Kohärenz-Gating.

- Gauß-Strahl: rms-Radius im Fokus w_f/2 auf 2 % (ein Vorzeichenfehler im komplexen
  Strahlparameter ließ den Strahl zuerst divergieren – vom Test gefunden).
- Leistung an der Übergabe: 1,000 (z_h = 0) bzw. 0,998 (z_h = 3 µm durch das Zufallsmedium).
- **Befund:** Der direkte Schätzer des modengekoppelten Signals ist zu verrauscht für Vergleiche
  zwischen z_h: Nur wenige Pakete treffen den Phasenraum der Mode (Taille 2 µm, NA ≈ 0,05), das
  Signal schwankt jenseits des Fokus um Größenordnungen. Nahe am Fokus liegen z_h = 0 und 3 µm in
  derselben Größe (10⁻³). Nötig ist ein Next-Event-Schätzer (bei jeder Streuung der Beitrag aller
  Richtungen in die Mode, mit Abschwächung bis zur Oberfläche) – S9b-4.

## T17b: Next-Event-Schätzer (S9b-4)

`hybrid.next_event_signal`: Bei jeder Streuung der erwartete Anteil, der in die Detektionsmode
gestreut wird und die Oberfläche ohne weitere Wechselwirkung erreicht – das Integral über die
Richtungen analytisch im paraxialen Grenzfall der Mode:
E = w p(u·(−ẑ)) e^(−μt z₀) (2π/k²) e^(−ρ₀²/2V)/V, V = Varianz der Modenintensität in der Tiefe z₀
(Reziprozität: Kopplung eines Punktstreuers ∝ |E_Mode(r₀)|²). Transport liefert die Streuereignisse
auf Wunsch stapelweise an eine Funktion (`record_scattering=callable`), damit große Läufe nichts
speichern.

- Gegen den direkten Schätzer mit einer Mehrmoden-Akzeptanz (10 µm, NA ≈ 0,14), bei der dieser
  konvergiert: Gesamtsignal auf 3–4 % gleich (paraxiale Näherung), im Hybridlauf 3,14 gegen
  3,00·10⁻⁶ und mittlere Tiefe 25,6 gegen 26,9 µm.
- Mit der beugungsbegrenzten Mode (2 µm) zeigt der lokale Vergleich z_h = 0 / 2,5 / 5 µm bei 2·10⁵
  Paketen gleiche Größenordnungen, aber in den tiefen Bins noch Schwankungen um den Faktor 2 – die
  Abnahme braucht 10⁷ Pakete und mehrere Seeds (Server-Auftrag).
  kohärenter Boden von einigen Prozent (endlich viele Speckles).
- **Kalibrierläufe auf dem Server** (`docs/bericht-s9b2.md`, `tools/calibrate_bwpm.py`; 512² × bis
  3450 Ebenen, 16 Stufen statt 1024²/64 aus Kostengründen, beides geprüft): μs aus der Extinktion
  0,96 ± 0,03 von Born (einzelne Läufe 0,89–0,99), g vorwärts innerhalb 3,5 % von Born, 1 − g
  bei l_c ≤ 0,5 µm oder D < 5 um 20–35 % unter Born. Status: Kalibrierläufe erledigt, volle
  Spezifikation (1024², 64 Stufen) offen.

## T16b: Gemessene Streuparameter für den Monte Carlo (S9b-5)

`Calibration.scattering()` liefert (μs, Phasenfunktion) der BWPM-Realisierung: vorwärts die
gemessene Winkelverteilung; **rückwärts** – von der Vorwärtsausbreitung nicht erfasst, aber für OCT
und Reflexion entscheidend – Born in absoluten Werten (σ_Born/μs,gemessen); unter der
Winkelauflösung des Gitters (zwei k⊥-Zellen) ebenfalls Born. Ein erster Ansatz, der die
Rückwärtshälfte nur in der Form von Born übernahm und bei 90° an die Messung anschloss, unterdrückte
die Rückstreuung um vier Größenordnungen (die Vorwärtsausbreitung sieht bei 90° fast nichts) – das
OCT-Signal fiel entsprechend. Mit dem jetzigen Ansatz liegt das Hybrid-Signal (z_h = 2,5 µm) mit
gemessenen Parametern in derselben Größe wie mit Born (lokal, 5·10⁴ Pakete).

## T16c: Vektorielle Born-Referenz und Leistung mit E_z (S9b-7)

Anlass: BWPM K17 und Hinweise zu S9b-2 (BWPM #29, `docs/hinweise-scattuccino-s9b2.md`). Die
Referenz war skalar, die Messgröße ließ E_z weg; beides zusammen verdeckte teilweise die Grenze der
FPWPM bei großen Winkeln (BWPM E9: s-Anteil um r = cos θ (3 − cos θ)/2 zu klein, p-Anteil um
r/cos²θ zu groß).

- `microstructure.born(vector=True)` (Standard): σ(θ)·(1 + cos²θ)/2, gegen die geschlossene Form
  der exponentiellen Korrelation (Quadratur) auf 10⁻⁶; rückwärts und vorwärts unverändert, 90°
  halbiert. `vector=False` ist die bisherige skalare Rechnung (T16 unverändert).
- `microstructure.plane_wave_power`: |F{E_x}|² + |F{E_y}|² + |k⊥·F{E⊥}|²/k_z², ebene Welle bei 40°
  p-polarisiert 1/cos²θ, s-polarisiert 1 (10⁻⁹). Verwendet in `calibrate_with_bwpm` und
  `tools/calibrate_bwpm.py`.
- BWPM-Realisierung (96² × 60, l_c = 0,5 µm, D = 4, Seed 1; Python-Fallback):

  | | skalar, nur E⊥ (bisher) | vektoriell, mit E_z |
  | --- | --- | --- |
  | μs gefittet / μs Born | 0,960 | 1,010 |
  | g vorwärts gemessen / Born | 0,934 / 0,946 | 0,913 / 0,955 |

  μs liegt 1 % über der vektoriellen Born-Näherung (Test: 10 %). Das ist **kein** Beleg für
  vektorielle Streuung: Die FPWPM erhält die Quernorm Σ|E⊥|², nicht den Fluss; ihr Kohärenzverlust
  ist deshalb ∫ σ r²/cos θ dΩ über die Vorwärtshälfte (BWPM K17, Punkt 8, je Realisierung auf
  0,2–1,7 % bestätigt), für dieses Medium 0,978 × skalares = 1,028 × vektorielles Born. Eine erste
  Annahme aus Scattuccino (gestreute Querleistung ∫ σ r² dΩ, 0,932) war falsch. Das μs aus der
  Extinktion der BWPM ist damit nicht das physikalische μs. g weicht um −0,042 ab: die FPWPM
  streut mit vollem Feld bei großen Winkeln zu viel. Die K17-Vorhersage für Einfachstreuung gibt
  durch E_z eine Verschiebung von −0,014 (gemessen −0,021), absolut liegen die Messungen etwa 0,03
  darunter (Mehrfachstreuung bei μs z = 0,4, nicht geprüft). Der g-Vergleich ist als
  `xfail(strict)` mit Verweis auf BWPM E9 markiert (Freigabe des Maintainers 2026-11-19); wird die
  FPWPM korrigiert (BWPM #37), schlägt er an.
- Folge für den Hybrid: Born als Standard im Monte Carlo ist jetzt vektoriell (μs −4,9 %, g 0,935 →
  0,948 für das Medium aus #29); die Rückwärtshälfte der gemessenen Phasenfunktion (T16b) ebenfalls,
  bei 180° unverändert.

## T18: Kohärente Übergabe ohne Husimi-Bias (S9b-6)

`hybrid.split_coherent`: Projektion des Feldes an der Übergabe auf den frei propagierten Strahl
(`hybrid.free_field`); der kohärente Anteil kommt aus der Wigner-Funktion des Gaußstrahls
(`hybrid.gaussian_beam_packets`), nur der Rest über Husimi (Entwurf D8, Ergänzung; #31).

- Wigner-Pakete: Varianzen von Ort, Wellenvektor und Korrelation gleich dem Gaußstrahl
  (w(z)²/4, 1/w_f², Δz·Var u) innerhalb 1,5 % (5σ) bei 4·10⁵ Paketen.
- Projektion: |c|² zweier versetzter Gauß-Strahlen exakt (10⁻⁹), Leistung ohne Kreuzterm (10⁻¹²).
- **Einfachstreuung im engen Fokus** (w_f = 2 µm, z_f = 12 µm, z_h = 0, μs = 1 /mm isotrop, 6 Seeds à
  4·10⁵): Hybrid-Signal / geschlossene Form 0,994–1,016 bei 0,4–2 % Standardfehler (5–50 µm). Mit
  der bisherigen Husimi-Übergabe (σ = 0,6 µm) 0,85 → 0,37 – der Befund aus #29 direkt sichtbar.
  Vorwärts gerichtete HG (g = 0,9) gibt dasselbe Mittel, aber schwere Ausläufer (rückgestreut, dann
  vorwärts in die Mode: p(1)/p(−1) ≈ 7000); der Test nimmt isotrope Streuung.
- Kohärenter Anteil hinter 3 µm Zufallsmedium (nx = 96): 0,537 gegen Born e^(−μs z_h) = 0,549.
- **Lokal, nur orientierend** (nx = 128, Born, 3·10⁵ Pakete, Verhältnis zu z_h = 0 bei 16–22 µm):
  Husimi 1,28 (z_h = 5 µm) / 1,31 (10 µm), Wigner 0,91 / 0,77 (3 Seeds); eine Wiederholung bei
  10 µm mit 4 Seeds ergibt 0,93 ± 0,10 / 1,06 ± 0,24 / 0,85 ± 0,10 für σ = 0,25 / 0,4 / 0,8 µm
  (z_h = 0 selbst ±9 %). Keine σ-Abhängigkeit nachweisbar; ob neben der Husimi-Glättung ein
  zweiter Effekt bleibt, entscheidet erst die Abnahme mit 10⁷ Paketen (T17).

## T10b: Erweitertes Huygens-Fresnel, geschlossene Form (S9c-1)

`scattuccino.ehf` nach thrane2000 (Entwurf `docs/entwurf-s9c-ehf.md`, D11 A); Test
`tests/test_s9c1_ehf.py`.

- Ψ(z) in geschlossener Form gegen die Quadratur von thrane2000 Gl. 21 auf 10⁻⁸ (Fokus auf der
  Grenzfläche, w_H²/w_S² zwischen 0,05 und 0,95, damit alle drei Terme zählen). Die Integration
  gibt beim Kreuzterm den Faktor 4 wie levitz2004 Gl. 2; thrane2000 Gl. 22 druckt 2 (#36, Label
  `paper`).
- ρ₀ gegen die Quadratur von Gl. A1 (B(s) = s/n im Gewebe, ohne Streuung bis zur Linse) auf 10⁻¹⁰.
- Grenzfälle, diffus und spiegelnd: μs → 0 gibt Ψ → 1; ohne seitliche Kohärenz (θ_rms → ∞)
  Ψ → e^(−2μs z) (je 10⁻⁶).
- thrane2000 Abb. 2 (μs = 20/mm, Fokus auf der Grenzfläche; diffus mit und ohne
  Shower-Curtain-Effekt, spiegelnd) und Abb. 3 (z = 0,5 mm, μs = 2–7,5/mm, Glasplatte nach Gl. 27):
  log₁₀ Ψ aus den Abbildungen digitalisiert (300 dpi, Rahmen als Achsen). Die Rechnung liegt
  durchweg 0,08–0,15 Dekaden darüber; die Abweichung ist bei allen Kurven gleich gerichtet und
  kommt vom Ablesen der Achsen. Toleranz 0,25 Dekaden.
- Bei den Parametern der Abbildungen ist w_H²/w_S² ≪ 1; ob dort mit 2 oder 4 gerechnet wurde, ist
  nicht zu unterscheiden.

## T10b: Allgemeines EHF mit der Phasenfunktion des Mediums (S9c-2)

`scattuccino.ehf` (`overlap_width`, `phase_correlation`, `spherical_mtf`, `general_efficiency`);
Test `tests/test_s9c2_ehf.py`. Geometrie des Hybrids: Medium mit angepasster Brechzahl ab der
Oberfläche, Strahl U = exp(−a p²) an der Oberfläche, Wellenzahl k im Medium. Das Signal ist das
Überlappintegral der mittleren Intensitäten von Beleuchtung und Detektionsmode (lutomirski1971
Gl. 10; thrane2000 Gl. 11, 15, 16), nach Parseval Ψ = ∫ M_s² e^(−ϱ²/A²) d²ϱ/(πA²) mit
M_s = exp{−μs z [1 − ∫₀¹ b(ϱu) du]} (lutomirski1971 Gl. 11, Kugelwelle) und
b(ρ) = 2π ∫ p(cos θ) J₀(kρ sin θ) sin θ dθ über die Vorwärtshälfte.

- A gegen w_H aus thrane2000 Gl. 13 (w_H² = 2z²/(k²A²)) auf 10⁻¹².
- Mit der Gaußschen MTF von thrane2000 Gl. 20 und ρ₀ aus S9c-1 reproduziert die allgemeine Form
  die geschlossene Form (n = 1, Linse an der Oberfläche) auf 10⁻⁶ – auch den Kreuzterm 4 (#36).
- b(ρ) einer schmalen Gaußkeule (θ_r = 0,05) gegen exp(−k²θ_r²ρ²/4) auf 2·10⁻³; M_s einer
  Gaußschen Korrelation gegen √π erf(c)/(2c) auf 10⁻⁶; ohne Streuung Ψ = 1.

**Befund Brechzahl:** In thrane2000 Gl. A1 steht in der Streuphase die Vakuum-Wellenzahl k mit
der reduzierten Länge B(s) = s/n. Ist θ_rms der Streuwinkel im Gewebe (so in den Beispielen,
θ_rms = arccos g), ist ρ₀ dort um den Faktor n zu groß; in der Herleitung hier (Wellenzahl im
Medium) tritt n nicht gesondert auf. Für das Medium unten ändert das Ψ um den Faktor 1,02 (5 µm) bis 1,8 (≥ 30 µm)
(Issue mit Label `paper`).

**Medium aus #29** (vektorielles Born: μs = 190,3/mm, g = 0,948; Vorwärtsanteil b(0) = 0,995),
Strahl der Hybridläufe (w_f = 2 µm, Fokus 12 µm, n₀ = 1,38, λ = 0,8 µm). Dokumentiert, ohne
Toleranz (D11):

| z (µm) | Ψ allgemein | geschlossen, θ_e = 0,148 (levitz2004) | dasselbe mit k₀ wie gedruckt | geschlossen, θ = arccos g = 0,325 | e^(−2μs z) |
| --- | --- | --- | --- | --- | --- |
| 5 | 0,822 | 0,955 | 0,976 | 0,821 | 0,149 |
| 10 | 0,423 | 0,641 | 0,769 | 0,285 | 0,022 |
| 20 | 0,090 | 0,176 | 0,288 | 0,043 | 4,9·10⁻⁴ |
| 30 | 0,034 | 0,083 | 0,148 | 0,018 | 1,1·10⁻⁵ |
| 40 | 0,020 | 0,057 | 0,103 | 0,012 | 2,4·10⁻⁷ |
| 60 | 0,011 | 0,038 | 0,070 | 0,008 | 1,2·10⁻¹⁰ |

Die Whittle-Matérn-Phasenfunktion (scharfe Spitze, breiter Fuß) liegt zwischen den beiden
Gaußschen Näherungen; keine trifft sie auf besser als einen Faktor 1,5–3. Ab z ≈ 10 µm dominiert
nach dem EHF das vorwärts mehrfach gestreute Licht (Ψ ≫ e^(−2μs z)). Das ist die Größe, gegen die
S9c-3 die Hybriddaten aus #34 vergleicht (D12).

## T16c (Fortsetzung): δn als Born-Inkrement in BWPM (S9b-8)

BWPM B7 (BWPM #29, 2026-11-19): δn je Scheibe als vektorielles Born-Inkrement auf einem Wirt n₀
statt als Brechzahlstufen; `microstructure.bwpm_medium(model=…)`, Option `model` in
`calibrate_with_bwpm`, `hybrid_oct` und `tools/calibrate_bwpm.py`. **Standard bleibt `"steps"`.**

- Dieselbe Realisierung wie oben (96² × 60, l_c = 0,5 µm, σ_n = 0,0572, dx = dz = 0,1 µm, Seed 1;
  BWPM fe2228e, Python-Fallback): Mit Stufen bleibt Σ|E⊥|² auf 0,5 % erhalten, mit dem Inkrement
  wächst sie auf 1,07 / 1,28 / 2,00 nach 2 / 4 / 6 µm. Die Leistung geht in Wellen nahe dem
  Grenzwinkel: |E⊥|² bei cos θ < 0,3 ist 0,031 / 0,198 / 0,865 (Stufen 0,0004 / 0,0006 / 0,0010),
  wo das vektorielle Born kaum hinstreut.
- Abhängigkeit (Zuwachs von Σ|E⊥|² nach 6 µm): σ_n = 0,0572 / 0,0286 / 0,0143 bei dz = 0,1 µm
  100 % / 4,8 % / 1,1 %; bei σ_n = 0,0572 und dz = 0,05 / 0,025 µm 18 % / 10 %. Kleineres dz dämpft
  den Zuwachs, beseitigt ihn nicht. BWPM hat das Inkrement bis σ_n ≈ 0,03 und dz ≤ 0,04 µm geprüft
  und Wellen nahe dem Grenzwinkel ausdrücklich offen gelassen.
- Test `test_bwpm_born_increment_conserves_power` (Leistung wie bei den Stufen auf 2 %) als
  `xfail(strict)` mit Verweis auf BWPM #29; der g-Test mit Stufen bleibt `xfail` (E9). Ist das
  Inkrement stabil, wird der Standard umgestellt und der g-`xfail` entfernt (Freigabe des
  Maintainers 2026-11-19).

## T10c: Rückstreuung der BWPM an δn gegen Born erster Ordnung (S9c-4a)

BWPM B8 (BWPM #29, #54): Rückstreuung an `fluctuation=δn`. `microstructure.born_reflection`:
Born erster Ordnung der Reflexion einer Schicht mit je Scheibe konstantem δn (Weyl-Entwicklung der
dyadischen Green-Funktion, bornwolf1999), dieselbe Realisierung wie in der BWPM. Test
`tests/test_s9c4_backscatter.py`.

- `born_reflection` gegen die exakte Reflexion einer homogenen Schicht (Airy, δ = 10⁻⁴, 1–30
  Scheiben): Abweichung ≤ 3δk₀d relativ zur Skala δ/n₀ – Born behält die Phase k₀n₀ des Wirts.
- BWPM (64², 40 Scheiben à 0,1 µm, l_c = 0,5 µm, iterations = 1, |cos θ| > 0,28), Seeds 1–3:

  | σ_n | Leistung BWPM/Born | Feldfehler |
  | --- | --- | --- |
  | 0,0018 | 0,995–1,006 | 0,031–0,039 |
  | 0,0036 | 0,996–1,019 | 0,062–0,079 |
  | 0,0143 | 1,009–1,105 | 0,243–0,311 |

  Der Feldfehler ist linear in σ_n (Verhältnis 1,99–2,01 bei 2×, 3,92–3,94 bei 4×): die
  Phasenmodulation durch Vorwärtsstreuung auf Hin- und Rückweg, die Born erster Ordnung nicht
  enthält. Die Leistung folgt im Born-Bereich auf 2 %. **Die Rückstreuung der BWPM ist damit für
  die Realisierung bestätigt.**
- Auswertung: Zellen nahe dem Grenzwinkel ausgenommen; dort hat das Born-Feld den Faktor 1/k_z,
  und im periodischen Kasten liegt zufällig ein Mode fast auf k⊥ = k (einzelne Zellen dominieren).
- **Offen (S9c-4b):** Ensemble gegen die kontinuierliche Born-Phasenfunktion und gegen den Monte
  Carlo. Erste Läufe (6,4 µm Kasten, 2–6 Realisierungen) geben 0,93 / 1,10 / 1,33 für h = 0,1 /
  0,05 / 0,025 µm – das Spektrum des Gitterfeldes bei q_z = 2k hängt von h ab (Abtastung, je Scheibe
  konstant), und die Statistik je Winkel ist zu klein. Dafür braucht es feine Gitter und viele
  Realisierungen (Server). Die Leistung nach Rückstreuung wird mit `plane_wave_power(Feld)`
  ausgewertet (die Funktion transformiert selbst).
