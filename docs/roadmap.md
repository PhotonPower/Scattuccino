# Plan und Stufen

> Erstfassung 2026-11-07. Entscheidungen des Maintainers (offen: D1–D6) stehen am Ende.

## Ziel und Prinzip

Ein Kern – polarisierter, zeit- und weglängenaufgelöster Photonentransport (Monte Carlo) in
Schichten, Voxeln und Tetraedernetzen – schreibt für jedes detektierte Photon ein **Protokoll**
(Austrittsort und -richtung, Gewicht, Gesamt- und Teilweglänge je Gewebe, Zahl der Streuungen,
maximale Tiefe, Impulsübertrag je Gewebe, Stokes-Vektor, bei Bedarf Phase). Die Methoden der
Medizintechnik sind **Detektoren**, die dieses Protokoll auswerten, ohne neu zu rechnen.
Volumenausgaben (Fluenz, absorbierte Energie) bedienen Photoakustik und Therapie.

## Methoden und Stufen

| Methode | Stufe |
| --- | --- |
| diffuse Reflexion/Transmission (CW), Pulsoximetrie/PPG | S1, S2 |
| NIRS/fNIRS (CW, Frequenz- und Zeitbereich), Teilweglängen | S2 |
| DRS, SFDI, hyperspektrale Bildgebung, Dermatoskopie, Narrow-Band-Imaging | S2, S10, S11 |
| Polarimetrie, Mueller-Matrix-Bildgebung | S6 |
| Fluoreszenz (ICG, PpIX, Autofluoreszenz), Raman, SORS | S7 |
| DCS, Laser-Speckle-Kontrast, Laser-Doppler | S8 |
| OCT (TD/SD/SS), OCT-Angiographie, PS-OCT; konfokal, Multiphoton | S9 |
| DOT, Fluoreszenztomographie (Jacobi-Matrizen) | S10 |
| Photoakustik (Fluenz → Anfangsdruck), PDT/Lasertherapie/LITT/Optogenetik (Dosis) | S2 (Fluenz), S12 (Kopplung) |
| Streulicht im Auge | S3, S12 |

## Stufen

| Stufe | Inhalt | Abnahme |
| --- | --- | --- |
| **S0 Spezifikation** | Strahlungstransportgleichung, Größen und Einheiten, Konventionen (Raytatouille), Format des Photonenprotokolls, Quellen- und Detektormodell, Gewebeparameter; Volltext der Kernquellen | Spezifikation freigegeben, Quellen geprüft |
| **S1 CPU-Referenz in Schichten** | Monte Carlo in Python/NumPy (Schichten, Fresnel an Grenzen, Henyey-Greenstein, Gewichtsverfahren mit Russischem Roulette); **eigene** unabhängige Referenz: Adding-Doubling (van de Hulst, Prahl) | Reflexion/Transmission gegen Adding-Doubling und veröffentlichte MCML-Tabellen, Energieerhaltung |
| **S2 Photonenprotokoll und Detektoren** | ortsaufgelöste Reflexion, Zeit- und Frequenzbereich, Teilweglängen, Fluenz und absorbierte Energie | Diffusionslösungen (Farrell, Kienle-Patterson für Schichten), Zeitbereich (Patterson), Fluenz gegen Diffusion |
| **S3 Phasenfunktionen und Mikrostruktur** | modifizierte HG, Gegenbauer, tabellierte Phasenfunktionen, Mie (eigene Reihe und CliffordBem), anisotrope Zylinderstreuer; Parameter γ, δ für den subdiffusiven Bereich | Mie-Momente, Querprobe CliffordBem, subdiffusive Reflexion gegen Literatur |
| **S4 Geometrie** | S4a Voxel (✅, zählerbasierter Zufall je Paket); S4b Tetraedernetze (✅, im C++-Kern; Netze aus Nurbforge/gmsh folgen mit der Fusion), gekrümmte Grenzen mit Brechzahlsprung | Schichtfall auf Voxeln/Netz = S1; analytische Einschlüsse (Störungstheorie) |
| **S5 C++/CUDA-Kern** | wie BWPM B4: Python-Referenz bleibt, C++-Kern mit CPU-Pfad (S5a ✅), CUDA (S5b ✅, S5c ✅) und FP32-Option (S5d ✅, GPU 21–36× schneller als 32 CPU-Kerne), Sollwertdateien | je Paket gleich der Python-Referenz (T13); Laufzeit |
| **S6 Polarisation** | S6a ✅ Stokes-Vektoren, Mueller-Matrizen je Streuung, mitgeführtes Bezugssystem (Python-Referenz); S6b ✅ Fresnel-Mueller an Brechzahlsprüngen; S6c ✅ C++-Kern (CPU, Voxel); offen: GPU (mit FP32), Konventionen von `rtt-polar`, Benchmarks Ramella-Roman | Benchmarks von Ramella-Roman, Prahl, Jacques; Mie-Mueller |
| **S7 Fluoreszenz und Raman** | S7a ✅ zweistufig in Schichten (Python): Quellen aus der Absorption, Quantenausbeute, Re-Absorption; offen: Emissionsspektren, Voxel/C++/GPU, Fluoreszenzlebensdauer im Zeitbereich | Fluoreszenz-Diffusionslösungen, Energiebilanz |
| **S8 Dynamik** | S8a ✅ Impulsübertrag je Schicht → g₁(τ), g₂, Speckle-Kontrast (Brownsche Bewegung, ungerichtete Strömung); offen: Dopplerspektrum (LDF), Voxel/C++/GPU, Gewebemodelle mit Gefäßanteil | Korrelationsdiffusion (Boas, Yodh) |
| **S9 OCT und Fokus** | S9a ✅ weglängenaufgelöste Rückstreuung, Kohärenz-Gating, Einfach-/Mehrfachstreuung (Intensitätsmodell); S9b Fokus aus Fieldtracing und Hybrid mit BWPM (Entwurf `docs/entwurf-s9b-hybrid.md`, D8–D10 freigegeben; S9b-1 Phasenraum ✅ in RTT-Fusion; S9b-2 Zufallsmedium ✅, Kalibrierläufe auf dem Server; S9b-3 Hybrid-OCT); S9c erweitertes Huygens-Fresnel, Speckle | Einfachstreumodell und EHF (Thrane, Schmitt), Querprobe BWPM im schwach streuenden Grenzfall |
| **S10 Inverse Probleme und schnelle Modelle** | adjungiertes und Störungs-MC (Jacobi-Matrizen), Skalierung einer Simulation auf viele (μa, μs) („White MC“), Nachschlagetabellen, SFDI-Inversion; Diffusion im Frequenzbereich über hp-FEM (Vorschlag) | Jacobi gegen finite Differenzen, Skalierung gegen direkte Läufe |
| **S11 Gewebe- und Phantomdatenbank** | Chromophorspektren (HbO₂, Hb, Wasser, Lipid, Melanin, Bilirubin, Kollagen) mit Quellen, Gewebemodelle (Haut, Hirn, Brust, Muskel, Auge), Phantomrezepturen (Lipid, TiO₂/Silikon, Mikrosphären, Tinte), inverse Adding-Doubling für die Charakterisierung | Quellenwerte reproduziert, IAD gegen eigene Vorwärtsrechnung |
| **S12 Fusion** | Adapter für RTT-Fusion: Strahlen ↔ Photonenpakete, radiometrischer Port (mit M9), BSSRDF-Tabelle für den sequenziellen Tracer; Datenports zu hp-FEM (Akustik, Wärme), BioOptoMD, Voxelspaetzle | Konformitätssatz der Fusion; gemeinsamer Pilot mit Voxelspaetzle |

**Gemeinsamer Pilot mit Voxelspaetzle:** gedruckte Mikrooptik auf der Faser als Endoskop- bzw.
OCT-Sonde, die in ein Gewebephantom abbildet (Druckmodell → Sonde → Gewebe → Signal).

## Entscheidungen (offen, Issues mit Label `frage`)

| Nr. | Frage | Vorschlag |
| --- | --- | --- |
| D1 | Sprache und Kern | Python-Referenz zuerst (S1–S4), dann C++/CUDA-Kern wie BWPM (S5); Python-Anbindung über pybind11 |
| D2 | Reihenfolge der Methoden nach S5 | nach klinischer Verbreitung: Pulsoximetrie/NIRS → DRS/SFDI → Fluoreszenz → OCT → DCS/LSCI → Polarimetrie |
| D3 | Polarisation | Stokes-Vektor von Anfang an im Protokoll (abschaltbar), Rechnung ab S6 |
| D4 | Vergleich mit GPL-Codes | nur über veröffentlichte Ergebnisse, kein Code |
| D5 | Diffusion, Akustik, Wärme | hp-FEM um skalare H1-Probleme erweitern (Vorschlag dort) statt eigener FEM |
| D7 | FP32 im GPU-Kern | freigegeben 2026-11-07: FP64 Standard und Referenz, FP32 als `precision="single"`, Abnahme statistisch |
| D8 | Feld → Pakete bei S9b | freigegeben: Husimi-Phasenraum (Befund: Übergabe an/hinter den Fokus) |
| D9 | Übergabetiefe, Konsistenz oben/unten | freigegeben: z_h ≈ 1–2/μs, Zufallsmedium in BWPM auf μs, g kalibriert |
| D10 | Detektion in S9b | freigegeben: Modenanpassung im Phasenraum (ballistisches Licht über das Feld) |
| D6 | Phantome und Messdaten | Validierung zunächst an Literatur; eigene Phantommessungen, sobald vorhanden |
