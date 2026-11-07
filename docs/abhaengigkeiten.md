# Abhängigkeiten

> Grundsatz des Maintainers (2026-11-07): so wenig externe Abhängigkeiten wie möglich; bevorzugt
> die Programme im eigenen GitHub-Konto (PhotonPower), notfalls werden sie erweitert.

## Erlaubt ohne Rückfrage

| Paket | Zweck |
| --- | --- |
| Python ≥ 3.10, NumPy, SciPy | Referenzimplementierung, Numerik |
| pytest, ruff, mypy | Entwicklung (`requirements-dev.txt`) |
| CMake, C++20-Compiler, CUDA ≥ 12.6 | ab dem C++/CUDA-Kern, wie BWPM |

## Programme des Maintainers (lesend; Erweiterungen als Vorschlag)

| Programm | Nutzung in Scattuccino |
| --- | --- |
| Raytatouille | Konventionen, Polarisation (`rtt-polar`, PRT/Mueller), Strahlbündel als Quellen und Detektoren, radiometrischer Port mit M9, BSSRDF-Tabellen im sequenziellen Tracer |
| RTT-Fusion | Port-Typen, Adapter-Schnittstelle, Querproben |
| Fieldtracing | Fokusfelder (OCT, Mikroskopie), Übergang Wellenoptik → Transport |
| CliffordBem | Phasenfunktionen und Mueller-Matrizen einzelner Streuer (Zellen, Organellen, Partikel), Mie-Querprobe |
| SpeckleBem | raue Gewebeoberfläche (BRDF), Speckle an der Oberfläche |
| BWPM | dünne, schwach streuende Proben und die ersten freien Weglängen im Fokus; Zufallsmedium (S9b-2, optional importiert) |
| hp-FEM | Vorschlag: skalare H1-Löser für Diffusion (Frequenzbereich = komplexe Helmholtz-Gleichung), Akustik (Photoakustik) und Wärme (Lasertherapie) |
| BioOptoMD | Brechzahl- und Absorptionsprofile molekularer Schichten (Datenport N) |
| BindOptiX (PhotoKinX) | Kopplung an Mikrofluidik und thermo-optische Effekte bei Biosensoren (prüfen) |
| Mikroptim | Erfahrung mit CUDA-Raytracing für den GPU-Kern (prüfen) |

## Bewusst nicht verwendet

- Fremde Monte-Carlo-Codes (MCX, MMC, PyXOpto, MCML …) nicht als Abhängigkeit und nicht als Code; ihre veröffentlichten Ergebnisse dienen als Vergleich. MCML-Tabellen und Adding-Doubling werden als unabhängige Referenz **selbst** implementiert.
- Fremde FEM-Pakete (FEniCS, NIRFAST, TOAST++): stattdessen hp-FEM erweitern.
- Akustik-Pakete (k-Wave): Photoakustik bis zur Druckverteilung hier, Ausbreitung über hp-FEM (Vorschlag).

Jede neue externe Abhängigkeit braucht eine Entscheidung im Plan (`docs/roadmap.md`).
