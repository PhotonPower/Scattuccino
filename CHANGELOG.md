# Changelog

Format nach [Keep a Changelog](https://keepachangelog.com/de/1.1.0/).

## [Unreleased]

### Hinzugefügt

- S0 begonnen: Plan mit Stufen S0–S12, Konventionen, Abhängigkeiten, Testfälle, Quellen, Regeln
  für Agenten, CI.
- S0: Spezifikation (`docs/spezifikation.md`).
- S1: Monte-Carlo-Referenz in Schichten (`scattuccino.layers`) und unabhängige Referenzen
  (`scattuccino.reference`: ohne Streuung, H-Funktion); Tests T1a, T1b, T2.
- S1b: Adding-Doubling (`scattuccino.adding_doubling`) als unabhängige Referenz für Schichten mit
  HG und Brechzahlsprüngen; Monte Carlo gegen Adding-Doubling (T1c).
- S2: Photonenprotokoll (`scattuccino.transport`), Volumengitter, mikroskopisches Beer-Lambert
  mit exakter Umgewichtung, Detektoren (`scattuccino.detectors`), Diffusionsreferenzen
  (`scattuccino.diffusion`, Kienle-Patterson); Tests T3a–T3d.
- S3: Phasenfunktionen (`scattuccino.phase`: HG, MHG, Gegenbauer, tabelliert, Mie) mit Momenten und
  γ, δ; Mie-Reihe (`scattuccino.mie`); `Layer(phase=…)` in Transport und Adding-Doubling; Tests
  T4, T5.
- S4a: Voxelgeometrie (`scattuccino.voxels`) mit zählerbasierten Zufallszahlen je Paket
  (`scattuccino.random`); Tests T6a.
- S5a: C++-Kern (`scattuccino._core`, scikit-build-core, pybind11, OpenMP) für den Voxeltransport,
  je Paket gleich der Python-Referenz; Fassade `scattuccino.fast`; CI baut den Kern, Job ohne Kern.
- S5b: gemeinsamer Paketcode für CPU und GPU (`cpp/include/scattuccino/photon.hpp`), CUDA-Kernel
  (`cpp/src/core_cuda.cu`, ein Thread je Paket), Backend `cuda` in `scattuccino.fast`, Test T14 und
  Anleitung für den Server (`docs/agent-server-s5b.md`).
- S5b auf dem Server bestätigt; S5c: persistente Threads, Absorption in Registern, Voxelwertung und
  Protokoll abschaltbar; Summen ohne Protokoll (`VoxelSimulation.totals`); Messauftrag
  `docs/agent-server-s5c.md`.
- S5d: FP32-Option `precision="single"` für CPU- und CUDA-Kern (Paketcode als Template auf die
  Genauigkeit, double bitgleich wie zuvor), Zähler `n_exit`/`n_roulette`, Test T15; auf der
  GPU 11- bis 13-mal schneller als in double.
- S4b: Tetraedernetze im C++-Kern (`photon_tet.hpp`, `tets.cpp`), Netzhelfer und Fassade
  (`scattuccino.tets`); Tests T6b (je Paket gleich dem Voxelgang, Glaskugel gegen geometrische Optik).
- S6a: polarisierter Transport (`scattuccino.polarisation`: Mueller-Tabellen für Mie und Rayleigh,
  Streuschritt, Schicht/Halbraum mit angepasster Brechzahl, Stokes im Laborsystem); Tests T7a, T7b.
- S6b: polarisierte Fresnel-Ereignisse an Brechzahlsprüngen (`interface`, `fresnel_amplitudes`), Schicht
  mit Brechzahlen und schrägem Einfall in `simulate_polarised`; Tests T7c. D7 (FP32-Option) freigegeben.
- S6c: Polarisation im C++-Kern (`photon_pol.hpp`, `simulate_voxels_polarised`), Python-Referenz auf
  zählerbasierten Zufall umgestellt, Fassade `fast.simulate_polarised_slab`; Test T7e (je Paket gleich).
- S7a: Fluoreszenz und Raman (`scattuccino.fluorescence`), innere isotrope Quellen und
  Absorptionsereignisse im Schichttransport (`transport.Sources`, `record_absorption`); Tests T8.
- S8a: Impulsübertrag je Schicht im Protokoll, `scattuccino.dynamics` (g1, g2, Speckle-Kontrast),
  Korrelationsdiffusion; Tests T9.
- S9a: OCT-A-Scans aus dem Protokoll (`scattuccino.oct`: Detektor, Tiefe aus optischem Weg,
  Kohärenz-Gating, exakte Einfach-Rückstreuung); Tests T10a.
- S9b-2: Zufallsmedium (`scattuccino.microstructure`: Whittle-Matérn, Born, Zufallsfeld, Kalibrierung
  mit BWPM); Tests T16; Auftrag für die Kalibrierläufe (`docs/agent-server-s9b2.md`).
- S9b-2: Kalibrierläufe des Zufallsmediums in BWPM (`tools/calibrate_bwpm.py`, Bericht
  `docs/bericht-s9b2.md`, Rohdaten unter `tools/out/`).
