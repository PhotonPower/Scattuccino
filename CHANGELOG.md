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
