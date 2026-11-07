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
