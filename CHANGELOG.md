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
- S9b-3a: Hybrid-OCT-Kette (`scattuccino.hybrid`), innere Quellen mit Richtung (`Sources.direction`);
  Tests T17a; Befund: Next-Event-Schätzer nötig.
- S9b-4: Next-Event-Schätzer der Modenkopplung (`hybrid.next_event_signal`), Streuereignisse im
  Transport (`record_scattering`, auch stapelweise), allgemeine Gauß-Akzeptanz; Tests T17b;
  Abnahmeauftrag `docs/agent-server-s9b4.md`.
- S9b-2: Kalibrierläufe des Zufallsmediums in BWPM (`tools/calibrate_bwpm.py`, Bericht
  `docs/bericht-s9b2.md`, Rohdaten unter `tools/out/`).
- S9b-5: Monte Carlo im Hybrid mit den in der BWPM gemessenen Streuparametern
  (`Calibration.scattering()`, `hybrid_oct(scattering=…)`); Raumwinkel-Jacobi in der g-Auswertung;
  Abnahmeauftrag mit Kalibrierschritt.
- S9c begonnen: Entwurf erweitertes Huygens-Fresnel als OCT-Referenz (`docs/entwurf-s9c-ehf.md`,
  D11–D13 zur Freigabe).
- S9c-1: geschlossene Form des EHF-Modells nach Thrane (`scattuccino.ehf`: `heterodyne_efficiency`,
  `beam_radii`, `lateral_coherence_length`; diffus und spiegelnd, Shower-Curtain-Effekt), Kreuzterm
  mit Faktor 4 (#36); Test T10b gegen Quadraturen und thrane2000 Abb. 2 und 3.
- S9c-2: allgemeines EHF mit der Phasenfunktion des Mediums (`ehf.general_efficiency`,
  `phase_correlation`, `spherical_mtf`, `overlap_width`) für die Geometrie des Hybrids; reproduziert mit
  Gaußscher MTF die geschlossene Form; Vergleich für das Medium aus #29 dokumentiert (T10b).
- S9b-6: kohärente Übergabe ohne Husimi-Bias (#31): Projektion des Feldes an der Übergabe auf den
  frei propagierten Strahl (`hybrid.split_coherent`, `hybrid.free_field`), kohärenter Anteil aus der
  Wigner-Funktion des Gaußstrahls (`hybrid.gaussian_beam_packets`), nur der Rest über Husimi;
  `hybrid_oct(coherent="wigner")` (Standard, `"husimi"` wie bisher), `HybridResult.power_coherent`;
  Test T18 (Einfachstreuung im engen Fokus gegen geschlossene Form).
- S9b-7: vektorielle Born-Referenz `microstructure.born(vector=True)` (Standard, Dipolfaktor
  (1 + cos²θ)/2; `vector=False` skalar wie bisher) und Leistung der gestreuten ebenen Wellen mit E_z
  (`microstructure.plane_wave_power`) in `calibrate_with_bwpm` und `tools/calibrate_bwpm.py`, nach
  BWPM K17/E9; Test T16c, g-Vergleich mit der BWPM als bekannte Abweichung (`xfail`). D9 ergänzt:
  Monte Carlo im Hybrid mit vektoriellem Born, gemessene Werte als Diagnose, solange BWPM #37 offen
  ist.
- S8b: Laser-Doppler mit gekennzeichneten bewegten Streuern (`Layer.moving_fraction`,
  `Layer.phase_moving`, `Protocol.momentum_moving`, `Protocol.n_moving`), `dynamics.g1_tagged`,
  `unshifted_fraction`, `doppler_spectrum`, `mean_frequency`; Referenz nach Bonner und Nossal
  (`reference.bonner_*`); Test T9b.
- S9b-8: δn in BWPM wahlweise als vektorielles Born-Inkrement (`microstructure.bwpm_medium`,
  `model="increment"` in `calibrate_with_bwpm`, `hybrid_oct`, `tools/calibrate_bwpm.py`; BWPM B7);
  Standard bleibt `"steps"`, weil das Inkrement bei σ_n ≈ 0,06 Leistung in Wellen nahe dem
  Grenzwinkel pumpt (T16c, `xfail`, BWPM #29).
