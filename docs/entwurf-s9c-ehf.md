# Entwurf S9c: erweitertes Huygens-Fresnel (EHF) als analytische OCT-Referenz

> Entwurf, Stand 2026-11-19. D11 (Option C) und D12 freigegeben 2026-11-19; D13: Volltexte ausstehend. Formeln hier nur in ihrer
> Struktur; Gleichungen und Gleichungsnummern werden erst nach dem Lesen der Volltexte
> übernommen (Regel 6). Alle genannten Quellen haben den Status **Z**.

## Kontext

Der Hybrid (S9b) wird bisher nur in sich geprüft: Unabhängigkeit von der Übergabetiefe z_h
(T17, Abnahme #34) und Einfachstreuung im engen Fokus gegen eine geschlossene Form (T18). Es fehlt
eine **unabhängige Referenz für das OCT-Signal mit Fokus und Mehrfachstreuung**. Das liefert das
EHF-Modell von Thrane, Yura und Andersen (thrane2000), aufbauend auf Schmitt und Knüttel
(schmitt1997) und der EHF-Ausbreitung durch Zufallsmedien (lutomirski1971):

- **Annahmen:** Streuung vorwärts unter kleinen Winkeln (paraxial, g nahe 1), Rückstreuung einmal
  an der Tiefe z (Rückstreukoeffizient), Mehrfachstreuung nur auf Hin- und Rückweg, die Probe
  liegt im Fokus, Kohärenzlänge kurz gegen die Weglängenstreuung der Vorwärtsstreuung.
- **Ergebnis (Struktur):** mittleres heterodynes Signal ⟨i²(z)⟩ = ⟨i²⟩₀ · Ψ(z) mit drei Termen:
  ungestreut hin und zurück (∝ e^(−2μs z)), Kreuzterm ungestreut/gestreut (∝ e^(−μs z)(1 − e^(−μs z))),
  gestreut hin und zurück (∝ (1 − e^(−μs z))²). Die Gewichte der beiden letzten Terme hängen vom
  Verhältnis der Strahlradien mit und ohne Streuung ab (w_S/w_H); dort geht die seitliche
  Kohärenzlänge ρ₀(z) ein, und über sie μs und die Winkelbreite der Streuung.

Der erste Term ist im Grenzfall μs z ≪ 1 genau die Einfachstreuung im Fokus, die T18 bereits als
geschlossene Form prüft. Das verbindet S9c mit vorhandenen Tests.

## D11: Welche Form des EHF

| Option | Inhalt | Bewertung |
| --- | --- | --- |
| A geschlossene Form (thrane2000) | Gaußsche Kleinwinkelnäherung der Phasenfunktion (θ_rms), Ψ(z) als Formel | direkt mit veröffentlichten Kurven vergleichbar; θ_rms ist für die Whittle-Matérn-Phasenfunktion (Spitze plus breiter Fuß) nur eine Näherung |
| B allgemeines EHF | gegenseitige Kohärenzfunktion mit der Phasenstrukturfunktion des Mediums, berechnet aus der tatsächlichen Phasenfunktion (Kleinwinkelteil), Integrale numerisch | gleiche Annahmen wie A, aber ohne θ_rms; passt zum Medium der Hybridläufe |
| **C = A + B (Vorschlag)** | A als Prüfstein gegen die veröffentlichten Kurven, B als Referenz für den Hybrid | B muss im Grenzfall der Gaußschen Phasenfunktion A exakt reproduzieren; das prüft die numerische Umsetzung unabhängig |

## D12: Vergleichsgröße Hybrid ↔ EHF

Der Hybrid liefert das in die Detektionsmode gekoppelte, kohärenzgefensterte Signal je mm,
normiert auf die eingestrahlte Leistung; das EHF liefert ⟨i²(z)⟩ bis auf die Systemkonstante.

- **Vorschlag:** Verglichen wird der **Verlauf** S(z)/S_einfach(z), also das Signal geteilt durch
  die geschlossene Form der Einfachstreuung (T18) bei derselben Tiefe, gegen Ψ(z)/e^(−2μs z). Die
  Systemkonstante, der Rückstreukoeffizient und die Fokusfunktion kürzen sich heraus; übrig bleibt
  der Beitrag der Mehrfachstreuung.
- Zusätzlich ein absoluter Vergleich bei kleinem μs z (dort ist Ψ ≈ e^(−2μs z) und der Hybrid ist
  durch T18 absolut geprüft).

## D13: Volltexte und Abgrenzung

- **Volltexte nötig** vor jeder Formel: thrane2000 (Ψ(z), w_H, w_S, ρ₀), schmitt1997
  (Modellrahmen, Rückstreuung), lutomirski1971 (EHF, Strukturfunktion), dazu levitz2004
  (Korrekturen und Anwendung des Modells auf Messungen). Der Maintainer stellt sie bereit oder gibt
  den Zugang frei; bis dahin bleibt S9c bei diesem Entwurf.
- **Grenzen des EHF**, die die Vergleiche bestimmen: nur Kleinwinkelstreuung (das Medium aus #29
  hat g ≈ 0,95; geeignet), Rückstreuung einmal, keine Mehrfachrückstreuung (diffuser Hintergrund
  fehlt in großer Tiefe), Probe im Fokus. Abweichungen des Hybrids **jenseits** dieser Annahmen
  sind kein Fehler, sondern werden als Befund dokumentiert.
- **Speckle** (Statistik der einzelnen Realisierung, zweiter Teil von S9c in der Roadmap) folgt
  nach S9c-3 in einem eigenen Entwurf.

## Stufen nach Freigabe

1. **S9c-1** geschlossene Form nach thrane2000 (D11 A) in `scattuccino.ehf`; Tests: Grenzfall
   μs → 0 gegen die Einfachstreuung von T18, Grenzfälle der Strahlradien (ohne Streuung w_S = w_H),
   Nachrechnung veröffentlichter Kurven aus thrane2000 (Werte aus den Abbildungen, Toleranz nach
   Ablesegenauigkeit).
2. **S9c-2** allgemeines EHF mit der Phasenfunktion des Mediums (D11 B); Test: Gaußsche
   Phasenfunktion reproduziert S9c-1 auf 10⁻⁴, Whittle-Matérn gegen S9c-1 mit passendem θ_rms
   (Abweichung dokumentiert, keine Toleranz).
3. **S9c-3** Vergleich mit dem Hybrid (D12) an den Daten aus #34 (kein neuer Serverlauf);
   Testfall T10b.
4. **S9c-4** Querprobe BWPM im schwach streuenden Grenzfall (Roadmap): mittleres Signal über
   Realisierungen des Zufallsmediums gegen S9c-2, sobald die FPWPM-Korrektur (BWPM #37) vorliegt.
