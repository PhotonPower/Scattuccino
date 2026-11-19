# Entwurf S9c: erweitertes Huygens-Fresnel (EHF) als analytische OCT-Referenz

> Stand 2026-11-19. D11 (Option C) und D12 freigegeben 2026-11-19. D13: thrane2000 und levitz2004
> im Volltext gelesen; lutomirski1971 liegt vor (für S9c-2); schmitt1997 fehlt noch.

## Kontext

Der Hybrid (S9b) wird bisher nur in sich geprüft: Unabhängigkeit von der Übergabetiefe z_h
(T17, Abnahme #34) und Einfachstreuung im engen Fokus gegen eine geschlossene Form (T18). Es fehlt
eine **unabhängige Referenz für das OCT-Signal mit Fokus und Mehrfachstreuung**. Das liefert das
EHF-Modell von Thrane, Yura und Andersen (thrane2000), aufbauend auf der EHF-Ausbreitung durch
Zufallsmedien (lutomirski1971) und dem Modell von Schmitt und Knüttel (schmitt1997):

- **Annahmen** (thrane2000, Abschnitt 2): Streuung vorwärts unter kleinen Winkeln (paraxial,
  bis etwa 0,5 rad), Volumenrückstreuung vernachlässigt, reflektiert wird an einer Grenzfläche in
  der Tiefe z (diffus oder spiegelnd), Absorption vernachlässigt, keine Polarisation. Hin- und
  Rückweg statistisch unabhängig; Phase normalverteilt mit Varianz μs z, Gaußsche
  Phasenkorrelation.
- **Ergebnis:** ⟨i²(z)⟩ = ⟨i²⟩₀ Ψ(z) mit ⟨i²⟩₀ = α² P_R P_S σ_b/(π w_H²) und

  Ψ(z) = e^(−2μs z) + 4 e^(−μs z)(1 − e^(−μs z))/(1 + w_S²/w_H²) + (1 − e^(−μs z))² w_H²/w_S²

  (thrane2000 Gl. 21 integriert; levitz2004 Gl. 2). Die Terme sind ungestreut hin und zurück,
  Kreuzterm und gestreut hin und zurück. w_H und w_S sind die **1/e-Radien der Intensität** in der
  Tiefe z ohne und mit Streuung (Gl. 13, 14; ABCD-Elemente A = 1, B = d + z/n); w_S enthält die
  seitliche Kohärenzlänge ρ₀(z) = √(3/(μs z)) · λ/(π θ_rms) · (n d + z)/z (Gl. A6, A7; Gl. 24 für
  den Fokus auf der Grenzfläche, d = f − z/n). Spiegelnde Reflexion: Gl. 25, 26.
- **θ_rms:** Halbwertsbreite bei e⁻¹ einer Gauß-Anpassung an die Vorwärtskeule der
  Phasenfunktion; das Modell verwendet g_eff = cos θ_rms ≥ g (levitz2004, Abschnitt 2.2; die
  Beispiele in thrane2000 nutzen θ_rms = arccos g).

**Druckfehler:** thrane2000 Gl. 22 hat beim Kreuzterm den Faktor 2. Die Integration von Gl. 21 gibt
4, wie levitz2004 Gl. 2 (Issue mit Label `paper`). Für die Kurven in thrane2000 Abb. 2 und 3 ist
der Unterschied ohne Belang, weil dort w_H²/w_S² ≪ 1 ist.

Der erste Term ist im Grenzfall μs z ≪ 1 die Einfachstreuung im Fokus, die T18 bereits als
geschlossene Form prüft. Das verbindet S9c mit vorhandenen Tests.

**Offener Punkt für S9c-2:** In Gl. A1 steht die Vakuum-Wellenzahl k mit der reduzierten
Länge B(s) = s/n. Daraus folgt eine Aufweitung des Halos w_S² − w_H² = θ_rms² μs z³/(3n²),
unabhängig von d. Ob θ_rms der Winkel im Gewebe ist (dann wäre der Faktor 1/n² zu prüfen), klärt
die eigene Herleitung in S9c-2 für die Geometrie des Hybrids (Mode im Medium, Brechzahl angepasst).

## D11: Welche Form des EHF (freigegeben: C)

| Option | Inhalt | Bewertung |
| --- | --- | --- |
| A geschlossene Form (thrane2000) | Gaußsche Kleinwinkelnäherung der Phasenfunktion (θ_rms), Ψ(z) als Formel | direkt mit veröffentlichten Kurven vergleichbar; θ_rms ist für die Whittle-Matérn-Phasenfunktion (Spitze plus breiter Fuß) nur eine Näherung |
| B allgemeines EHF | gegenseitige Kohärenzfunktion mit der Phasenstrukturfunktion des Mediums, berechnet aus der tatsächlichen Phasenfunktion (Kleinwinkelteil), Integrale numerisch | gleiche Annahmen wie A, aber ohne θ_rms; passt zum Medium der Hybridläufe |
| **C = A + B** | A als Prüfstein gegen die veröffentlichten Kurven, B als Referenz für den Hybrid | B muss im Grenzfall der Gaußschen Phasenfunktion A exakt reproduzieren; das prüft die numerische Umsetzung unabhängig |

## D12: Vergleichsgröße Hybrid ↔ EHF (freigegeben)

Der Hybrid liefert das in die Detektionsmode gekoppelte, kohärenzgefensterte Signal je mm,
normiert auf die eingestrahlte Leistung; das EHF liefert ⟨i²(z)⟩ bis auf die Systemkonstante.

- Verglichen wird der **Verlauf** S(z)/S_einfach(z), also das Signal geteilt durch die geschlossene
  Form der Einfachstreuung (T18) bei derselben Tiefe, gegen Ψ(z)/e^(−2μs z). Die Systemkonstante,
  der Rückstreukoeffizient und die Fokusfunktion kürzen sich heraus; übrig bleibt der Beitrag der
  Mehrfachstreuung.
- Zusätzlich ein absoluter Vergleich bei kleinem μs z (dort ist Ψ ≈ e^(−2μs z) und der Hybrid ist
  durch T18 absolut geprüft).

## D13: Volltexte und Abgrenzung

- **Volltexte:** thrane2000 und levitz2004 gelesen (S9c-1). lutomirski1971 für S9c-2;
  schmitt1997 (Rückstreuung im Volumen statt an einer Grenzfläche) vor S9c-3.
- **Grenzen des EHF**, die die Vergleiche bestimmen: nur Kleinwinkelstreuung (das Medium aus #29
  hat g ≈ 0,95; geeignet), Reflexion an einer Grenzfläche statt Rückstreuung im Volumen, keine
  Mehrfachrückstreuung (diffuser Hintergrund fehlt in großer Tiefe). Abweichungen des Hybrids
  **jenseits** dieser Annahmen sind kein Fehler, sondern werden als Befund dokumentiert.
- **Speckle** (Statistik der einzelnen Realisierung, zweiter Teil von S9c in der Roadmap) folgt
  nach S9c-3 in einem eigenen Entwurf.

## Stufen

1. **S9c-1** geschlossene Form nach thrane2000 (D11 A) in `scattuccino.ehf`; Tests: Gl. 22 gegen
   die Quadratur von Gl. 21, ρ₀ gegen die Quadratur von Gl. A1, Grenzfälle (μs → 0: Ψ → 1; ohne
   seitliche Kohärenz: Ψ → e^(−2μs z)), Nachrechnung von thrane2000 Abb. 2 und 3 (Werte aus den
   Abbildungen, Toleranz nach Ablesegenauigkeit).
2. **S9c-2** allgemeines EHF mit der Phasenfunktion des Mediums (D11 B) für die Geometrie des
   Hybrids; Test: Gaußsche Phasenfunktion reproduziert S9c-1, Whittle-Matérn gegen S9c-1 mit
   passendem θ_rms (Abweichung dokumentiert, keine Toleranz).
3. **S9c-3** Vergleich mit dem Hybrid (D12) an den Daten aus #34 (kein neuer Serverlauf);
   Testfall T10b.
4. **S9c-4** Querprobe BWPM im schwach streuenden Grenzfall (Roadmap): mittleres Signal über
   Realisierungen des Zufallsmediums gegen S9c-2, sobald die FPWPM-Korrektur (BWPM #37) vorliegt.
