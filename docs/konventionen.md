# Konventionen

Verbindlich; dieselben wie Raytatouille, Fieldtracing, BWPM und RTT-Fusion.

- **Einheiten:** Längen in mm, Wellenlängen in µm (Vakuum), Zeiten in ps, Koeffizienten μa, μs
  in 1/mm, Winkel in rad.
- **Koordinaten:** rechtshändig; Gewebeoberfläche bei z = 0, Gewebe bei z > 0, Licht von z < 0
  (Normalfall). Abweichende Anordnungen über die Pose im Fusion-Modell.
- **Felder** (wo Wellenoptik eingeht): ∝ exp(i(k·r − ωt)), Fresnel Convention A.
- **Polarisation:** Stokes-Vektor (I, Q, U, V) und Mueller-Matrizen in den Konventionen von
  `rtt-polar` (Raytatouille); Bezugsebenen werden im Protokoll mitgeführt.
- **Phasenfunktion:** normiert auf ∫p dΩ = 1; g = ⟨cos θ⟩; μs′ = μs(1 − g).
- **Gewicht:** Photonenpakete tragen ein Gewicht (Leistungsanteil); Ergebnisse je eingestrahlter
  Leistung 1 und mit Standardfehler.
- **Zufall:** Generator und Seed in jedem Ergebnis.
- **Ergebnisse:** Format `raytatouille-result` (ADR 0023), wie RTT-Fusion.
