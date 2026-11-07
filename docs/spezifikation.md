# Spezifikation (S0)

> Stand 2026-11-07. Grundlage für S1–S4; spätere Stufen ergänzen ihre Abschnitte.

## 1. Modell

Strahlungstransportgleichung (RTE) für die Strahldichte L(r, ŝ, t) in einem Medium mit
Absorptionskoeffizient μa, Streukoeffizient μs (1/mm), Phasenfunktion p(ŝ·ŝ′) und Brechzahl n:

    (n/c) ∂L/∂t + ŝ·∇L = −(μa + μs) L + μs ∫ p(ŝ·ŝ′) L(r, ŝ′, t) dΩ′ + Q

mit ∫p dΩ = 1, g = ⟨cos θ⟩, Albedo a = μs/(μa + μs), optische Dicke τ = (μa + μs) d. An Grenzen
mit Brechzahlsprung gilt Fresnel-Reflexion (unpolarisiert: Mittel aus |r_s|² und |r_p|², ab S6
vektoriell) und Snellius-Brechung.

## 2. Monte-Carlo-Verfahren (S1)

Photonenpakete mit Gewicht W (Leistungsanteil), wie MCML (wang1995), vektorisiert über viele
Pakete:

1. Start mit W = 1 in der Einfallsrichtung; spiegelnde Reflexion an der ersten Grenze wird
   deterministisch abgezogen (R_sp, Fresnel), das Paket tritt mit W = 1 − R_sp ein.
2. Schrittlänge s = −ln ξ/μt; trifft der Schritt eine Grenze, wird bis dorthin gelaufen und mit
   der Fresnel-Wahrscheinlichkeit reflektiert, sonst gebrochen (in die Nachbarschicht oder hinaus;
   hinaus = Beitrag zu R_d bzw. T_d mit dem aktuellen Gewicht); der Rest des Schritts wird in der
   neuen Schicht mit deren μt fortgesetzt.
3. Wechselwirkung: Absorption ΔW = W μa/μt (implizites Einfangen), Streuung nach der
   Phasenfunktion (Henyey-Greenstein: cos θ = (1 + g² − ((1 − g²)/(1 − g + 2gξ))²)/(2g), g = 0:
   cos θ = 2ξ − 1; Azimut gleichverteilt).
4. Russisches Roulette unter W < W_th (10⁻⁴): Überleben mit 1/m (m = 10), dann W ← mW.

Ergebnisse je eingestrahlter Leistung 1: R_sp, diffuse Reflexion R_d, Transmission T_d,
Absorption je Schicht, jeweils mit Standardfehler über Pakete; Energiebilanz
R_sp + R_d + T_d + A = 1 im Erwartungswert. Generator: `numpy.random.Generator(PCG64(seed))`.

## 3. Photonenprotokoll (ab S2)

Je austretendem Paket: Austrittsort (x, y), Richtung, Gewicht, Gesamtweglänge, Teilweglängen je
Schicht/Gewebe, Zahl der Streuungen, maximale Tiefe, Impulsübertrag je Gewebe (S8), Stokes-Vektor
(S6), bei Bedarf Phase (S9). Felder sind optional; der Kern führt nur die angeforderten mit.

## 4. Unabhängige Referenzen für S1

- **Ohne Streuung** (μs = 0): Beer-Lambert mit inkohärenten Mehrfachreflexionen an beiden Grenzen,
  geschlossen lösbar.
- **Halbunendliches Medium, isotrope Streuung, gleiche Brechzahl, senkrechter Einfall:**
  R = 1 − H(1)√(1 − a) mit Chandrasekhars H-Funktion (H(μ) = 1 + (a/2) μ H(μ) ∫₀¹ H(μ′)/(μ + μ′) dμ′,
  iterativ gelöst) – exakt, unabhängig vom Monte Carlo.
- **Adding-Doubling** (prahl1993) für Schichten mit HG-Streuung und Brechzahlsprüngen (S1b).
- Energieerhaltung.
