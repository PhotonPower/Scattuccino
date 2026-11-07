# Entwurf S9b: Übergang Wellenoptik → Photonentransport (BWPM → Monte Carlo)

> Entwurf zur Freigabe (Entscheidungen D8–D10), Stand 2026-11-07. Betrifft auch die Fusion
> (Ports, ADR-Entwurf Fusion, Punkt „Darstellung P“).

## Kontext

OCT, konfokale und Multiphotonen-Mikroskopie fokussieren Licht ins Gewebe. In den ersten freien
Weglängen bestimmen Wellenoptik und Mikrostruktur das Signal (Fokusform, Aberration durch
Brechzahlschwankungen, Kohärenz, Speckle); tiefer ist das Licht diffus und der Photonentransport
(S1–S9a) die passende Beschreibung. Dazwischen liegt eine Übergabeebene z = z_h:

- **oben** (0 < z < z_h): BWPM mit einer **Realisierung der Mikrostruktur** (Brechzahlfeld
  n(r) = n₀ + δn(r) auf dem BWPM-Gitter), angeregt vom fokussierten Feld aus Fieldtracing
  (Port G in RTT-Fusion);
- **unten** (z > z_h): Monte Carlo mit den **statistischen** Größen desselben Gewebes (μs, g bzw.
  Phasenfunktion, μa, n).

Drei Fragen sind zu entscheiden: wie aus dem Feld bei z_h Photonenpakete werden (D8), wo z_h liegt
und wie oben und unten zueinander passen (D9), und wie das zurückkommende Licht detektiert wird
(D10).

## D8: Feld → Pakete (Phasenraum)

Ein Paket braucht Ort **und** Richtung; das Feld E(x, y) bei z_h legt die gemeinsame Verteilung
nur über eine Phasenraumdarstellung fest.

| Option | Verteilung von (x, y, k⊥) | Bewertung |
| --- | --- | --- |
| A Strahlnäherung | Ort ∝ \|E\|², Richtung = lokaler Poynting-Vektor (Phasengradient) | einfach; falsch in Fokusnähe und für Speckle (Phase nicht glatt) |
| B getrennt | Ort ∝ \|E\|², Richtung aus dem Gesamtspektrum \|Ẽ(k⊥)\|² | verliert die Korrelation Ort–Richtung (ein Fokus wird zur Kugelwelle aus jedem Punkt) |
| **C Husimi** | Wigner-Funktion von E, geglättet mit einem Gauß-Fenster minimaler Unschärfe (Breite σ, σ_k = 1/(2σ)) – immer ≥ 0 | **Vorschlag**: positiv, energieerhaltend (Integral = Leistung), richtige Grenzfälle (Strahl, Fokus, diffuses Feld); σ wählbar zwischen Auflösung in Ort und Richtung |

Umsetzung C: Q(x, y, k⊥) = |⟨φ_{x,k}|E⟩|² mit Gauß-Wellenpaketen φ, berechnet als
Kurzzeit-Fourier-Transformation von E auf dem BWPM-Gitter (FFT je Fensterposition, für
572²-Gitter einige Sekunden); Pakete per Rejection oder Alias-Tabelle aus Q, Gewicht = Leistung/N,
Richtung k⊥/(n k₀) (evaneszente Anteile verworfen und als Verlust gebucht), Polarisation aus den
Feldkomponenten des Fensters (Stokes-Vektor, für S6).

**Prüfungen:** Leistung erhalten; Fokus ohne Streuung: Pakete aus C, ballistisch weitergeführt,
treffen den Fokus wie die BWPM-Fortsetzung des Feldes (Strehl-Ort, Breite innerhalb σ);
diffuses Zufallsfeld: Richtungsverteilung lambertsch.

## D9: Übergabetiefe und Konsistenz oben/unten

- z_h vom Nutzer, Vorschlag **z_h ≈ 1–2 / μs** (die ballistische Leistung ist auf e⁻¹…e⁻² gefallen,
  Mehrfachstreuung beginnt). Ergebnisse dürfen von z_h nur innerhalb ihrer Toleranz abhängen –
  **Abnahmekriterium**: z_h = 1/μs und 2/μs liefern dasselbe OCT-Signal und dieselbe Reflexion.
- **Konsistenz:** Das Zufallsfeld δn(r) oben (Gaußsches Zufallsfeld mit Whittle-Matérn-
  Korrelation, Rogers/Backman) wird so gewählt, dass seine Streuung die Parameter unten trifft.
  Prüfung in BWPM selbst: Extinktion der ballistischen Leistung → μs, Winkelspektrum des
  gestreuten Feldes → Phasenfunktion und g. Erst wenn beide auf wenige Prozent passen, wird
  gekoppelt. (Born-Näherung gibt die Startwerte: μs und p(θ) aus dem Spektrum von δn.)
- Rückwärts gestreutes Licht **oberhalb** von z_h: BWPM liefert die Reflexion des oberen Bereichs
  selbst (bidirektionale BWPM, B3); sie wird zum Monte-Carlo-Ergebnis addiert. Licht, das von
  unten wieder durch z_h nach oben läuft, durchquert den oberen Bereich als Monte Carlo mit den
  statistischen Größen (keine Rückkopplung ins Feld; die Mikrostruktur oben wirkt auf dem Rückweg
  nur statistisch).

## D10: Detektion des zurückkommenden Lichts

| Option | Vorgehen | Bewertung |
| --- | --- | --- |
| A geometrisch | Pakete mit Ort/Winkel im Fänger (wie S9a) | ignoriert die Modenanpassung der Faser; überschätzt Mehrfachstreuung |
| **B Modenanpassung im Phasenraum** | Gewicht jedes austretenden Pakets × Husimi-Dichte der **Detektionsmode** (Rückwärtsfortsetzung der Fasermode durch Objektiv und oberen Bereich, also derselbe Strahl wie die Beleuchtung bei Reziprozität) am Austrittsort/-richtung | **Vorschlag**: entspricht der heterodynen Effizienz im erweiterten Huygens-Fresnel-Modell (S9c-Referenz); ohne Streuung exakt die Kopplung in die Mode |
| C vollständig rückwärts | zurückkommende Pakete in Feld zurückverwandeln und mit BWPM/Fieldtracing propagieren | nicht eindeutig (Pakete haben keine Phase); verworfen |

Kohärenz-Gating wie S9a über den optischen Weg; oben addiert die BWPM den Weg ihres Bereichs
als Phase, unten der Monte Carlo als OPL.

## Folgen für die Fusion

- Neuer Übergang **G → P** (Feld → Phasenraum-Pakete) und **P → Kopplungsgewicht** in eine Mode;
  beides als Funktionen in `rtt-fusion` (Python-Prototyp zuerst), nicht in BWPM oder Scattuccino.
- Port P (radiometrisch, M9) bekommt damit seine erste konkrete Form: Pakete mit Ort, Richtung,
  Gewicht, optischem Weg, Stokes.

## Stufen nach Freigabe

1. **S9b-1** Husimi-Übergang (D8) mit Prüfungen; ohne Gewebe gegen Fokus aus Fieldtracing/BWPM.
2. **S9b-2** Zufallsmedium in BWPM, Kalibrierung μs, g gegen Born und Monte Carlo (D9).
3. **S9b-3** Hybrid-OCT mit Modenanpassung (D10); Unabhängigkeit von z_h.
4. **S9c** erweitertes Huygens-Fresnel als analytische Referenz für 3.
