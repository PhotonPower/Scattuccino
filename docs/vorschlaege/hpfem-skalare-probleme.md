# Vorschlag für hp-FEM: skalare H1-Probleme für Diffusion, Akustik und Wärme

> Entwurf aus Scattuccino (D5); in hp-FEM nur lesend. Zur Entscheidung durch den Maintainer.

hp-FEM hat bereits H1-Räume (für die konische Gitterrechnung) neben Nédélec-Elementen. Drei
Aufgaben der biomedizinischen Optik sind skalare H1-Probleme derselben Form

    −∇·(a ∇u) + c u = f   (a, c komplex, ortsabhängig)

| Anwendung | a | c | Bemerkung |
| --- | --- | --- | --- |
| Diffusion im Frequenzbereich (NIRS, DOT, SFDI) | D = 1/(3(μa + μs′)) | μa + iω/v | Robin-Rand (extrapolierte Grenze) |
| Photoakustik, Ausbreitung im Frequenzbereich | 1/ρ | −ω²/(ρc²) | PML wie im optischen Fall |
| Wärme (Pennes, stationär oder Zeitschritt) | κ | ω_b ρ_b c_b (+ ρc/Δt) | Quelle = absorbierte Energie aus Scattuccino |

**Gewünscht:** skalarer Pfad mit Materialfeldern auf Elementen, Robin- und PML-Rand, Quelle als
Volumenfeld (Datenport), Ausgabe von u und Fluss; Python-Anbindung wie die vorhandene. Abnahme:
analytische Lösungen (Punktquelle im unendlichen Medium, Halbraum mit extrapolierter Grenze,
Kugelwelle). Vorteil: keine fremden FEM-Pakete (FEniCS, NIRFAST, TOAST++, k-Wave).
