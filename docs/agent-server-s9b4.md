# Auftrag für einen Agenten: Abnahme des Hybrid-OCT – Unabhängigkeit von z_h (S9b-4)

> Für einen Programmieragenten auf dem Rechner des Maintainers. Es gelten `AGENTS.md`, die
> Entscheidungen D8–D10 (`docs/entwurf-s9b-hybrid.md`) und die Testfälle T17a/T17b.

## Ziel

Das Abnahmekriterium aus D9: Das modengekoppelte, kohärenzgegatete OCT-Signal in Tiefen jenseits
der Übergabe darf nicht davon abhängen, wo die Wellenoptik an den Monte Carlo übergibt. Lokal
(2·10⁵ Pakete) schwanken die tiefen Bins mit dem Next-Event-Schätzer noch um den Faktor 2; mit
10⁷ Paketen und mehreren Seeds soll die Streuung klein genug sein, um z_h = 0, ½/μs, 1/μs, 2/μs zu
vergleichen.

## Regeln

- Branch `s9b4-nee`; nur ein neues Skript unter `tools/` und der Bericht; keine Änderungen an
  Code, Tests, Toleranzen. RTT-Fusion (mit `phasespace.field_to_packets(periodic=…)`, PR #19 dort)
  und BWPM mit CUDA-Kern installiert.
- Läufe über 30 min vorher abschätzen.

## Läufe

`scattuccino.hybrid.hybrid_oct(medium, n0=1.38, mua=1.0, wavelength_um=0.8, beam, z_h, z, l_c,
packets, seed, nx, dx, max_path, sigma_husimi, backend="cuda")`, Next-Event-Schätzer (Standard).

| Größe | Wert |
| --- | --- |
| Medium | `WhittleMatern(0.0572, 0.5e-3)` (Born μs = 200 /mm, g = 0,935; 1/μs = 5 µm) |
| Strahl | `Beam(w_f = 2e-3, z_f = 0.012)` (Fokus 12 µm tief) |
| Gitter | `nx = 256`, `dx = 1.5e-4` (38 µm seitlich) |
| z_h | 0; 0,0025; 0,005; 0,010 mm |
| A-Scan | `z = np.linspace(0, 0.06, 31)`, `l_c = 0.004` |
| Pakete | 10⁷ je Lauf, `max_path = 0.15`, `sigma_husimi = 0.6e-3` |
| Seeds | 1, 2, 3 je z_h (Mittelwert und Standardfehler je Tiefe) |

## Auswertung

Für Tiefen z ≥ 15 µm (hinter allen Übergaben): Verhältnis zu z_h = 0 mit Standardfehler aus den
Seeds. Abnahme, wenn alle Verhältnisse innerhalb 3σ bei 1 liegen oder um höchstens 15 % abweichen
(Born-Parameter unten gegen eine BWPM-Realisierung oben: bekannte Abweichung einige Prozent,
T16). Zusätzlich die Laufzeit je Lauf und die Zahl der Streuereignisse.

## Bericht

PR „S9b-4: hybrid OCT acceptance runs“ (Basis `main`): Tabelle Tiefe × z_h (Mittel ± Fehler,
Verhältnis), Umgebung, Skript, Befund.
