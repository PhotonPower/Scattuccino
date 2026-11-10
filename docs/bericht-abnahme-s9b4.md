# Bericht: Abnahme S9b-4, Hybrid-OCT (Unabhängigkeit von z_h)

Auftrag: Issue 29. Skript `tools/abnahme_s9b4.py`, Daten `tools/out/abnahme_s9b4/`.

## Ergebnis: nicht bestanden

Von 69 Vergleichen (23 Tiefen ab 16 µm × 3 Werte von z_h) verfehlen 4 das Kriterium (Verhältnis
innerhalb 3σ bei 1 **oder** höchstens 15 % Abweichung):

| Tiefe µm | z_h µm | Verhältnis zu z_h = 0 |
| --- | --- | --- |
| 16 | 2,5 | 1,164 ± 0,037 |
| 18 | 10 | 1,316 ± 0,077 |
| 20 | 10 | 1,289 ± 0,084 |
| 50 | 10 | 0,744 ± 0,072 |

Bewertung (über das Kriterium hinaus, nur beschreibend):

- Unmittelbar hinter den Übergaben (16–22 µm) liegt das Signal mit z_h > 0 systematisch über dem
  mit z_h = 0: mittleres Verhältnis über diese vier Tiefen 1,10 (z_h = 2,5 µm), 1,20 (5 µm),
  1,25 (10 µm); der Effekt wächst mit z_h. Tiefer (24–60 µm) liegen die mittleren Verhältnisse bei
  1,15 / 1,08 / 1,07, dort dominiert das Rauschen (Fehler 5–40 %).
- Innerhalb von 15 % liegen 15, 16 und 13 von 23 Tiefen (z_h = 2,5 / 5 / 10 µm).
- Die Fehler sind Standardfehler aus nur 3 Seeds; bei so wenigen Seeds ist 3σ ein schwaches
  Kriterium (t-Verteilung mit 2 Freiheitsgraden), und die Fehler selbst streuen stark. Die Abnahme
  ist deshalb statistisch nicht sauber entscheidbar: die tiefen Bins sind nach 10⁷ Paketen und
  3 Seeds noch zu verrauscht, die flachen (16–22 µm) zeigen dagegen einen Trend jenseits des
  Rauschens. Mehr Seeds oder Pakete bräuchte eine Freigabe; ich habe nichts darüber hinaus
  gerechnet. Ursache des Trends habe ich nicht untersucht.
- `power_handoff` (Vorwärtsleistung an der Übergabe, einfallende Leistung 1) liegt zwischen 0,989
  und 1,021, also nicht streng ≤ 1; Ursache nicht untersucht.

## Umgebung

RTX 6000 Ada, nvcc 12.9, Python 3.13.5, BWPM mit CUDA-Kern (`True True`), RTT-Fusion `51bba21`
(Merge PR #19), Scattuccino auf `main` `52cc941`. Tests `test_s9b2_microstructure`,
`test_s9b3_hybrid`, `test_s9b4_next_event`: 10 bestanden, keine Skips. Beam w_f = 2 µm, Fokus
12 µm, n₀ = 1,38, μa = 1 /mm, 10⁷ Pakete, `max_path` = 0,15 mm, σ_Husimi = 0,6 µm, 256² bei
dx = 0,15 µm.

## Kalibrierung (Schritt 0, 256² × 134 Ebenen, Seed 11, Laufzeit 64 s)

| Größe | Wert |
| --- | --- |
| μs Born | 200,12 /mm |
| μs gefittet (BWPM) | 190,73 /mm (Verhältnis 0,953) |
| g vorwärts Born | 0,9462 |
| g vorwärts BWPM (`cal.g_forward`) | 0,9506 |
| g der gemessenen Phasenfunktion (`scattering[1].g`, ganze Kugel) | 0,9398 |
| max \|total − 1\| | 0,0032 |

Rohdaten (kohärentes Feld, Tiefen) in `kalibrierung.json`.

## Laufzeiten

Pilot: 165 s für 10⁶ Pakete (Hochrechnung 6,5 h für alle Läufe, unter der Grenze von 12 h). Die 14
Läufe mit 10⁷ Paketen brauchten 613–764 s (z_h = 0), 769–989 s (2,5 und 5 µm), 1029–1159 s (10 µm);
Einzelwerte in den `npz`-Dateien (`runtime_s`), Summe etwa 3,3 h. Satz B: 674 s (z_h = 0) und 819 s (5 µm).

## Signal je Tiefe (Mittel ± Standardfehler, 3 Seeds, Satz A)

| z µm | z_h=0 µm | z_h=2.5 µm | z_h=5 µm | z_h=10 µm |
| --- | --- | --- | --- | --- |
| 16 | 0.000648 ± 1.4e-05 | 0.0007543 ± 1.8e-05 | 0.0007906 ± 4.9e-05 | 0.0007692 ± 5.6e-05 |
| 18 | 0.000512 ± 1.4e-05 | 0.0005768 ± 1.1e-05 | 0.0006183 ± 4.2e-05 | 0.0006739 ± 3.5e-05 |
| 20 | 0.0003976 ± 1.1e-05 | 0.0004328 ± 5.8e-06 | 0.0004885 ± 3.2e-05 | 0.0005124 ± 3e-05 |
| 22 | 0.0003229 ± 1.2e-05 | 0.0003348 ± 7.4e-06 | 0.0003665 ± 1.6e-05 | 0.0003885 ± 2.4e-05 |
| 24 | 0.0002604 ± 1.2e-05 | 0.0002689 ± 1.4e-05 | 0.0002773 ± 1.7e-05 | 0.0003015 ± 1.4e-05 |
| 26 | 0.0002026 ± 1.4e-05 | 0.0002152 ± 1.1e-05 | 0.0002259 ± 1.1e-05 | 0.0002401 ± 9.2e-06 |
| 28 | 0.0001674 ± 8.9e-06 | 0.0001724 ± 1.2e-05 | 0.0001816 ± 9.6e-06 | 0.0001831 ± 5e-06 |
| 30 | 0.0001374 ± 8.8e-06 | 0.0001436 ± 1.8e-05 | 0.0001447 ± 8.6e-06 | 0.000138 ± 4.2e-06 |
| 32 | 0.0001079 ± 6.5e-06 | 0.0001195 ± 1.8e-05 | 0.0001172 ± 9.2e-06 | 0.0001154 ± 5.5e-06 |
| 34 | 8.702e-05 ± 7.3e-06 | 9.83e-05 ± 1e-05 | 9.306e-05 ± 6.4e-06 | 9.667e-05 ± 4.8e-06 |
| 36 | 7.949e-05 ± 2.3e-06 | 8.547e-05 ± 2.9e-06 | 7.524e-05 ± 4.4e-06 | 7.553e-05 ± 7.5e-06 |
| 38 | 7.032e-05 ± 9.6e-07 | 7.62e-05 ± 2.9e-06 | 6.385e-05 ± 4.3e-06 | 6.294e-05 ± 7e-06 |
| 40 | 5.596e-05 ± 2.9e-06 | 6.43e-05 ± 4.2e-06 | 6.002e-05 ± 1.5e-06 | 5.667e-05 ± 3.4e-06 |
| 42 | 4.449e-05 ± 2.7e-06 | 5.426e-05 ± 1.9e-06 | 5.157e-05 ± 2.5e-07 | 4.785e-05 ± 6.8e-07 |
| 44 | 3.756e-05 ± 3.3e-06 | 4.881e-05 ± 5.8e-06 | 3.946e-05 ± 3e-06 | 4.251e-05 ± 9.8e-07 |
| 46 | 3.307e-05 ± 2.4e-06 | 4.384e-05 ± 7.2e-06 | 3.417e-05 ± 2.6e-06 | 4.062e-05 ± 3.1e-06 |
| 48 | 3.274e-05 ± 2.3e-06 | 3.67e-05 ± 5e-06 | 3.553e-05 ± 5.6e-07 | 3.233e-05 ± 4.6e-06 |
| 50 | 3.328e-05 ± 1.1e-06 | 3.14e-05 ± 3.6e-06 | 3.437e-05 ± 2.4e-06 | 2.477e-05 ± 2.3e-06 |
| 52 | 2.642e-05 ± 1.2e-06 | 2.973e-05 ± 9.5e-07 | 3.073e-05 ± 2.8e-06 | 2.432e-05 ± 1.8e-06 |
| 54 | 2.147e-05 ± 2.9e-06 | 2.598e-05 ± 2e-06 | 2.851e-05 ± 3.1e-06 | 2.288e-05 ± 1.8e-06 |
| 56 | 1.799e-05 ± 2.7e-06 | 2.141e-05 ± 3.7e-06 | 2.288e-05 ± 8e-07 | 2.052e-05 ± 2.9e-06 |
| 58 | 1.497e-05 ± 5.9e-07 | 2.117e-05 ± 6.4e-06 | 1.513e-05 ± 1.3e-06 | 1.94e-05 ± 1.6e-06 |
| 60 | 1.3e-05 ± 2e-06 | 1.774e-05 ± 3.3e-06 | 1.299e-05 ± 1.2e-06 | 1.696e-05 ± 1.6e-06 |

## Verhältnis zu z_h = 0 (Fehler aus den Seeds; ✗ = Kriterium verfehlt)

| z µm | z_h=2.5 µm | z_h=5 µm | z_h=10 µm |
| --- | --- | --- | --- |
| 16 | 1.164 ± 0.037 ✗ | 1.220 ± 0.080 | 1.187 ± 0.090 |
| 18 | 1.126 ± 0.036 | 1.208 ± 0.088 | 1.316 ± 0.077 ✗ |
| 20 | 1.088 ± 0.033 | 1.229 ± 0.086 | 1.289 ± 0.084 ✗ |
| 22 | 1.037 ± 0.044 | 1.135 ± 0.064 | 1.203 ± 0.085 |
| 24 | 1.032 ± 0.071 | 1.065 ± 0.080 | 1.158 ± 0.075 |
| 26 | 1.062 ± 0.094 | 1.115 ± 0.096 | 1.185 ± 0.096 |
| 28 | 1.030 ± 0.088 | 1.085 ± 0.081 | 1.094 ± 0.065 |
| 30 | 1.045 ± 0.145 | 1.053 ± 0.092 | 1.004 ± 0.071 |
| 32 | 1.107 ± 0.184 | 1.086 ± 0.108 | 1.069 ± 0.082 |
| 34 | 1.130 ± 0.150 | 1.069 ± 0.117 | 1.111 ± 0.109 |
| 36 | 1.075 ± 0.048 | 0.947 ± 0.062 | 0.950 ± 0.098 |
| 38 | 1.084 ± 0.043 | 0.908 ± 0.062 | 0.895 ± 0.101 |
| 40 | 1.149 ± 0.095 | 1.073 ± 0.061 | 1.013 ± 0.080 |
| 42 | 1.220 ± 0.086 | 1.159 ± 0.070 | 1.076 ± 0.067 |
| 44 | 1.300 ± 0.193 | 1.050 ± 0.123 | 1.132 ± 0.104 |
| 46 | 1.326 ± 0.238 | 1.033 ± 0.109 | 1.228 ± 0.130 |
| 48 | 1.121 ± 0.172 | 1.085 ± 0.079 | 0.988 ± 0.157 |
| 50 | 0.944 ± 0.112 | 1.033 ± 0.081 | 0.744 ± 0.072 ✗ |
| 52 | 1.125 ± 0.063 | 1.163 ± 0.120 | 0.920 ± 0.081 |
| 54 | 1.210 ± 0.187 | 1.328 ± 0.231 | 1.066 ± 0.167 |
| 56 | 1.190 ± 0.272 | 1.272 ± 0.195 | 1.141 ± 0.234 |
| 58 | 1.414 ± 0.434 | 1.010 ± 0.096 | 1.296 ± 0.120 |
| 60 | 1.365 ± 0.330 | 0.999 ± 0.178 | 1.305 ± 0.232 |

## Satz B (Born) gegen Satz A (gemessene Streuung), Seed 1

Verhältnis B/A je Tiefe 16…60 µm (gleicher Seed, aber verschiedene Streuung, also verschiedene
Pfade; die tiefen Bins sind rauschdominiert):

- z_h = 0: 1,007 0,972 0,950 0,909 1,049 1,157 1,015 0,946 0,951 1,081 0,890 0,785 0,790 1,000
  1,280 1,067 0,822 0,915 1,409 1,955 2,300 1,796 1,236
- z_h = 5 µm: 0,978 0,944 0,950 0,980 0,951 0,920 0,967 0,975 0,946 1,113 1,480 1,549 1,195 1,089
  1,357 1,398 0,937 0,656 0,697 0,810 0,917 1,089 0,997

Bis 34 µm liegt B/A innerhalb etwa 10 % (meist 0,95–1,05, einzelne 1,16), die gemessene Streuung
ändert das flache Signal also nur wenig; ab 36 µm sind die Abweichungen größer als der
Seed-Fehler der Tabelle oben (einzelner Seed). Eine Aussage zur Verbesserung durch `scattering`
ist mit einem Seed nicht möglich.
