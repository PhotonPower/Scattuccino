# Regeln für Programmieragenten

Diese Datei liest jeder Agent zu Beginn jeder Sitzung. Sie gilt zusammen mit `docs/roadmap.md`,
`docs/konventionen.md`, `docs/abhaengigkeiten.md`, `docs/quellen.md` und `docs/testfaelle.md`.

## Projekt in einem Satz

Scattuccino simuliert die Lichtausbreitung in Gewebe und Gewebephantomen für alle gängigen optischen Verfahren der Medizintechnik: ein polarisierter, zeit- und weglängenaufgelöster Photonentransport (Monte Carlo) als gemeinsamer Kern, die Methoden als Auswertungen seines Photonenprotokolls, angebunden an Raytatouille und die Fusion (`PhotonPower/RTT-Fusion`).

## Bauen und testen

```bash
pip install -r requirements-dev.txt && pip install -e .
pytest
ruff check . && ruff format --check . && mypy
```

## Arbeitsregeln

1. **Ein Auftrag = ein Ziel einer Stufe** aus `docs/roadmap.md`. Stufen nicht vorziehen.
2. **Erst Test, dann Code.** Sollwerte aus analytischen Lösungen, unabhängigen Rechnungen oder
   veröffentlichten Messungen, nie aus dem geprüften Code. Tests und Sollwerte nur mit Freigabe
   des Maintainers ändern.
3. **Konventionen sind Gesetz** (`docs/konventionen.md`): die von Raytatouille und der Fusion.
4. **Abhängigkeiten minimal** (`docs/abhaengigkeiten.md`): NumPy und SciPy; sonst bevorzugt die
   Programme des Maintainers (PhotonPower). Fehlt dort etwas, wird ein Erweiterungsvorschlag unter
   `docs/vorschlaege/` geschrieben (dort nur lesend). Neue externe Abhängigkeiten nur mit Freigabe.
5. **Fremder Code:** GPL- oder unklar lizenzierter Code wird nicht übernommen, nur als Vergleich
   über veröffentlichte Ergebnisse genutzt.
6. **Formeln mit Quelle** (`docs/quellen.md`); Gleichungsnummern nur aus dem Volltext. Status
   jeder Quelle festhalten (Volltext gelesen / nur zitiert).
7. **Determinismus:** Zufall nur mit explizitem Seed; Monte-Carlo-Ergebnisse mit Standardfehler.
8. **Keine Geheimnisse** in Dateien, Commits, Remote-URLs.
9. **Unklarheit = Frage**: Issue mit Label `frage`; Abweichung von einer Quelle: Label `paper`.

## Sprache und Stil

- Code, Kommentare, Docstrings, Commit-Messages englisch; Projektdokumentation deutsch.
- ruff (Zeilenlänge 100), mypy strict; Commits im Imperativ, ein Thema pro Commit.

## Definition of Done (PR-Checkliste)

- [ ] CI grün
- [ ] Referenztests mit unabhängigem Sollwert
- [ ] `docs/testfaelle.md` aktualisiert (Status, Messwerte)
- [ ] neue Quellen in `docs/quellen.md`
- [ ] Docstrings mit Einheiten, Vorzeichen und Verweis
- [ ] `CHANGELOG.md` unter „Unreleased“
