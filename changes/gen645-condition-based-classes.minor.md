**Headline:** A planet or moon's class can now be worked out from its finished conditions (air, temperature, water, size and gravity), and the class pages show the rule.

### Added
- **The classifier (GEN.645, GEN.647 first slice).** `planetgen.physics.planet_classes.classify(body)` returns a class letter from a planet or moon's stored mass, gravity, mean surface temperature, surface pressure and gas mixture, water cover and share, and starlight variation. The rules are the ordered rows of `physics/data/planet_classes.json` (a one-sentence reading and a `basis` each); the last row of each kind takes the rest, so every body gets one letter. All 20 letters stay, H, P, O, J, M, I, C, D and N are marked core, and R, U, W, X, Y and Z are kept free for new classes. Design: `docs/design/condition-based-planet-classes.md`.
- **The shadow sweep.** `python src/tests/class_shadow_cli.py` compares the derived letter with the drawn one (98.3 percent of 32,487 bodies from 2,000 systems keep their letter) and prints the confusion table. Generation still draws a class; nothing it stores changes.

### Changed
- **The planet class pages give the rule.** `/classes/planet/<letter>` shows the class's name, sentence, the rule in plain words, which classes are checked before it and whether it is a core class, in place of the table of radius ranges, zones and pigments the zone draw used.
