**Headline:** The materials table now holds the clays, amphiboles, feldspathoids, pyroxenes and oxides the recipes and the norm name, each with its density, hardness, crystal system and colour.

### Added
- **45 silicate and oxide minerals and 3 groups in the materials table (GEN.1355, part 1 of the missing minerals).** Talc, brucite, the smectites (montmorillonite, saponite, nontronite, beidellite), the chlorites (clinochlore, chamosite), the serpentine kin (antigorite, lizardite, chrysotile), illite, palygorskite, sepiolite, the amphiboles (tremolite, hornblende), epidote, nepheline, leucite, hypersthene, augite, pigeonite, jadeite, analcime, wollastonite, the CAI phases (hibonite, grossmanite, melilite, akermanite, gehlenite), the iron oxides and hydroxides (maghemite, lepidocrocite, akaganeite, ferrihydrite, schwertmannite, ulvospinel, magnesiochromite) and the silica and lime phases (cristobalite, tridymite, opal, lime, portlandite, lonsdaleite), plus the `smectite`, `chlorite` and `amphibole` groups. Each value was read from the mineral's Webmineral page on 2026-11-07; wollastonite is 2.84 g/cm3 (the research file's recalled 2.9 is replaced), hornblende reads the magnesiohornblende page, lepidocrocite 4 is confirmed.
- **Mohs hardness, crystal system and colour on a materials row** (`mohs` property in the table, `Material.crystal_system`, `Material.colour`); the loader refuses a Mohs value off the scale 1 to 10 or a crystal system it does not know.
- **A flag on every ground row without the three visible-band albedos** (the note says `ground_albedo` raises `MissingAlbedo` for it), and a test that fails a ground row with neither the bands nor the flag.
- The 45 minerals and 3 groups are tagged vocabulary-only (TEST.685) except nepheline, leucite, hypersthene and wollastonite, which the CIPW norm already names; the norm module joins the reachability scan.

`generation-output: unchanged`
`TODO items resolved: GEN.1355`
