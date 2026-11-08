**Headline:** Each of the 84 old asteroid component names now has a host mineral in a primitive body's mix or is tagged vocabulary-only, and a belt records the names its minerals can carry.

### Added
- **Asteroid vocabulary hosts (GEN.1369, `physics/asteroid_vocabulary.py`).** Every name of the old `tuning.ASTEROID_COMPONENTS` list is either hosted (the platinum-group elements and gold in kamacite and taenite, the osmium-iridium alloy in the inclusion phases, titanium in ilmenite and perovskite, the polycyclic aromatic hydrocarbons in the organic matter, the major elements in the silicates, metal and sulfides) or vocabulary-only with the setting that forms it and the item that builds it (magmatic sulfide ores, carbonatites and pegmatites, supergene zones, evaporites, igneous crusts). The host table is a recalled stand-in (section 6.214, ACV1). A test enforces both lists and that every key the mixes use is a row id and not an alias.

### Changed
- A belt's stored makeup detail gains `hosted`: the vocabulary names its mixes can carry, each with its largest host row. No mass is drawn for them, so the shares, porosities and densities are unchanged; the stored belt detail moves, so the next epoch is needed.

`generation-output: changed`
`TODO items resolved: GEN.1369`
