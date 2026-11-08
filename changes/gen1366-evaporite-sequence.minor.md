**Headline:** The salts on a dried ocean floor, and in the ocean itself, now come from the chemistry of the body's own water instead of a fixed halite, gypsum and calcite mix.

### Added
- **The evaporite sequence (GEN.1366, `physics/evaporite_sequence.py`).** A body's brine is built from the cations of the crust it weathers (the melt of its mantle, with sodium and potassium put back because the stored rock carries no alkalis) and the anions of its ocean class (neutral, soda, chloride brine, acid sulfate). Its ions go to salts in the order of the chemical divides: jarosite and alunite (acid water), calcite, the alkaline path (magnesium carbonate, then nahcolite, trona, thermonatrite or natrite), the calcium path (gypsum or anhydrite, epsomite, hexahydrite or kieserite, mirabilite or thenardite, calcium chloride), perchlorate in thin oxidizing air, then halite, carnallite, sylvite, bischofite and antarcticite. Gypsum or anhydrite follows Hardie's measured equilibrium of temperature and water activity (R5143), so a halite-saturated brine makes anhydrite above about 18 degrees C. The other hydrate boundaries and the ion budget are registered stand-ins (section 6.212, rows EV1 to EV11).
- Two table rows the sequence needs: natron and antarcticite (Webmineral pages read on 2026-11-07).
- The Hardie source on the gypsum and anhydrite equilibrium is in the references registry as R5143.

### Changed
- The evaporite ground and the ocean's dissolved salt are read from the sequence (`surface_composition.dried_salts`, `ocean_salts`): the old fixed mix (four tenths halite, three tenths gypsum, three tenths calcite) and the halite-only ocean salt are gone. A body with no silicate rock still gets halite. The areal ground mix, the draws and the albedo rows are unchanged; the stored mass composition changes, so the next epoch is needed.
- Twenty-one rows of the materials table that GEN.1355 and GEN.1357 tagged vocabulary-only are now read by the sequence and lose the tag (TEST.685).

`generation-output: changed`
`TODO items resolved: GEN.1366`
