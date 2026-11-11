**Headline:** a belt and every bound comet of a star now pull on its planets and moons as tide sources, so the effect ledger lists them beside the other bodies that raise a tide.

### Changed

- `physics/tide_sources.py` (GEN.1464, tide side): `belt_source` makes a belt one source of kind `belt`: its mass cut into the eight rings of `secular_chain.belt_rings` (the rings nearer the receiver than `BELT_ALPHA_LIMIT` left out, as the secular pull does), summed so the accelerations of the rings add, with the rings' torques and heats added apart and the tide going round at the mass-weighted rate of the rings. `comet_source` makes a bound comet a source on its own orbit. `planet_sources` and `moon_sources` list both among the siblings of a planet (a moon takes them at its planet's place); `register` gives each body the belts of its own host's list, and `sources_of`/`moon_sources` take the planet's list when the caller has it.
- `generation/system.py`: the tide pass works host by host, so a belt or a comet pulls on the planets and moons of its own star only. The circle means in `tide_sources` take arrays of radii (same values as before, faster).
- Registry: TS5 rewritten (belt and comet sources, reusing BR1 and BR2), TS7 added (dwarf planets and belt members are in their belt's mass; the tide on a belt member's own moons is not built).
- Census, 60 systems (554 planets and moons; 30 have a belt, 56 have comets): the largest belt pull on any body is 1.8e-4 of the largest source on it (under the ledger floor of 1e-3, so it is folded into each body's summary row) and the largest comet pull is 4e-13; the spring and neap tide range of a body moves by at most 3.4e-9 of itself, none by 1e-6. The effect ledger gains the belts and comets as tide sources, and the stored tide ranges differ in their last digits.

generation-affects: fill
