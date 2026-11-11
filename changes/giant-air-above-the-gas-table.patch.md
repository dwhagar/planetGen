**Headline:** A giant lit by a very hot star no longer stops the build of its system; its air is worked out at the hottest temperature the gas table reaches.

### Changed

- GEN.1514 (bug): `mixture_giant.giant_mixture` solved the layer's gas at the layer's own temperature, and the gas table ends at 20,000 K. A giant or brown dwarf near an O star (zeta Ophiuchi's example, seed 152 of a seeded sweep, 27,560 K at the quoted layer) found no gas to make and the system failed with "no species can be made of these elements at this temperature". The layer is now equilibrated at no more than `mixture_giant.CHEMISTRY_CEILING_K` (the table's top, 20,000 K), as a colder one already was at `CHEMISTRY_FLOOR_K`; the temperature drawn is kept everywhere else (the quench level, the record). Ionization is not modelled, so the stand-in is the neutral and molecular gas at 20,000 K, tier "from needed output", tweak L0, until the table is extended.
- This was the cause of `test_example_star_age_never_exceeds_lifespan[zeta_ophiuchi_system.json]` failing on a rare draw (the test builds the example three times with unseeded draws, and a build that stops fails it). Test: `tests/test_atmosphere_mixture_giant.py` (20,000, 27,560 and 50,000 K).

generation-affects: none
