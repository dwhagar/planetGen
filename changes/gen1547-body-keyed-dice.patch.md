**Headline:** The dice of every planet and moon end stage (air, haze, aerosols, strength, oddities, mineral rows, isotopes, giant settle, Kozai) are now keyed to the body's birth, so remaking or moving a body never re-rolls them.

### Changed

- GEN.1547: `draw.body_stream_of` (kind, mass, the distance the body formed at, its star's mass) replaces the streams that held the body's current distance, radius, surface temperature, air pressure or the star's age: the first-pass air (`calculate_atmospheric_conditions` no longer draws from the run's stream), `refinish-rocky`, `settle-giant`, `settle-giant-air`, `late-accretion-end`, haze, aerosols, body strength, the oddity key, the ore, sediment, carbide and hydrothermal rows, the isotope record and the Kozai draw. `test_end_stage_replay.py` (Lane 07's prebuild) now passes without its xfail marks: it replays an end stage on a body moved by a hair in distance, age and temperature and requires the same streams. Four tests pinned to seed draws were repinned (W19, W20, argon, the escape sample) and the aerosol test now asserts the new rule.

generation-output: changed
generation-affects: fill
