**Headline:** The read-back test compares a planet's irregular moons to a rounding, the way it already compares every other float.

### Fixed
- `test_planets_moons_belts_and_comets_are_read_back_as_the_dicts_they_were_written_from` compared the stored irregular set (a dict of moons) for exact equality, while it compares every other float to 1e-12. An irregular moon is stored in km and g/cm3 and read back in AU and kg/m3, so a distance or density comes back one rounding off (59.96918180368469 AU read back as 59.969181803684684). The test now recurses into the set and holds each value to the same 1e-12, and takes its comet from a system asked to hold one (the physics gives few stars comets, so 60 seeds may hold none); the store and the generator are unchanged, and a missing or extra key still fails.

generation-output: unchanged
