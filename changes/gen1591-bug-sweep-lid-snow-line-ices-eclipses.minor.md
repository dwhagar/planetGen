**Headline:** Five small bugs from the bug sweep are fixed: a remade planet's air is solved under the same lid and surface as a first build, a planet's tectonic regime is the body's own throw, a migrated planet's snow-line flag follows where it formed, a moon with ice draws its other ices, and a fully eclipsed planet keeps no light.

### Changed

- GEN.1591: `calculate_atmospheric_conditions` clears the lid regime and the surface temperature the pass before left, so a remade body solves its air under the stagnant lid and default surface of a first build, and `generate_lid` takes its throw from a stream of the body's own (no draw of the run's stream, the same throw whenever the body is made again).
- GEN.1596: `born_beyond_snow_line` follows where the body formed (`formed_beyond_snow_line`), and is NULL for a star with no snow line instead of failing.
- GEN.1600: the secular exchange compares the spacing problems a system already had without the distances their messages quote, so a hair's move does not make an old problem look new.
- GEN.1603: a moon that holds its origin and has ice draws its other ices from the comet volatiles as an icy planet does (it stored the median species); a body given both a mass and a radius draws the iron branch for its rock and core nickel as one drawn from its star does.
- GEN.1606: a giant's ices follow where it formed, not where it orbits now; a light factor of exactly 0 keeps no light (`eclipses.kept_share`), where `or 1.0` turned it into all of it.
- GEN.1595 was checked and changes nothing: 25 of 25 rocky planets that came out airless stayed airless on a full replay with their system's late accretion (a comment in `remake_air` says so).

generation-output: changed
generation-affects: fill
