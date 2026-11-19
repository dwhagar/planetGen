**Headline:** A giant can now hold a whole planet at its L4 or L5, grown from the solids of its 1:1 zone, on the giant's orbit 60 degrees ahead of it or behind it.

### Added
- GEN.1151 slice A (lock item, migration 0270 on this branch, 0272 at the train: the revision id must equal the schema version in order, so the branch cannot carry 0272 until the earlier ones are on main): `generation/trojan_planet.py`. For each giant that grew where it orbits, the solids of its two tadpole zones (`trojan_bodies.zone_solids_earth`, the nebula the planets and belts were made from), a share of them (TB3, 0.02) over the share of the growth time the gas disk allowed (TB4, a ramp from 0.3 to 3 Myr) and the room the mass rule leaves make a ceiling; a ceiling over the lightest planet draws a mass log-uniform to it, and the giant's neighbours and companion stars decide whether the body stays (the survival of the swarm, GEN.1455). No rate is set.
- The planet takes the giant's plane, node, periapsis and eccentricity and stands 60 degrees from it in mean anomaly (`Planet(trojan_of=...)`, `build.trojan_of`), and is stored with `planets.lagrange_point` and `planets.lagrange_host_id`; the loader links it to its giant.

### Changed
- TB3 `FORMATION_EFFICIENCY` 0.1 to 0.02 (tweak L2) and TB4 from a flat 1 Myr minimum to a ramp (0.3, 3.0) Myr (tweak L1), because the measure (`reports/gen1151/`) found the flat step on the median G dwarf disk life and a third of the giants holding a planet-mass trojan at 0.1; new rows TB5 (zone half width), TB6 (azimuth share), TB7 (migration share).

generation-output: changed
generation-affects: fill
