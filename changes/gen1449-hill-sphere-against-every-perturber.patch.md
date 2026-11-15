**Headline:** A planet's moons now stay inside the part of its Hill sphere that its companion stars and sibling planets leave it, not the host star's alone.

### Changed

- GEN.1449 (slice 1): a moon's outer stability limit, the irregular moons' zone and the Hill escape of the moon history were worked out against the planet's host star alone. They now read the Hill sphere shrunk by the mean tide of the host's companion stars and the other planets of the star (`physics/hill_perturbers.py`, reusing `tide_sources.planet_sources` of GEN.1446, `orbits.calculate_hill_sphere` and the Domingos fit of `clearance_ledger`). A planet with no companion star and no sibling gets exactly the sphere it had, so its moons are unchanged; a planet whose orbit crosses a massive sibling's, or one near a close companion star, loses its far moons. Stand-ins HP1 and HP2 are in GEN.725 (section 6.299). Stored values change for a planet with a companion or sibling: its moons' distances, by a hair in a spaced system (under 0.1% in 200 sampled systems), and the draws after them.

generation-output: changed
generation-affects: fill
