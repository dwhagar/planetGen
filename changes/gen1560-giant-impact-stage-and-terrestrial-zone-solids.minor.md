**Headline:** The rocky planets of a Sun-like star's temperate zone are now Earth-sized instead of Mars-sized: their zone holds more solids and its oligarchs merge in a giant-impact stage, so the habitable zone has about half a planet per star where it had nearly one, and more of them are Earth-mass.

### Changed

- GEN.1560: from 0.45 AU to the snow line the disk's solids the oligarchs are built from are 4 times (`PLANET_ROCKY_ZONE_SOLIDS_BOOST`, stand-in PA16) what they were; then the oligarchs of that zone collide, the closest pair in mutual Hill radii first, until every pair is 25 Hill radii apart (`GIANT_IMPACT_SPACING`, `planet_architecture.giant_impacts`, stand-in PA17), mass kept to `COLLISION_MASS_RETENTION`. The inner system and everything past the snow line are as they were.
- Round a Sun-like star's architecture (1,000 drawn): 0.57 planets of 0.5 to 1.5 Earth radii in the habitable zone (0.95 to 1.7 AU) against 0.93 before (the survey, R0572: 0.37 to 0.60), and 0.63 planets of 0.75 to 2 g at 0.7 to 2.5 AU against 0.46. In 90 whole F, G and K systems the temperate rocky planet has a median of 1.4 Earth masses against 0.5.
- Every FGK system with a terrestrial zone changes; M dwarfs, whose snow lines sit inside 0.45 AU, are unchanged.

generation-output: changed
generation-affects: fill
