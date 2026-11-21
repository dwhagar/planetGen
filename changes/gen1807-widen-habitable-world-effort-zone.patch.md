**Headline:** The zone that habitable-world effort goes to (the Class M tier dials and the request for a habitable world) is now twice as wide as the habitable zone, half its width more on each side.

### Changed

- GEN.1807: `orbits.candidate_zone` widens a star's habitable zone by `orbits.CANDIDATE_ZONE_WIDENING` (0.5, registry row PA20, tier T9, tweak L3) of its width on each side; for the Sun 0.95 to 1.68 AU becomes 0.59 to 2.04 AU. The tier dials (`abundance_tier.candidate`: water class, envelope, quiet interior, carbon cycle, desert share) and the request for a habitable world (`StarSystem._generate_planets`, the ecosphere check of the request) read that zone. The habitable zone itself, a body's zone letter (`h`, `e`, `c`) and the pages keep the published edges, and the class of a body still comes last from its own conditions. Science (the default) draws nothing for the dials and no request, so its output does not move.
- The test fakes of the dials carry the zone and the distance a real planet carries.

generation-output: changed
