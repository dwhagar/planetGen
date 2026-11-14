**Headline:** Belts are born the spacing apart that the spacing pass keeps, so far fewer belts are nudged outward after the fact.

### Changed

- GEN.1673: the belt-moved-outward rate (`belt.spacing_push`) was 0.17 per system, not because of the chain rule (GEN.1521; it was 0.16 with the rule off) but because of how belts are placed. The unfinished tail and the outer disk both meet at the growth limit, so the outer disk started 0 AU from the tail and was moved 0.05 AU out in 70 percent of the pushes; and a belt's inner edge was worked out from the star's mass at formation while the planet inside it keeps its clearance from the star's present mass (white dwarfs), so it started inside that clearance. `belt_origin.build` now starts a stretch `MIN_ASTEROID_BELT_SEPARATION` past the stretch before it when they meet, and a belt keeps out of a planet by the planet's own clearance (`Planet.min_orbit_distance`). On 400 seeded systems the rate is 0.018 per system, the belts per system are unchanged (0.87).
- Systems that have a tail belt and an outer-disk belt, or a belt beside a white-dwarf planet, are drawn the same but their belt edges move by up to a few percent.

generation-output: changed
