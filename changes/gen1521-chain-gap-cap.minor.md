**Headline:** the chain-spacing rule no longer pushes a planet so far from its neighbour that the gap reads as a belt, which brings the belt count of Train 24 back down.

### Changed

- `clearance_ledger.chain_distance_au` (GEN.1521) leaves the pairwise rule when the nearest distance that holds the chain's effective spacing would open a gap of `CHAIN_MAX_GAP_K` (40) mutual Hill radii or more, the width at which the belt rule reads a gap as a belt (`belt_origin.GAP_SPACING_K`, registry CH5, section 6.272). The effective spacing is the mean less half the spread, so beside one wide gap the next planet was pushed to gaps of 43 to 60, the share of neighbour pairs at 40 or more rose from 7 to 16 percent, and the extra pairs were drawn as belts (G2V: 1.21 belts per system on main, 1.73 on Train 24).
- Effect, on the census's own 100 systems per star type: belts per system 1.43, 1.43, 1.25 and 1.89 (A0V, G2V, K5V, F5V) before, 0.83, 1.16, 0.91 and 0.94 after, planets per system unchanged to within 0.1 to 0.3 of the draw; the cap leaves the pairwise rule in place of a push that opened a gap of 40 or more between semi-major axes.

generation-output: changed
generation-affects: plan
