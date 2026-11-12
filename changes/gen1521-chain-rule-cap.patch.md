**Headline:** the chain-spacing rule no longer flings a planet to many times its distance, and a system whose planets were left too close by a late eccentricity change is spaced again.

### Changed

- `clearance_ledger.chain_distance_au` (GEN.1521) moves a planet at most `CHAIN_MAX_PUSH` (2) times its pairwise distance out for the chain rule, looking for the nearest distance that holds the effective spacing in steps of `CHAIN_SCAN_STEP` (registry CH4). A window beside a giant can reach the effective spacing only far out or never (the gap in mutual Hill radii saturates with the distance), and the doubling search that stood here bracketed a root 9 times out for one such window (M3V seed 1317: 2.2 AU became 18.5 AU after an eccentricity changed by a part in a thousand). A window that no such distance closes keeps the pairwise rule, as one that no distance closes already did.
- After the settle pass, a star's planets whose spacing the final check would refuse are spaced again with `validation.space_orbits`; a star whose planets pass is not touched. An eccentricity the oddities force after the spacing pass (resonant chain W5, `finish_orbit_oddities`) had closed a belt's gap by 0.13 percent (M3V seed 1273).
- Effect: on 120 systems (40 each of M3V, G2V, K5V; seeds 4000 to 4039), 63 of the 67 with a chain-rule evaluation had it bind; in 23 of those 67 an evaluation would have pushed a planet more than twice out (212 of 16,735 binding evaluations, 117 of them five times out or more, up to 633 AU). Those systems come out different (the planets stay at the pairwise distances); the rest keep their distances to within 1e-12, because the search finds the same nearest distance. The two failures of the M3V validation test (400 seeded systems, seeds 1000 to 1399) validate.

generation-output: changed
generation-affects: plan
