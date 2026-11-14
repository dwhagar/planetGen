**Headline:** The kind of ice a body is made of now sets its strength: a comet of methane clathrate is stronger than one of water ice, and a large cold icy moon holds ice II to VI at its center.

### Changed

- The roundness engine takes `ice_factor`, the strength of the kind of ice (`ice_types`, GEN.1177) over water ice Ih (`roundness.yield_strength_pa`, `round_radius_m`, `round_mass_kg`, `roundness_ratio`, new `roundness.ice_strength_pa`). `body_strength` names the kind: a comet, a belt member and the members of a field from their makeup (`ice_type_of_makeup`), a planet or a moon from the phase of water ice at its center (`ice_type_of_body`, none for a body whose shape was set above water's triple point). `body_status.status`, `roundness` and `round_mass_kg` take it; the belt page and the field members pass it.
- The status word and the pages call a body round by `roundness.is_round` (the one definition of a ratio of 1 or more), not by comparing with 1 by hand.
- Wiring register (GEN.1684): the four roundness constants (RN1, RN2, RN4, RN11) and `compaction.taper_weight` are marked by design with their reasons; `roundness.is_round` and the `ice_types` module are read and their rows dropped.
- Measured on 40 seeds (12,831 bodies and members, before and after): no belt member changes its status word and no planet, moon or comet changes between round and not round; 4 of 573 planets, 4 of 480 moons and 2 of 18 comets change their stored roundness ratio and yield strength (strength 0.1 to 7.2 times what it was, the median 7.2 for the cold icy bodies in ice II to VI).

generation-output: changed
generation-affects: fill
