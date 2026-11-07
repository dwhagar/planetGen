**Headline:** The planets a star's birth cluster pulls off their orbits now leave the system and become rogue planets that say their birth cluster unbound them.

### Changed
- A star that belongs to a cluster (`stars.cluster_id`, GEN.359) is built with the cluster's history: its density and its relative speed from the end of the disk to the cluster's age, worked from the stored mass, half-mass radius and age (`cluster_unbinding.history_of`; the radius grew from a birth radius as `open_cluster_evolution` has it, the speed is the virial one of a Plummer sphere). Each planet of the star's architecture takes a draw; a planet the cluster unbinds (`unbinding_epoch_myr`, the rate of [R4184] eqs. 36 and 37 in the pieces of the history where the orbit is soft against the encounters, so a tight orbit is never lost) is taken out of the system, recorded as a `stripped` loss of the architecture (`planet_architecture.strip`), and made again as a rogue by the same chain as a thrown-out planet.
- The rogue leaves at the speed that just unbinds it (`sqrt(2) - 1` of its circular speed, added to its orbital velocity), keeps its moons out to the Domingos fraction of its own Hill radius, and its note reads "Unbound by its birth cluster from the ... system ...". An `EjectionEvent` has a `cause` ("thrown" or "cluster").
- The same history is given to a scattered star's system built on its own (the uncharted-systems route), so both ways of building it agree. A system asked for by planet count or slot list is left as asked.
- Stand-ins CU1 (soft pieces only) and CU4 to CU7 (section 6.113). The wiring guard line for `cluster_unbinding` is gone.
- Generated stars of clusters change their planets (golden epoch).

generation-output: changed
