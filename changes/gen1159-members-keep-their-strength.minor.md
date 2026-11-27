**Headline:** Every body of a belt or an asteroid field now keeps its roundness ratio, yield strength and binding energy.

### Added
- A stored belt member and a stored field member carry three more numbers (`roundness_ratio`, `yield_strength_pa`, `binding_energy_j`, schema 0281), made by the one strength module (`body_strength.member_values`) from the member's own mass, size, density, ice, metal and porosity, with the kind of its ice (clathrate or water) and the cold strength a comet has. The strength draw is keyed on the member's own numbers, so the run's stream does not move and the same body has the same strength each time it is made again (a body drawn on demand, GEN.1674, included).
- A planet keeps whether its tilt wanders chaotically (`body_history.tilt_chaotic`, the quantity oddities W6 and W61 read) and the moons its tide brought inside its Roche limit (the new `fallen_moons` table: mass and the years after the planet formed that each fell, which its rings are made from), and a restored planet gets both back (`body_history.attach`, `attach_fallen`), where it had neither after a load (GEN.1874). A planet saved before has no such rows and reads as it did.
- A belt page body and a field member's page data show them in a `strength` block; a member stored before this release says `null` for each until the belt is filled again.

generation-output: changed
generation-affects: fill
