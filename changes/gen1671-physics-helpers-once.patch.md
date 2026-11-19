**Headline:** The physics helpers that were copied across modules (rule-file checks, number checks, 3-vector arithmetic, a body's spin axis and orbit normal) are written once; the generator's output is unchanged.

### Changed
- GEN.1671: `elements.check_rule_entry` (new: the tier, tweak and basis check of one rule entry) and `elements.load_rule_numbers` replace the copies in `trojans`, `sulfuric_vapor`, `sulfur_cloud_set`, `life_supply` and `life_air`; `elements.positive_number` and `elements.finite_number` replace the copies in `comet_population` and `materials`; `vec3` (with a new `sub`) replaces the copies in `galaxy/navigation`, `galaxy/system_position`, `physics/collision_outcome`, `physics/lagrange_geometry` and `db/orbital_events`; `spin.spin_axis_of` and `spin.orbit_normal_of` replace the copies in `planets`, `system`, `kozai_bodies` and `eclipses`. A rule file's refusal now words its message the same way in every module, and `test_comet_population` matches "must be a positive number".

generation-output: unchanged
generation-affects: none
