**Headline:** two bug-sweep findings settled with no change to generated output: the carbon and nitrogen inventories that follow the water are the documented rule and are now pinned by a test, and a stored row that is a super-puff now reads as one.

### Changed

- GEN.1598: `habitability_world.inventories` credits carbon and nitrogen with at least the water's inventory because the habitability-index design says so (volatiles arrive with the water); it is not a copy slip. The docstring says so and a water-world case joins `test_carbon_and_nitrogen_follow_the_water_or_the_air`. No stored value changes.
- GEN.1601: `planet_type.planet_type` given a column-named mapping now returns super-puff when the mapping's stored `planet_type` says so and its envelope is transitional (a stored row has no star or radius law to judge by); a stored type never overrides the envelope share. No caller in the generator passes a mapping, so no stored value changes.
