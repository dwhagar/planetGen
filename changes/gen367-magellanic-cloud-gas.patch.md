**Headline:** The Magellanic Clouds hold about as many giant molecular clouds as their gas would: the Large Cloud's count falls from 1,085 to 338 and the Small Cloud's from 340 to 216.

### Changed
- The nebula field's gas tracer for a satellite is its young stars' density (GEN.528), and a Cloud's young stars fill a ball out to six half-mass radii, so the field kept far more clouds than the gas holds. `ComponentField.gas_density` now scales a named satellite's tracer by `tuning.SATELLITE_NEBULA_GAS_SCALE` (Large Cloud 0.435, Small Cloud 0.724; any other satellite keeps 1). The class M clouds the field keeps come to 338 for the Large Cloud and 216 for the Small, the Milky Way's 8,107 times each Cloud's molecular mass over 1.2e9 solar masses: 5e7 for the Large Cloud (272 clouds catalogued above the survey limit) and 3.2e7 for the Small.
- Three sources are on the reference list: [R5123] the Large Cloud's CO catalogue, [R3102] the Milky Way's, [R5124] the Small Cloud's gas. Stand-in rows MG1 and MG2 (section 6.171) hold the two scales with their tier and tweak level.
- The unscaled counts are 1,085 and 340 on the main that has GEN.1246 (the thin disk's infall mass term, `tuning.DISK_INFALL_MASS_TERM_PER_GY`); with that term at 0, as before it, they are 847 and 268. The Milky Way reference of the nebula field moved by design with GEN.1246, and the scales were solved again on the new reference.

### Added
- `tests/test_gen367_cloud_gas.py`: the sourced numbers, each Cloud's expected class M clouds against the target within 15 percent, the unscaled count (1,085) and the solved scale, and a satellite that is not named keeping the young-star tracer.

generation-output: changed
TODO items resolved: GEN.367 (the epoch and golden ride train 11; only galaxies with the satellite galaxies switched on are touched)
