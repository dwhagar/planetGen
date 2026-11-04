**Headline:** The Magellanic Clouds hold about as many giant molecular clouds as their gas would: the Large Cloud's count falls from 847 to 338 and the Small Cloud's from 270 to 216.

### Changed
- The nebula field's gas tracer for a satellite is its young stars' density (GEN.528), and a Cloud's young stars fill a ball out to six half-mass radii, so the field kept far more clouds than the gas holds. `ComponentField.gas_density` now scales a named satellite's tracer by `tuning.SATELLITE_NEBULA_GAS_SCALE` (Large Cloud 0.52, Small Cloud 0.85; any other satellite keeps 1). The class M clouds the field keeps come to 338 for the Large Cloud and 216 for the Small, the Milky Way's 8,107 times each Cloud's molecular mass over 1.2e9 solar masses: 5e7 for the Large Cloud (272 clouds catalogued above the survey limit) and 3.2e7 for the Small.
- Three sources are on the reference list: [R5096] the Large Cloud's CO catalogue, [R3102] the Milky Way's, [R5097] the Small Cloud's gas. Stand-in rows MG1 and MG2 (section 6.171) hold the two scales with their tier and tweak level.

### Added
- `tests/test_gen367_cloud_gas.py`: the sourced numbers, each Cloud's expected class M clouds against the target within 15 percent, the unscaled count (847) and the solved scale, and a satellite that is not named keeping the young-star tracer.

generation-output: changed
TODO items resolved: GEN.367 (the epoch and golden ride train 11; only galaxies with the satellite galaxies switched on are touched)
