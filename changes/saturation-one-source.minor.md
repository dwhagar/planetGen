**Headline:** The frozen nitrogen air of a rogue planet and the boiling point of water now come from the one vapour-pressure function.

### Changed
- **One saturation function (GEN.632).** `rogue_surface.nitrogen_vapor_pressure_pa` and its constants are replaced by `frozen_nitrogen_pressure_pa`, which reads `physics/saturation.py` (nitrogen frost pressure at 40 K is 6.6 Pa, as before, but follows the real curve up to the triple point and the liquid curve above it, capped at the critical pressure); `hydrosphere.water_boiling_k` finds the temperature at which the water curve reaches the pressure, and `WATER_TRIPLE_POINT_PA`, `WATER_BOILING_K`, `WATER_CRITICAL_K` and `WATER_VAPORIZATION_J_MOL` are gone from `hydrosphere`. A test fails if one of the old names comes back. Stored values that change: a rogue planet's surface nitrogen pressure and the temperature at which a world's water boils.

generation-output: changed
