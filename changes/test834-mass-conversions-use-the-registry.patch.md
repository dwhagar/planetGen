**Headline:** Two hand-written mass conversions use the unit registry again

### Fixed

- TEST.834: `test_mass_conversions` was red on main: the GEN.1493 slice 3 late supply (`late_supply.eject_share`, a planet's mass over the star's mass in kilograms) and the GEN.1814 melt sulfur (`melt_sulfur`, a body's mass in Earth masses) each multiplied or divided by a mass constant by hand. Both now call `unit_registry.converter`; the numbers are the same to the last bit (200,000 draws compared), so generation output is unchanged (35 tests pass across `test_mass_conversions`, `test_late_supply`, `test_melt_sulfur`).

generation-output: unchanged
