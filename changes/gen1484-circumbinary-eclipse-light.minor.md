**Headline:** A planet that circles a close pair of stars is slightly cooler: its light is cut by the share the two stars hide of each other, 0.1 to 1 percent.

### Changed
- **Eclipses of a close pair take light from its planets (GEN.1484).** `eclipses.light_factor_of` gives a planet of a close pair the share of its light it keeps after the two stars eclipse each other, from the pair's separation, period and plane and the planet's own height above it (the eclipse functions of GEN.1462, registry rows EC1 to EC4). It is stored as `light_factor` (`planets`, `moons`, `rogue_planets`, migration 0241; NULL for every other body) and applied to the light in the air and temperature chain, a giant's light and radius, the mixture's equilibrium temperature, the thermal tide and the energy flux of the habitability score. A median circumbinary planet loses 3e-4 of its light and the most lose 1.9 percent (56 planets in 400 systems), which is 0.1 to 0.5 percent in temperature. A lone star, a moon and a rogue planet read as before. New epoch, `--affects fill`.

generation-output: changed
generation-affects: fill
