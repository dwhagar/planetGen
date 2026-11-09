**Headline:** The stellar wind a planet's magnetic field holds off is now the wind of every star of its system, and a moon inside its planet's magnetosphere holds off the planet's field instead of the wind (GEN.1452).

### Changed

- `magnetism.wind_pressure_at` sums the host's wind (the number of before) and each companion's wind at 1 AU times the mean `1/d^2` between the body and the companion (the geometry GEN.1441 added for the light). The magnetopause, the field the star's wind induces in a rocky planet and the radiation belts (and so the dose of their moons) read the sum. A lone star gives the same numbers to the last digit.
- A moon whose orbit lies inside its parent's magnetopause has its own magnetopause set against the pressure of the parent's field at its orbit (`magnetism.apply_parent_magnetosphere`, called from `generate_induced_field`), not the wind's. Jupiter and Ganymede give 2.4 Ganymede radii.
- Stand-in section 6.240 (WN1 to WN4) in the registry; design note updated in `habitability/activity-magnetism-radiation-hydrosphere.md`.
- Census (120 seeded systems, binary and wide-binary mix, `/tmp` script not kept): 164 planets store a magnetopause; with a companion attached, 2 of them move by more than 0.1 percent, by -0.17 to -0.28 percent (the companions of the draw sit far out, so their wind adds little); no planet of a lone star moves by any amount (0 changes). Of 81 moons, 2 have their own field with a field-holding parent and none lies inside the parent's magnetopause, so the moon rule changes none of them yet; it is covered by the Ganymede test. The golden digests of the phenomenon case move because of this change; the sector cases already moved with GEN.1441.

generation-output: changed
