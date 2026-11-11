**Headline:** The time the disk's gas takes to carry the vapor of an ice line to a giant now comes from the disk's own drift speed (its temperature, the star's mass and the turbulence), not a fixed 0.05 Myr per AU.

### Changed

- `disk_drift.lag_myr` was a fixed 0.05 Myr per AU for every star. It is now the distance over the speed a steady accretion disk carries its gas inward, `v = 3 nu / (2 r)` with `nu = alpha c_s^2 / Omega` (new `viscosity_at_1_au_au2_per_myr` and `drift_speed_au_per_myr`, stand-in DE5): 10.5 AU per Myr for the Sun (0.095 Myr per AU), 8.8 for half a solar mass, 13 for 1.4 solar masses. The old 0.05 was half the steady speed's lag for the Sun, so lags roughly double for Sun-like stars and the gas at a giant's place stays the disk's own gas for longer (GEN.1517).
- Giants of a Sun-like star that gather their gas at 3 Myr no longer see the CO2 and argon pebble ice from the lines 30 AU away (their weights go from about 0.12 to 1); the water weight of a 300 K giant moves from 0.80 to 0.87, and the carbon-to-oxygen ratio of a 200 K giant from 0.94 to 0.84. Giants of dwarf stars (all gather at the disk's end) and of heavier stars hardly change.

generation-output: changed
generation-affects: fill
