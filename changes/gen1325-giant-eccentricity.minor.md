**Headline:** Giant planets now draw their eccentricity from the Beta law fitted to the radial-velocity giants (mean 0.22) instead of the Kepler small-planet laws, so generated giants average 0.18 to 0.20 and not 0.08 (GEN.1325).

### Changed

- `orbit_elements.draw_planet_eccentricity(multiple, e_max, giant=False)`: a giant (`body_type` `g`), alone or with other planets, takes Beta(0.867, 3.03) cut at its periapsis limit; every other planet keeps the laws of OE1 to OE3. `planets.generate_orbit_shape` passes the flag. Stand-in OE16 in section 6.12.
- Measured on 900 systems (94 giants): mean eccentricity 0.18 (median 0.10, 90th percentile 0.44), 0.20 beyond 1 AU; before, 0.064 over 300 systems (50 giants). The tides still circularise the close giants and the placer still keeps the orbits apart.
- The orbit of every giant changes, so the planets laid out after it change too (epoch).
