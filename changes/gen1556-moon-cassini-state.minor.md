**Headline:** A moon locked to its planet now holds its axis at the tilt its own shape and the pull on its orbit give it (Europa a tenth of a degree, Io almost none, the Moon a few degrees), instead of a tilt drawn from a made-up spread.

### Changed

- GEN.1556 (part 2 of GEN.1536): `planets.settle_moon_tilts` gives each locked regular moon the tilt of its first Cassini state, `i |s| / (V + |s|)`, from the regression of its node (the planet's bulge `(3/2) n J2p (Rp/a)^2` plus the star's `(3/4) n_p^2 / n`; `cassini_state.moon_regression_rad_yr`) and its own precession constant (the hydrostatic figure for an icy or hot moon, a fossil one for a cold rocky moon, as for planets). The axis is turned in the planet's equator frame and back. Checked on the Galilean moons (0.01 to 0.14 degrees against Io under 0.04 and Europa 0.1) and the Moon (3.7 to 4.7 against the 6.7 of its state 2). Registry OB13 and OB14; SP22 now replaced (burn-down line set to replaced).
- Locked moons' tilts change, so the golden moves; no schema change.

generation-output: changed
generation-affects: fill
