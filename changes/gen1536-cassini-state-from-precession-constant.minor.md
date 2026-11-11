**Headline:** A planet locked to its star or in a spin-orbit resonance now holds its axis at the tilt its own shape and its neighbours' pull give it, a few arcminutes to a fraction of a degree for a cold rocky world (Mercury's figure gives 2.4 arcmin against the 2.04 measured), instead of a tilt drawn from a made-up spread.

### Changed

- GEN.1536: new `physics/cassini_state.py`. The precession constant of a locked or resonant rocky planet comes from its figure (`J2`, `C22`): the hydrostatic figure for an icy body or a rocky one with a heat flow of 0.5 W/m2 or more, and for a cold rocky body that has despun a fossil figure 22.7 to 94.2 times the hydrostatic `J2` (the Moon to Mercury, Table 2 of [R5228]), drawn log-uniform from a stream of the planet's own. The tilt is the first Cassini state against the nodal modes of its orbit (the siblings' and a companion star's), at most the orbit's inclination. `planets.settle_obliquity` sets it at the end of the settle pass, replacing the share `tidal_spin.cassini_tilt_deg` drew early (moons keep that share: GEN.1536 part 2). Registry OB9 to OB12 added; SP22 narrowed to moons.
- Over 60 sampled systems (80 locked or resonant rocky planets) the median tilt goes from 0.158 to 0.017 degrees and the 95th percentile from 1.43 to 0.80 (Mercury is 0.034).
- Locked and resonant planets' tilts change, so the golden moves; no schema change, no new column (the figure is derived).

generation-output: changed
generation-affects: fill
