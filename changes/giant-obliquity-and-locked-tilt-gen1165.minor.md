**Headline:** A giant planet can now lie on its side after a giant impact, and a locked planet or moon keeps a small tilt instead of exactly none.

### Added
- **The tilt a giant impact gives an ice giant (GEN.1165, from the GEN.1053 module).** A giant of 5 to 30 Earth masses has a recorded giant merger with the stand-in chance of 0.55 (GT3); the impact's mass, impact parameter and direction are drawn and the planet's axis turns by what the impact's angular momentum adds to its spin, on top of the Rayleigh tilt it had before (`giant_tilt.tilt_after_merger_deg`, row GT4). Over 200 sample systems 4 of 25 such giants stand past 60 degrees (the sweep had none in 312 giants, with the real range out to 98 degrees); a heavier giant is not drawn a merger, since an impact of a quarter of its mass cannot turn its spin, and keeps its Rayleigh tilt (largest 56 degrees). Not changed here: the spin period after the impact and the plane of the moons and rings, which stay GEN.1058.
- **The first Cassini state of a locked or resonant body (GEN.1165, row SP22).** The axis of a body the tide has locked, or caught in a resonance, lies at a share of its orbit's inclination from the orbit normal (log-uniform from 0.003 to 1, `tidal_spin.cassini_tilt_deg`); it was exactly 0. The median is a fraction of a degree and the p95 under 5, as the sweep expected (Mercury, 0.034 degrees against an inclination of about 8).

### Changed
- Generator epoch bump (numbered at merge): the tilt of every giant planet of 5 to 30 Earth masses and of every locked or resonant planet and moon changes, and so does every later draw of the system.
generation-output: changed
