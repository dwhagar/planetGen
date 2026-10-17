**Headline:** A supernova remnant inside a superbubble is now sized by the cavity's thin gas, so the remnants of the cavities are far larger than those of the disk.

### Changed
- **A remnant's size reads the superbubble its sector lies in (GEN.298).** A remnant's radius (the Sedov-Taylor blast wave, held at the snowplow onset, GEN.151) took the density of its own class however thin the gas around it was, so a remnant inside a cavity was sized as if in the disk. Where the sector's profile places it in a superbubble, the remnant (rolled by the fill or built from the scatter) takes the cavity's interior density, 0.001 to 0.01 cm^-3, in its size law (`SystemConfig.AMBIENT_DENSITY_CM3`, `run_sector._bubble_density`), so it grows to a few hundred light-years where the disk's would have stopped at tens. The class contents and every draw are unchanged, and outside a bubble nothing changes. The remnant rate keeps its factor: the cavity is where the explosions went off.

generation-output: changed

TODO items resolved: GEN.298. Added: none.
