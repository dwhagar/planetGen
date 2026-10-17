**Headline:** Stars under 10 Myr are protostars and pre-main-sequence stars with a disk, instead of main-sequence stars that happen to be young.

### Added
- **The pre-main-sequence state (GEN.242).** A star under 8 solar masses and under 10 Myr is a Class 0 or I protostar (under 0.5 Myr) or a pre-main-sequence star, with luminosity, temperature and radius from a contraction track worked out from the energy it releases (a solar mass reaches the main sequence in about 40 Myr, 0.1 solar masses in about 0.9 Gyr). The star page names the stage ("Class I Protostar", "Pre-Main-Sequence Star") and says the star is still contracting.
- **The disk descriptor.** A young star holds a Class 0, I or II disk (age under 0.15, 0.5 and 3 Myr; none after) of 0.01 to 0.3 of its mass, from the dust sublimation radius out to 10 to 300 AU, with no bodies in it. It is stored in `stars.disk_kind`, `disk_mass_solar`, `disk_inner_au` and `disk_outer_au` and shown on the star page.
- **Companions of protostars.** A Class 0 protostar has a companion in 57% of systems and a Class I protostar in 23% (the Perseus survey), in place of the field's share by spectral class.

### Changed
- The companions of a young system (`multistar`) take their stage and radius from the same track; the old `31 Myr x m^-2.5` law is gone (it gave 10 Gyr at 0.1 solar masses).
- A star of 8 solar masses or more is never a protostar: it is on the main sequence by the time it shows.
- A young star counts as saturated in X-rays, and has no planets inside the engulfment radius of a giant.

generation-output: changed

TODO items resolved: GEN.242. Added: none.
