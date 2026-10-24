**Headline:** The stellar-particle shielding and radiation-weight kernels for the derived surface radiation are in place (GEN.444, slice 1).

### Changed

- New `planetgen.physics.surface_dose`: the Band spectrum of a stellar particle event with the share that gets past an atmospheric and a magnetic cutoff, the shielding of a mixed column by nuclear interaction length, and the ICRP 103 radiation weights. Nothing calls it yet and no stored value changes; the wiring, the epoch and the removal of the proxies are the second slice, after GEN.443, GEN.441 and GEN.459.
