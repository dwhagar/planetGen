**Headline:** The hand-set shape of a star's own comets is registered with a tier and a tweak level.

### Changed
- GEN.595: the period cuts, outer period ends, eccentricity ranges, inclination caps and family weights of `COMET_PERIOD_CLASSES`, the nucleus diameter range, the perihelion range, the chance a star holds comets and how many are now rows CMT1 to CMT8 of the stand-in registry (section 6.303), each with its tier, tweak level, knobs, falsifier and the data that would replace it, and tagged in `tuning.py`. Only the parabolic share was registered before (AUD20). The period cuts and the flat perihelion are derived from the dynamical class definitions and from a loss cone; the rest are named as hand set, with a derived draw (a power-law nucleus size, a supply from the outer belt) left as its own item. Nothing generated changes.

generation-output: none
