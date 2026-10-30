**Headline:** The meteorite groups (CI, CM, ordinary, enstatite and the rest) can now be told apart from a body's own composition instead of being fixed inputs (not yet used by the generator).

### Changed

- New library `physics/meteorite_groups.py`: the chondrite groups as points of a four-axis composition space (olivine oxidation, metal, refractory inclusions, matrix; Weisberg et al. 2006) with a nearest-group label and a `none` for bodies far from every group. Nothing in generation reads it yet; the wiring step (labelling belts and replacing the fixed within-family splits of the taxon kernel) follows with the epoch train.
