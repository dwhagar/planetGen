**Headline:** A rock's kind (basalt, granite, anorthosite and so on) can now be worked out from the oxides of the rock and how much of it melted, instead of copied from the Earth's crust (not yet used by the generator).

### Changed

- New library `physics/rock_types.py`: the CIPW norm of a set of oxides, the melt of a rock by stages with the enrichment of K, Th and U, and the ground type of the result by the standard silica bands. Nothing in generation reads it yet; the wiring step (the surface draw taking its crust from the stored mantle) follows with the epoch train.
