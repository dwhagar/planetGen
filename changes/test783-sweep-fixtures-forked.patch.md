**Headline:** the two sweep tests that built their galaxy one system at a time now build it with four workers at once, so they finish inside the hang guard.

### Fixed

- TEST.783: `test_star_composition_sweep` filled its five region sectors one after another, and `test_radius_valley_sweep` made its 60 G2V systems one after another on one run seed; both ran past the hang guard (900 s, 1,500 s in a gate) on a 4-core box and errored instead of testing. The composition sweep now fills each region in a forked worker (each fill was already bound to its own seed, so the stars and their keys are the same as before). The valley sweep makes 8 blocks of 15 systems (120 systems, twice the sample) four at a time, each block on a run seed of its own (blocks 0 to 7, fixed in advance), since the systems of one run seed draw from one another's stream; the comment on `BLOCKS` gives the measured share and its spread. What the sweeps check and the bounds they check are not changed.
