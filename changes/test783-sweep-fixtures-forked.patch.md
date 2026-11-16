**Headline:** the two sweep tests that built their galaxy one system at a time now build it with four workers at once, so they finish inside the hang guard.

### Fixed

- TEST.783: `test_star_composition_sweep` filled its five region sectors one after another, and `test_radius_valley_sweep` made its 60 G2V systems one after another on one run seed; both ran past the hang guard (900 s, 1,500 s in a gate) on a 4-core box and errored instead of testing. The composition sweep now fills each region in a forked worker (each fill was already bound to its own seed, so the stars and their keys are the same as before). The valley sweep makes 4 blocks of 15 systems at the same time, each on a run seed of its own, since the systems of one run seed draw from one another's stream; the blocks are numbered 4 to 7 and the new comment on `FIRST_BLOCK` says why. What the sweeps check and the bounds they check are not changed.
