**Headline:** two tests that the job engine turned red now say what a held run really does, and the finish-run comparison no longer mistakes a row number for content.

### Fixed
- OPS.277 (bug, test only): `test_maintenance`'s map rebuild test named a database ("x") that does not exist. Since every unit's attempt is written to the `engine_items` table of the unit's own database (OPS.217), a run on such a name is held with a connection error, which is right: a run that cannot reach its galaxy cannot fence its attempts. The test now uses a real, empty database and reaches the step it is about (reopening the site) again.
- OPS.277 (bug, test only): `test_finish_run_generation` compared the knots of the sector paths by the path's row number, which is not content: a run settled in two stretches numbers its paths differently from one settled whole, so a galaxy finished after a hold looked different from one made without a stop. The knots are now compared by the object each path is of. The three tests were not hung: each takes about 15 minutes alone, and the shared gate's 1,500 second cap ended them under load.

generation-output: unchanged
generation-affects: none
