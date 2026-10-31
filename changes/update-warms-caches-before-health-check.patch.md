**Headline:** `update.sh` says what it is doing and what it waits for, shows progress bars on the long steps, and pauses 5 seconds before the reload and 10 seconds after it.

### Changed
- Every wait names what it waits for and for how long at most: the warm-up's opening phase (with the last line of its log every 30 seconds), the pause before the reload (5 s, `UPDATE_SETTLE_SECONDS`), the pause after it (10 s, `HEALTH_SETTLE_SECONDS`) and each health probe retry. The math check, the app import, the Apache configuration test and the cluster-mass step say what they are about to do.
- The cluster-mass step (`planetgen.cli.warm_shape`) shows a progress bar on a terminal and a line a minute otherwise.
