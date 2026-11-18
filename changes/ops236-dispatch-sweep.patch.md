**Headline:** the dispatch rule is proved over ten seeds and backlogs of 30, 300 and 3,000 units, with kills, small claims and a decider that fails.

### Added
- OPS.236: `tests/test_fuzz_dispatch.py` runs the simulator (`tests/enginesim.py`) over seeds 0 to 9 at 30 and 300 units and two seeds at 3,000. Each run has a visitor's jobs through the stage, a claim of 7,000 MB, a restarted unit that runs alone and a cut that needs the lane to itself. It checks in every turn that the claims admitted stay within the host less its reserve (a claim larger than that runs alone), that a token is held by one job, that the lone job runs alone and that a kernel kill only follows a claim that was too small or memory taken from outside; at the end that every job finished once, that the visitors' jobs wait under 1.5 s at the 99th percentile (5 s with a decider that fails one turn in ten), and that the big claim, the lone unit and the cut all start within two minutes of the end of the work queued before them.
- The simulator gained a claim that is too big for its neighbours (`big`), a decider that fails in a share of the turns (`crash`), the share of units whose claim is too small (`small`), memory held from outside the engine for a while (`pressure`), and a retry that runs alone on a claim grown by half after a kernel kill.

generation-output: unchanged
generation-affects: none
