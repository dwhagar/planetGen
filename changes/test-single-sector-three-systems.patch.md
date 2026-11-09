**Headline:** The deferred-sectors test of a single sector that exists generates 3 systems, not the 353 the wide galaxy's density gives (a test-only change).

### Changed
- **`test_a_single_sector_that_exists_is_not_cut_again_and_one_beyond_the_reach_is_refused` sets `args.num_systems = 3` before it generates the sector.** The test only needs the sector to exist so that the cut leaves it out, but at the density the wide galaxy gives address (1, 0, 0) it built 353 full systems, about 15 minutes on a gate machine, over the gate's 900 s timeout. It was not a hang and no code changed: it passes alone on main in 888 s, and now in 15 s. The command line refuses `--num-systems` with an address, so the test sets it on the parsed arguments after the cut.

generation-output: unchanged
generation-affects: none
