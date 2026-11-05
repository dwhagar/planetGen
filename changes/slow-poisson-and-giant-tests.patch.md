**Headline:** The halo and block-path Poisson tests, the giant-tilt test and the wide-pair velocity test now run in two to five minutes each instead of five to nine; no generated output changes.

### Changed
- **Sampling tests drawn only as far as their checks need (TEST.417).** `test_gen_halo.py` takes every second cell of the halo and edge set at 150 times the rate, still more than 100 stars to compare against the Poisson mean; `test_backfill_block_path.py` draws only as many seeds (ten at most) as it takes the expected stars to pass 100, since each seed costs about 30 seconds; the giant-tilt test in `test_tidal_spin.py` stops once it holds 8 ice giants and 6 gas giants (it asked for 4 of each from 60 systems that cost 7 seconds apiece); the wide-pair test in `test_system_velocity.py` places each of 10 wide pairs at three positions (30 checks, not 30 pairs built). Every property and every threshold is kept; none is skipped or loosened.

generation-output: unchanged

TODO items resolved: none (TEST.417 continues).
