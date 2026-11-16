**Headline:** The test that built 12 real large systems to check the sector separation law is retired, about 1,100 s off a full run.

### Removed

- `tests/test_bughunt_sector_placement.py` (3 tests, 1,147 s recorded): it built a pool of real large systems to check that placed pairs keep the Hill-sphere separation, that a full sector raises, and that an explicit position skips the check. The law reads only a body's Hill sphere, and `tests/test_fuzz_sector_placement.py` tests the first two with random Hill sizes (`test_random_placement_stays_inside_and_separated`, `test_full_sector_raises_after_bounded_attempts_not_hang`).

### Added

- `test_explicit_positions_skip_the_separation_check_as_documented` in `tests/test_fuzz_sector_placement.py`: the third retired check, which the fuzz file did not make (two systems given overlapping explicit positions are both kept), with stand-in systems instead of real ones.

generation-output: unchanged
generation-affects: none
