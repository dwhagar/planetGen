**Headline:** The test suite loses its slowest duplicate: the real-system sector placement bug hunt, about 1,100 seconds a full run, whose three laws the fuzz file already tests.

### Changed
- **TEST.730: `test_bughunt_sector_placement.py` is retired.** It built up to 12 real large systems to check three laws (placed pairs keep the Hill-sphere separation, a full sector raises, an explicit position skips the check), and the separation law reads only a body's Hill sphere size. `test_fuzz_sector_placement.py` tests the same three with random Hill sizes in 81 seconds: `test_random_placement_stays_inside_and_separated`, `test_full_sector_raises_after_bounded_attempts_not_hang` and `test_explicit_positions_on_the_cube_boundary_are_inside`. The recorded 1,147 seconds leave `scripts/data/gate_durations.json`. Tests only; no generated output moves.

generation-affects: none
