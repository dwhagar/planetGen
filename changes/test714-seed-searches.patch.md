**Headline:** two more tests that pinned a seed for a rare case now search seeds for it, so a change to what the seeds draw stops breaking them.

### Changed

- TEST.714: `test_circumbinary_light` pinned run seed 100039 for a close pair with planets and `test_refinish_rocky` pinned seed 38 for a moon the planet's history had replaced, and both would break whenever a train moved the draws. The first finds its pair with `first_seed_with` (the bounded search of `tests/fuzz_support.py`) over 200 seeds and says so when none has one; the second sweeps 12 seeds and checks every moon of every one, which holds for the replaced moons too. The two `pytest.skip` lines of the same file ("no moon in the draw", "no giant in the draws") became searches that fail with the number of seeds tried, so a case that vanished is a finding, not a skipped test.

generation-output: unchanged
