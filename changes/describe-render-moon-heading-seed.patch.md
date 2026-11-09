**Headline:** The heading test of the describe renderers no longer needs the first seeded system to have a moon.

### Changed

- `tests/test_describe_render.py::test_headings_follow_the_tree_in_each_format` took the first of three seeded systems and asserted it had a moon section. `MOONS` forces a chance at moons, not a moon, and the Train 21 merge changed what the seed draws (the system of run seed 100 has five planets and no moon now; one moon before). Measured over 80 systems with `MOONS` on, 46 had a moon before Train 21 and 48 after (340 and 354 moons), so the moon rate did not move. The test takes the first of the three seeded documents that has a moon, and fails if none has.

generation-affects: none
