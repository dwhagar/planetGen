**Headline:** The giant fixed-point test accepts a cold giant that starts on its fixed point and settles in the one step that shows it.

### Changed

- GEN.1797 (test only): `test_a_giants_temperature_is_a_fixed_point_and_does_not_depend_on_where_the_solve_began` asserted that each of the first six sampled giants took more than one step of the temperature solve. The seeds now give a cold giant (start 111.8 K, fixed point 111.9 K) whose air adds no greenhouse, so the start is already the answer and it settles in one step; the solve itself is unchanged (the other giants take 3 to 5 steps). Each giant now either iterates or ends where it began within `TOLERANCE_K`, and at least one of the six must iterate, so the one-step shortcut of GEN.1716 is still caught.
