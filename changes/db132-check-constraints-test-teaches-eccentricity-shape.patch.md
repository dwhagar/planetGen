**Headline:** The CHECK-constraint test knows the shape of the two eccentricity columns DB.132 added to body_history, so it breaks and sees refused a bad value for each.

### Fixed

- `tests/test_db_check_constraints.py`: `_breaking_value` reads `(a is null) or ((a >= 0) and (a < 1))` (`held_eccentricity`, `forced_eccentricity`) and breaks it with 1.0, which the CHECK must refuse (a fix for Train 39).

generation-output: unchanged
generation-affects: none
