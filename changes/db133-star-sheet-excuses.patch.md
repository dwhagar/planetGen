**Headline:** The star property sheet says why it does not show the two columns DB.133 added to the stars table.

### Fixed

- `describe/sheet.py`: `NOT_SHOWN["star"]` gains `initial_mass_sol` and `phase_end_age_gy` with their reasons, so `tests/test_property_sheet.py` finds every column of `stars` shown or excused (a fix for Train 39).

generation-output: unchanged
generation-affects: none
