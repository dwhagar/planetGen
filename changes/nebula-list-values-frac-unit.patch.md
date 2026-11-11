**Headline:** The nebula composition range test of the object lists passes the unit the columns are stored in.

### Changed

- `tests/test_object_list_values.py::test_the_nebula_list_filters_and_sorts_by_the_stored_composition` gave the nebula composition ranges (`mass_fraction_h`) the unit `ratio`, but those columns are stored as `frac` (`db/column_units.py`) and `search_values.convert` refuses `ratio` to `frac` ("The two measure different things"). It has been red on main since the test came in with GEN.1276; the test now passes `frac`. No code changes.

generation-affects: none
