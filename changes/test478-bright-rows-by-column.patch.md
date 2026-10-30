**Headline:** The bright-star rows of the API tests are built by column name, so a new column cannot leave them short again.

### Changed

- TEST.478: `tests/test_api.py` builds its `bright_stars` rows with one helper that fills every `store.BRIGHT_STAR_COLUMNS` name from defaults and the values a test names, in the table's order. The rows had been one value short since GEN.846 added the rotation period (23 tests were failing on main with column count 1136 in the bright-star routes).
