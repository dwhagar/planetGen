**Headline:** The near-search test builds its scattered-star and black-hole rows by column name, so a new column no longer breaks it.

### Fixed
- **`test_near.py::_scatter_rows` built its `bright_stars` row as a 19-value tuple.** GEN.846 (schema v176, commit 04ed93b4d) added `rotation_period_hours` as the twentieth column of `store.BRIGHT_STAR_COLUMNS`, so the insert failed with "Column count doesn't match value count" and `test_the_search_counts_what_is_not_generated_yet` was red on main. The `phenomenon_scatter` tuple was still the right length (14) by luck.
- **Both rows are now built by name** from `store.BRIGHT_STAR_COLUMNS` and `store.PHENOMENON_SCATTER_COLUMNS` through a small `_row_of` helper: a column the test does not set is `None`, and a name that is not a column fails the test with its name. A column added to either table later cannot shift the values after it again.

Tests only. No generator output changes.

generation-output: unchanged

TODO items resolved: none (TEST.534 closes when this merges; the TODO thread owns the file).
