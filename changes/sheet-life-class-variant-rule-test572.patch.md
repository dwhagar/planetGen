**Headline:** The property sheet knows that only some life classes have a variant letter, so a body whose class has none is no longer flagged as an empty column.

### Changed

- `life_class_variant` is empty on a body whose life class comes in one form (classes 2, 6, 7, 8 and 9 have no pigment or habitat letter); only classes 1, 3, 4 and 5 have variants. The body sheet rules now say so, derived from the life class table, so `test_a_sample_has_no_empty_column_that_applies_and_is_not_filed` no longer fails on a class 2 or class 6 body (TEST.572). The other half of TEST.572, `stars.nuclear_arrival_age_gy` empty in every row, was already fixed by TEST.556 (#4756): both checks passed on main in two runs.

generation-output: unchanged
