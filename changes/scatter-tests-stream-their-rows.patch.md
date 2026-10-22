**Headline:** Five scatter tests that held a million rows or more in lists now stream them, so none holds more than a gigabyte.

### Fixed
- **Scatter tests of the toy galaxy held gigabytes (TEST.258).** Since GEN.404 the toy galaxy holds a real galaxy's stellar mass, so a shell of its scatter is a million rows and a pass over the far shells nearly two million. Measured alone, `test_the_far_pass_has_its_own_stream_and_is_reproducible` held 2.8 GB, the black hole stream test 1.4 GB, `test_the_same_seed_gives_the_same_rows_and_another_seed_others` 1.3 GB, `test_sectors_already_filled_get_no_row` 1.3 GB and the 1500 to 1540 case of `test_the_far_pass_agrees_with_the_shells_one_by_one` 1.1 GB. They now read the rows once and keep a count and a hash (`tests/row_digest.py`), and the largest of them holds 0.84 GB. The other 20 test files that scatter measured under 1 GB alone, the largest being `test_the_classes_every_component_has_keep_their_draws` at 0.98 GB.

TODO items resolved: TEST.258. Added: none.
