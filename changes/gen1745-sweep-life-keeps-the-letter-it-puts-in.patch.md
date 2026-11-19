**Headline:** Two Train 47 gate reds are fixed: the life sweep audit keeps the class letter it puts in (so its "same conditions as another class" check still compares something now that the class is read once at the end of the life pass), and the stand-ins of the water-limited weathering (WW1, WW2) have their registry section.

### Fixed

- `test_a_generator_that_skips_life_for_one_letter_is_reported` (a Train 47 gate red from GEN.1745 slice 3, PR 7211): the audit (`tests/audit/sweep_life.py`, `_run_life`) puts another class letter on a copy of a world and runs the life step to see whether the outcome follows the letter. Since the class is read once at the end of `life.apply_life_data`, the letter it put in was overwritten by the world's own class before the check, so every letter gave the same outcome and the audit compared nothing. While the audit runs the life step on a copy, `planets.apply_class` is now held off, so the copy keeps the letter it was given (or the world's own letter), as it did before the class moved to the end of the life pass. No generator code changes.
- `test_every_registry_table_row_of_the_plan_is_in_the_registry` and `test_every_tag_of_the_code_names_its_section_and_an_id` (PR 7218, GEN.1403): the plan's section 6.309 (rows WW1 and WW2) was not in the registry's section list (`scripts/build_stand_in_values.py`), and the two tags in `surface_pressure.py` named no section. The section is added to the list, the tags read `[Stand-in 6.309 WW1, ...]` and `[Stand-in 6.309 WW2, ...]`, and the stand-in list and reference pages are regenerated. No value changes.

generation-output: unchanged
generation-affects: none
