**Headline:** Small test clean-up from the retirement scan: the serialization fuzz runs five seeds instead of twenty, duplicate registry and sweep checks share helpers, and tests that only guard a removed thing are gone.

### Changed
- **TEST.734: `test_bughunt_serialization_fuzz.py` runs 5 seeds, was 20 (about 400 s off a full run), and its two identical phenomenon and compact-remnant round-trip tests are one parametrised test.**
- **TEST.735: the registry checks of the atmosphere and young-group draws share `tests/standin_registry_support.py`, and the three sweep files (atmosphere, bulk, class) share `sweep.assert_rules_obeyed`, `assert_every_rule_can_fail` and `assert_rules_are_complete`.** The tests keep their names and what they check.
- **TEST.736 (the part that touches no code): the wiki-routes-gone and no-wiki-url API tests (they would block the 9.1 spakr link plan), the retired-K-window notes test, the equipment-tags-no-levels-link test, the class-table-has-no-age-ranges test, the comparison to the retired activity table (the Skumanich band check stays) and the one-time black-hole print are retired; three tests lose only their "is gone" line or name.** The two back-compat tests that guard a shim (`test_belt_dose`, `test_benchmark`) are left for TODO, who checks them against the layered-fix scan. Tests only; no generated output moves.

generation-affects: none
