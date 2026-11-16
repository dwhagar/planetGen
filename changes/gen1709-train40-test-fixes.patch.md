**Headline:** The inner-system spread (GEN.1709) is cited by its registry ID, and four tests that its new draws moved (citations, companion gas, Urey range, Neptune desert) are fixed.

### Changed
- GEN.1709: the occurrence survey is written `[R2020: Table 2]` in `tuning.py`, the stand-in row PA19 (and `stand-in-values.md`, regenerated), the architecture test and the release note, where it was named by author and year.
- Tests only, with the reasons in each: `test_close_companion_gas` builds 1,000 architectures per factor (200 gave 36 and 23 giants, a cut of 0.36 within the counting noise of 23; 1,000 give 178 and 84, a cut of 0.53); `test_a_generated_body_has_a_urey_ratio_above_one_only_while_its_mantle_is_still_warming` leaves out a body the tide heats more than its own heat does (the ratio is then a difference of two large heats; the same value with the new conduction off); `test_few_neptunes_sit_inside_three_days_and_a_few_survive_there` compares with 3 to 30 days (3 to 10 days held 93 Neptunes of the 1,500 chains when the first planet sat at 5 days and holds 25 with it spread; inside 3 days went from 8 to 13).
- No generation output changes.
