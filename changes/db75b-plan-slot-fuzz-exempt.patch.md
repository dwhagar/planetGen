**Headline:** The serialization fuzz test knows a planet's plan slot is stored through its history row, not as a field of the planet.

### Fixed
- `test_starsystem_children_have_no_unlisted_attributes` listed `plan_slot` (DB.75b) as an unlisted attribute. It is stored as the body's `born_index` and given back on load, so it joins the test's runtime-only list for planets, with the reason.

generation-output: unchanged
