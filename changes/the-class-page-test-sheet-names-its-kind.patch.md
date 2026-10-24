**Headline:** The class page tests' fake property sheet names its kind, so the six phenomenon page tests pass again.

### Fixed
- **Six phenomenon-page tests in `test_class_pages.py` were red on main (TEST.365).** The phenomenon page now builds its "At a glance" panel from the property sheet's `kind` (UX.188), and the test's fake `get_object_sheet` answered `{"groups": []}` with no `kind`, so the page raised a `KeyError`. The fake now answers with the kind of the reference it was asked for, as the real sheet does.
- **TEST.366 (`test_memwatch`) passes on main already.** All nine tests pass; the budgets that test passes leave only the step's own peak to judge (PERF.187 #3125 changed the test), so no code change was needed.
TODO items resolved: TEST.365, TEST.366. Added: none.
