**Headline:** Two tests that were red on main are green: the job ledger no longer uses `VALUES()` in its upserts, and no item code shows in the text of a page or a log line.

### Fixed
- `test_describe_colour` no longer needs a gas giant among twelve default systems (they now almost never appear there); the stale decimal-time examples in docstrings read `4m 7s`.
- `ledger.begin_attempt`, `ledger.checkpoint` and `resumed._copy_checkpoint` write the new value (or take it from a derived table) instead of `VALUES(col)` in `ON DUPLICATE KEY UPDATE`, which MySQL deprecated.
- Item codes are cut from strings that reach a page, a sheet or a log line (the sheet's bookkeeping reasons, the redraw kinds, the asteroid vocabulary, a comet recipe label, a carbon-dioxide warning, the habitability table); the comments written into `schema.sql` by `body_fields.py` are registered as not shown.

generation-output: unchanged
