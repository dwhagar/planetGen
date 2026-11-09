**Headline:** The repair test reads the inserted-row count as a whole number, so "90 inserted" no longer fails the check for "0 inserted".

### Changed

- `tests/test_repair.py`: the check that the reason does not say "0 inserted" matched inside "90 inserted" now that the repaired sector reports 90 rows; it uses a word boundary (`\b0 inserted`).

generation-output: unchanged
