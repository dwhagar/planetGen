**Headline:** Every generator epoch now says which stages of generation it changed, the first step to bringing an existing galaxy up to date without regenerating all of it.

### Changed
- **Each epoch lists the stages it changed (GEN.764).** `src/tests/golden/epochs.json` gives every epoch after the first an `affects` list from `plan`, `scatter`, `fill`, `population` (the deepest stages its change altered), and lists the stages at the top. The 96 existing epochs were classified by reading their reasons, conservatively: when in doubt a stage is listed, and an epoch that changed the galaxy's shape or density counts as `plan`.
- **`python -m tests.golden_update --reason "..." --affects fill`.** A new epoch needs `--affects` (checked before the ten-minute rewrite of the golden files starts); a second run in the same change can change the pending epoch's list. `scripts/bump_version.py --check-pr` refuses an epoch without a valid list.
- **`planetgen.util.generator_epoch`** gains `stages()`, `affects(epoch)`, `stages_between(recorded, current)` (what a result made under an older epoch has to redo) and `needs_full_regeneration(recorded, current)`, and one `format_epochs` that every writer of the file uses so the short lists stay on one line.
- No generated output changes.

TODO items resolved: none (GEN.764 step 1 of 4). Added: none.
