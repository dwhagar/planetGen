**Headline:** The tier plan shows a progress bar for every tier: items finished out of items in the tier, with the percentage, rebuilt on every retire.

### Changed

- `scripts/build_todo_docs.py` reads a new `docs/plan/progress-ledger.md` (the phase of each finished item, read back from the git history of the phase files) and the tier plan page shows finished, open and percentage per tier at the top and in each phase heading. The build fails when a retired item has no ledger row.

TODO items resolved: none / Added: none
