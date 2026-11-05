**Headline:** A new docs test fails when an item's Version line in the TODO page disagrees with the version file its row sits in.

### Added
- `tests/test_todo_docs.py` checks that the last `Version 8.2.` style line in an item's text in `docs/TODO.md` names the version whose file (`docs/plan/releases/version-*.md`) holds the item's row, and fails with the item, the version it states and the file found. It reads only an item's own statement (a capital `Version`, a number, then a full stop, comma or similar), not sentences that mention another version.
