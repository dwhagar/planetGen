**Headline:** Epoch and migration numbers are renumbered by one script when a branch has to follow main, instead of by hand.

### Changed
- **`scripts/allocate_numbers.py` renumbers a branch's migrations and epochs after main took the numbers it had.** Run on a tree that holds the branch's changes (merged or rebased onto main, even with `epochs.json` in conflict), it renames the branch's migrations to follow main's newest and rewrites their revision ids, appends the branch's epochs to main's `epochs.json` with the next numbers, makes the `schema_v<head>.sql.gz` fixture the first migration is tested from out of main's `schema.sql` and moves its README line, and changes `v<old>`, `<old>` and `epoch <old>` in the lines the branch added. `--check` reports whether the numbers already follow main's and `--dry-run` shows the changes; it refuses when two own migrations share a number or main already holds the fixture. A test builds a main and a branch in a throwaway git repository and covers both the renumbering and the conflicted `epochs.json` (OPS.71).

TODO items resolved: OPS.71.
