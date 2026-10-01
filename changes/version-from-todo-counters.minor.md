### Changed
- **Version numbers are now MAJOR.REVISION.BUILD** (OPS.1). A `major`
  release note bumps MAJOR, any other note bumps REVISION, and BUILD is
  the sum of the TODO category counters in the new "Next free IDs" table
  of `docs/design/todo-number-map.md`, so it tracks how many TODO items
  have ever been filed. `bump_version.py --check` fails when `docs/TODO.md`
  uses an ID that table hasn't counted. See `changes/README.md`.
