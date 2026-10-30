**Headline:** The rules of designation history are written and tested as pure code: a retired ID is never reused, an old ID resolves to the current object, and a system's ID change moves every body.

### Added
- **Designation rules (GEN.1138, pure part).** `planetgen/galaxy/designations.py`: the closed list of reasons a designation is retired, `check_new_uid` (an ID that is live or was ever an ID is refused), `cascade` (the new IDs of a system's bodies when its ID changes), `resolve_uid` (current, redirect, gone or unknown), `resolve_name` (the live holder first, then the latest past holder that still exists), `resolve_member` (a belt member's or captured body's old designation follows its parent) and `former_lines` (the "Formerly" list, newest first). Nothing calls it yet; the `designations` table, its migration and the writers follow in the same item.

No generator output changes.

generation-output: unchanged

TODO items resolved: none (GEN.1138 stays open for the table, the migration and the writers).
