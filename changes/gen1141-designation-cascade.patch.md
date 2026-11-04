**Headline:** When a system's ID changes, every body in it takes the new ID in the same step, and each keeps its old ID on record, so an old link still finds the right body.

### Added
- `designation_rows.renumber_system` gives a system, or a free planet, a new ID and its stars, planets, moons, belts, comets and facilities the new ID with their own body numbers, in one transaction. Each body writes its own row in `designations` with the reason `parent change`, the system's row keeps the caller's reason, and a new ID that is live or reserved for any of them stops the whole change. Nothing calls it yet; the 8.6 edit flow will.
- `designation_rows.member_now` turns the old text of a belt member or a captured body (`<parent ID>.<rank hex>`) into its current text after the parent's ID changed. Members are stored nowhere, so they follow the belt's or the planet's ID by themselves and need no row of their own.
- Tests in `tests/test_designation_rows.py`: the bodies follow, each keeps its old ID reserved, a member's old text redirects to the new one, and a refused cascade changes nothing.

generation-output: unchanged
TODO items resolved: GEN.1141 cascade half (the tombstone on delete waits on GEN.1123's audit of the delete paths)
