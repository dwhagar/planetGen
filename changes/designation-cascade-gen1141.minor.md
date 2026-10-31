**Headline:** When a system's ID changes its bodies follow, and when an object is deleted its IDs are kept as gone: an old link says so (410) and the ID is never given to anything else.

### Added
- **Cascade and tombstones (GEN.1141, part e).** `designation_history.renumber` of a system renumbers its stars, planets, moons, belts and comets in the same transaction, each keeping its body number under the new system ID and getting its own history row; a body may only take an ID of its own system; an object with an edit lock keeps its ID until the lock goes. A belt member's or captured body's designation (`<belt ID>.<rank>`) is derived from its parent, so it follows without a row, and an old one still names the member (`resolve_member`). Every delete (a system with its bodies, a sector with its contents, a phenomenon, a sector, a body removed by an edit) now writes a tombstone first, so the old link answers 410. When a system's content is written again (`replace_star_system_content`) a body's history stays with the body and a body that did not come back is a tombstone; the stand-in names of that write are not kept. The `deleted` reason joins the closed list.

No generator output changes.

generation-output: unchanged

TODO items resolved: none (TODO closes GEN.1141 when the branch merges)
