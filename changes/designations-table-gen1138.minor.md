**Headline:** An object that is renamed keeps the name it had, and the table that will keep every old ID is in place: a renumbered object gives its old ID up for good.

### Added
- **The `designations` table (GEN.1138, migration 0197).** One row for each ID or name an object used to have: the retired text, the object that held it (the table and its row number), the epoch it began and stopped in, why (collision bump, rename, regeneration, parent change or admin renumber) and, for a deleted object, when it was deleted. A retired ID is unique in the table, so it can never be given to an object again; a retired name is not reserved, because the collision bump hands a stored object's name to another by design. Nothing is backfilled: a name or ID that changed before the table has no row.
- **The writers (`planetgen/db/designation_rows.py`).** `set_name` is now the one place a `name` of a sector, system, star, planet, moon, belt, comet or phenomenon changes after it is saved, and it keeps the old name with its reason; `renumber` gives an object a new ID, refuses one that is any object's or ever was (even the object's own old one), keeps the old one, and moves the object's `object_index` row and its constellation and lock pointers with it; `history`, `retired_uid` and `retired_names` read the table back as the rules of `galaxy/designations.py` expect.

### Changed
- `store.rename_star_system`, `rename_star`, `rename_body` and `rename_phenomenon`, the names that follow a system (its stars, planets, moons and comets), the Greek and Roman decoration of a sector or system that collides, the collision bump of an object's numbered name and the healing of a system name now keep the name they replace. A person's change is recorded as `rename`, every automatic one as `collision bump`. A system built again keeps nothing for the placeholder names its new bodies carry for a moment. The name an object ends up with is the same as before.

No generator output changes.

generation-output: unchanged

TODO items resolved: GEN.1138.
