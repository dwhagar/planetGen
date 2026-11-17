**Headline:** An edit to a sector now updates its parity group from a before-image instead of reading every other sector of the group, and three writers that marked a sector without telling the parity file now do.

### Fixed
- **Three writers left the parity file behind (DB.109, part 1).** `store.touch_sector` (a deleted system, a rogue top-up cleared), `store.set_sector_uid` and the sector `PATCH` route set the sector's `parity_dirty` flag with a plain `UPDATE` and never told the connection to bring the parity file up to date after the commit, so the file stayed behind until `planetgen refresh-parity`. They now mark through `storage_checksum.mark_sector`, like every other writer. `reserve_sector_name`'s rename of a sector no longer sets the flag a second time. `tests/test_parity.py` runs each way a sector's rows change (rename, facility, edit, replace, insert, rogue top-up, empty and delete a sector, the uid writers) with parity on and checks the file holds every sector, none is marked and every record matches its rows.

### Changed
- **An edit is a parity delta (DB.109, part 2).** The admin and property-sheet edit paths (`edits.save_system_edits`, `save_star`, `replace_phenomenon_content`, `delete_phenomenon`, `empty_sector`; `store.rename_star_system`, `rename_star`, `rename_body`, `rename_phenomenon`, `replace_star_system_content`, `add_facility`, `delete_facility`; the system and sector routes) take the sector's export before they write (`Connection.parity_before`, `parity.before_image`), and the commit hands it to the parity file, which updates the group's parity from the change. Before, every edit or deletion of a sector in a group read the other members to make the group's parity again. A sector that is still marked, damaged or not in the file gives no image and is made again from its rows as before. Nothing happens while `parity.enabled` is off, which is the default.

generation-output: unchanged
generation-affects: none
