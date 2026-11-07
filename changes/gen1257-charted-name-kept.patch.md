**Headline:** When charting gives an uncharted star its own name, the name its Object ID spelled is kept on the Formerly line, and a search for it still finds the star.

### Added
- A new retired-name reason, `charted` (migration 0221 widens `chk_designations_reason`; "named when the star was charted" on the Formerly line). The name an uncharted star's Object ID spells (a local name, `naming_key.local_name`, GEN.1255) is computed and never stored, so when the sector fill, or `POST /api/uncharted-systems/<id>/generate`, builds the system round a scatter star and the registry gives it a word-salad name, `store.keep_charted_name` writes the spelled name as the system's first retired name. `api.lookup` follows it, so an old link or a past-name search lands on the system. The Object ID does not change, so nothing redirects by ID.
- Nothing is kept for a system the fill made without a scatter star, when no naming key is in force, or twice for the same system. Scattered phenomena keep no word-salad name to replace (their names come from their IDs, GEN.1256).

### Changed
- Stored designation rows gain one `charted` row for each charted scatter star; no other stored value changes.

generation-output: changed
TODO items resolved: GEN.1257 (migration 0221, epoch and golden ride the train; the new rows are in `designations` only)
