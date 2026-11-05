**Headline:** Every list and message of the API names a system, sector or planet by its Object ID instead of a row number.

### Fixed
- **The object lists (`GET /api/object-lists/<kind>`) gave row numbers for a sector, a system and the source row.** `sector_id`, `system_id` and `source_id` of each item are now printed IDs (the list's "Sector" column and its link use the sector ID), and the `sector_id` filter takes a sector ID. The "no species on planet" answer names the planet by its ID. An audit test (`tests/test_api_link_ids.py`) now walks the answers that feed system and body links and fails on any row number.
