**Headline:** The Sectors list now tells the cells rogue planets landed in from generated sectors: they are named from their cell, marked uncharted, placed in their region and given no density.

### Changed
- A sector row a drifting rogue planet made in a cell nothing was generated in (`store.uncharted_sector_id`) is an uncharted placeholder (no cell in the chart ledger). `GET /api/sectors` and `GET /api/sectors/<id>` mark it `uncharted`, name it from its cell under the naming key (`naming_key.uncharted_sector_name`: two words, never the raw hex designation, which stays as `designation`), and give it its `region`, named whether or not any of that region is charted. `GET /api/galaxy/uncharted` carries the same `name` and `region`. `GET /api/sectors` adds `uncharted_total`.
- The Sectors page lists a placeholder as "Uncharted" in the Systems column with no density, and the header counts the generated sectors with the uncharted number beside it. A placeholder's page is the uncharted-sector page, with its name, an Uncharted badge and its region.
- New tests: `test_uncharted_sector_rows.py` (database and API) and `test_web_uncharted_sectors_browser.py` (browser, with the accessibility and phone-fit checks).

generation-output: unchanged
