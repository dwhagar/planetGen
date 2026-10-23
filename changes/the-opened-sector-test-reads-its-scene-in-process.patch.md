**Headline:** A star-selection browser test no longer times out on a loaded machine, and the test databases' reset notices a foreign key given another rule.

### Fixed
- **A browser test of selecting a star in an opened sector timed out at 30 s (TEST.272).** `test_a_star_in_an_opened_sector_is_selectable` fetched the sector's scene JSON through the served site, a request that waited behind the page's own work while two browsers drew with the software renderer. The sector browser tests now ask the app for the scene in the test process, like the sector designation already is (TEST.195).
- **A test that gave a foreign key another delete rule left the next test a damaged schema (TEST.270).** The key kept its name and its table's creation time, so the worker's database reset took the schema for unchanged. The structure signature now includes every foreign key's update and delete rule, and the database is thrown away when one differs.

TODO items resolved: TEST.270, TEST.272. Added: none.
