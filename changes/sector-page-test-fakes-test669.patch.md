**Headline:** The tests of the sector page and the section pages pass again: their stand-ins for the web layer's data calls had fallen behind two calls the pages learned.

### Changed

- `src/tests/test_web_sector_nav.py` and `src/tests/test_web_pages.py` (TEST.669): 35 tests answered a 404 page ("No such database: 'planetgen_web_test'") or a TypeError on main. Neither was a route or an ID-format regression. MAP.381 made the sector page ask for the globular clusters (`apiclient.get_containers`) to draw the ones that reach into the sector, and the sector test's fake data never stubbed that call, so it went through the in-process transport to the real API, which found no such database. MAP.403 gave `get_sectors` the `coverage` and `show_empty` keywords, and the section pages' fake `get_sectors` did not accept them. The sector test's fake now answers `get_containers` with no clusters and the other fake's `get_sectors` takes both keywords. The pages and routes are unchanged.

generation-output: unchanged
