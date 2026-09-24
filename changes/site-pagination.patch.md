### Changed
- **Every list on the site now pages 50 rows at a time with the same
  pager.** Browse's sectors and standalone systems, Phenomena, a sector's
  systems and nearby phenomena, a Galaxy Map Quadrant's sector list, each
  Search result panel and the admin API key list all share one control
  (`src/html/lib/pagination.py`): a "Showing X-Y of Z" summary, First/Prev,
  numbered pages, Next/Last. Phenomena no longer stops at 500 rows and
  Search no longer stops at 300 matches per panel; every match is
  reachable a page at a time.
- `GET /api/search` pages each result panel: `limit` (default 300, as
  before) plus `sectors_offset`/`systems_offset`/`stars_offset`/
  `planets_offset`/`moons_offset`/`belts_offset`, and each panel now
  reports `total`, `limit` and `offset` alongside `truncated`.
