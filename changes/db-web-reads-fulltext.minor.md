### Changed
- **Search matches whole words (PERF.16).** Name search uses a FULLTEXT
  index (schema v46) instead of scanning every row for the text, so
  "ara" no longer finds "Kemaral", while "Mu" or "Tobar IV" still find
  "Ossiran Mu" and "Tobar IV". The Galaxy Map's locate box still finds a
  name as you type the start of its last word. Result counts stop at
  300 and show "300+".
- **Web pages run fewer queries (PERF.15).** System pages load moons in
  one query, sector pages load their stars in one query, system lists
  read star types once per page, and search's filter lists are cached
  until a sector or system changes.

### Added
- **A time limit on web database queries (PERF.17).**
  `mysql.statement_timeout_seconds` in `config.json` (default 10, 0
  turns it off) stops any one query on the web interface and API; the
  page says "Took too long" (HTTP 504) instead of hanging.

Run `update.sh` (or `update.ps1`) after updating: it migrates the
database to schema v46, adding full-text indexes to five tables. Writes
to those tables pause while each index builds, which can take a few
minutes on a large galaxy.
