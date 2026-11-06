**Headline:** Four tests that went red when rogue planets, the chart ledger and the designation history arrived are fixed; no code of the generator changes.

### Changed

- The fingerprint's table check now names every table of the schema: `chart_ledger` and `designations` are left out on purpose (which run charted a cell and when; what edits retired), and the four rogue planet child tables are left out with a note that reading them is a lock item, because it changes every digest of a galaxy with a rogue and so the golden galaxy.
- The fresh-after-CLI tests count the sectors the CLI charted (the rows in the chart ledger), not the uncharted rows rogue planets create (GEN.1111, part 7), and look for a CLI-added cloud in the sector's own node, since a root node keeps only the biggest clouds and the seed's field clouds (MAP.178) fill it first.
- The database check's fixture sector has planets again (`_sector(..., planets=True)`): since GEN.1111 the shared helper makes none, so the tests that damage a planet changed nothing.
