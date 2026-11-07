**Headline:** The Sectors list and the search show charted sectors by default, with a choice of all sectors or only uncharted ones, and a separate option for the uncharted sectors nothing has been detected in.

### Changed
- `GET /api/sectors` takes `coverage=charted|all|uncharted` (default `charted`) and `empty=1` (the uncharted sectors nothing is known to be in; default left out). `total` counts the sectors listed; the answer adds `charted_total`, `empty_total`, `coverage` and `show_empty` beside `uncharted_total`. The default listing no longer includes uncharted placeholders (the API version is unchanged: no route moved).
- The Sectors page has a "Show" menu (Charted sectors, All sectors, Uncharted sectors only) and an "Include empty uncharted sectors" box, kept through sorting and paging; the header counts charted sectors, the uncharted ones with something detected in them and the empty ones, each linking to its list.
- The search takes `sector_coverage` and `sector_empty` for its Sectors panel (`GET /api/search`, and a menu and box under "Search by object and size"), with filter chips; by default only charted sectors are searched.
- `query.list_sectors`/`count_sectors` take `coverage` and `show_empty` (keeping every row when not given); new `count_empty_uncharted_sectors`.
- New tests: `test_sector_coverage_filter.py`; `test_uncharted_sector_rows.py` asks for the placeholders it lists.

generation-output: unchanged
