**Headline:** Population is a ranking slot on the server: worlds with people rank by head count, the rest by habitability score, and the charted-regions list can carry a region's population once head counts are stored.

### Changed
- **Population as a metric slot (MAP.318).** New `db/population_slot.py`: `rank_key` and `rank` put a world with a population (largest first) before every world without one, and rank the worlds without one by habitability score (so the inhabited label slice is the habitable worlds until head counts exist), with ties by name and id; `population_of` treats a missing, zero, negative or NaN value as no population; `region_totals(conn)` is the head count of each region, empty while nothing stores head counts (the world population simulation, POP.30, fills it by summing per 9-block). `GET /api/galaxy/charted-regions` gives a region a `population` only when it has a count, so a list column for it stays hidden until then. Nothing changes on the map today. `docs/reference/api.md` says it.

TODO items resolved: MAP.318. Added: none.
