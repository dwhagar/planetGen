**Headline:** Each star system now stores what its rows add up to (star, planet and moon mass, luminosity, body counts), so the gravity map and the sector paths stop summing them again.

### Changed
- **`system_summary` holds the figures (GEN.1478, migration 0237).** For each star system: the stars' mass, the stars', planets' and moons' mass, the brightest star's and all the stars' luminosity, and the counts of stars, planets, moons, belts and comets (`db/system_summary.py`). The system save and the editor write them; the migration fills every system already stored, in runs of 20,000, with progress shown.
- **Two readers use them.** The gravity sources (`db/gravity_sources.py`) read `system_mass_kg` and the sector paths (`db/sector_paths.py`) read `stars_mass_kg`, instead of grouping the `stars`, `planets` and `moons` rows of every system in the neighbourhood each time. The object index and the sector totals still sum their own rows (the sector total is an exact decimal sum and must not go through a stored double).

generation-affects: none
