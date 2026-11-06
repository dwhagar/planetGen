**Headline:** A sector's gravity grid is now available from the API, as the pull, the well depth or the tidal strength of each of its 4,096 zones.

### Added
- **`GET /api/sectors/<id>/gravity` (API.24).** Answers one number per zone of the sector's 16 x 16 x 16 grid for the mode asked (pull in m/s^2, well depth in m^2/s^2, tidal strength in 1/s^2), with the zone size, the corner of the grid, the range, the saddle zones between wells and the stamp of the grid. The default is the total field (the generated stars with the smooth galaxy added); `field=local` gives the stars' own field, which the smooth galaxy otherwise swamps in the pull. A grid that is not worked out yet, or that a moved star made stale, is queued as a job and the route answers 202 with its status URL.

### Changed
- **The part of the gravity grid job that needs no compiled code moved out of it.** The sector reach, the cache stamp and the zone centers are in `db/gravity_stamp.py`, and the hierarchy's settings in `galaxy/gravity_settings.py`, so the web app can tell whether a grid is cached without importing numba. Grids already cached keep their files (the stamp is the same).
