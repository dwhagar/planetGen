**Headline:** A system now has a sky route that gives the stars it sees, the nebulae and supernova remnants in its sky, and how much of the sky is built, in one cached file.

### Added

- `GET /system/<id>/sky?limit=` (the page's route) and `GET /api/systems/<id>/sky?limit=` (VIEW.31) build the sky a system sees from the stored stars (`api/system_sky.py`: `galaxy.sky.sky_from` over the sphere query) and send it packed (`PGK1`, `web/lib/skypack.py`). The limiting magnitude is `?limit=`, 2 to 9, default 7.5. The packed list's `meta` carries the nebula and supernova remnant sprites and a report of how much of the sky is built (`built_to_pc`, `complete_above_sol`, `wanted_pc`, `filling`, the counts per tier). The API route is additive, so the API version stays at 12.
- The sprites are one record per nebula, remnant or cloud that is wider than 0.002 radian and either shines under the naked-eye limit or dims what lies behind it: `{id, kind, family, direction, distanceLy, angularRadiusRad, inside, extinctionMag, surfaceBrightness, visible, seen}`, the output of `nebula_sprites.sprite_data` with the direction taken from where the object was when its light left it (`galaxy/sky_sprites.py`, `db/sky_sprites.py`). A remnant whose explosion has not been seen is left out. This is the shape `skysprites.js` reads.
- The web layer keeps each sky on disk in the Galaxy Map's cache (`web/lib/skycache.py`), keyed by system, limit and a version key made from the galaxy's tile stamp, so it is rebuilt when stars are stored and not on an orbital update. It is dropped with the cache's generation, counted in `cache.max_mb`, and `skycache.clear()` is the call for the cache page's clear (PERF.134). The API answer is not put in the shared page cache.
- The wire format, the sprite records and the report are written down in `docs/design/interface/system-sky-wire-format.md`; stand-ins SPR1 to SPR3 and SKV1 to SKV3 are registered (6.179).

generation-output: unchanged
