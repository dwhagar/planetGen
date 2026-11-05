**Headline:** Every scattered star and scattered phenomenon now has an Object ID from the moment it is placed, and a system or object built from it keeps that ID.

### Added
- Migration 0212 adds a `uid` column with a unique index to `bright_stars` and `phenomenon_scatter`. `db/scatter_uids.py` gives each new scatter row its birth sector, a serial that is its rank in that sector (brightest first, stars from 0, phenomena from the upper half of the serial space) and body number 0. A later write for the sector, such as a backfill of dimmer stars, continues the numbers. The work is set-based in the database, one window query and three joined updates for each batch the scatter writes. About 24 bytes a row: 0.65 GB at 26.9 million rows.
- `tests/test_gen1254_scatter_uid.py`: ranks brightest first in the right sector, a later write continues the numbers, stars and phenomena never share a serial, the same rows give the same IDs in any insert order, and a sector's fill gives the system built round a star the star's ID and starts its own serials past the stars'.

### Changed
- A system the sector fill builds round a scatter star takes the star's ID, a built scattered phenomenon takes its row's, and the fill's own systems start past the stars' serials. The route that builds one uncharted star's system does the same. The first system of a sector that has scatter stars no longer has serial 0.
- `clear_bright_stars` and `clear_phenomenon_scatter` forget the sector counters, so a scatter made again starts from serial 0.

generation-output: changed
TODO items resolved: GEN.1254 (the epoch and golden ride train 11: only the IDs of systems and objects in sectors that have scatter stars move)
