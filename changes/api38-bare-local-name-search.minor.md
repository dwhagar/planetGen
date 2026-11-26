**Headline:** A local name typed without its sector now finds the objects that have it, each with its sector, nearest sectors first.

### Added

- `GET /api/lookup?text=Zaidonr+Jaisonn[&near=<sector>]`: a bare local name (GEN.1255) names an object in every sector at once, so the lookup undoes the scramble of each sector in turn and reads the Object ID it gives. One match opens as an Object ID typed in full does (`current`, with its `sector`); several come back as `matches`, each with its sector; none is `unknown` and says how many sectors were searched. Every answer carries `searched`, `of` and `complete`.
- A search that is not done inside the request goes on as a job and answers `202 searching` with the matches of the nearest sectors so far. `GET /api/lookup/jobs/<id>` reports `searched` of `of` sectors with a progress line (`Searched 1,200 of 480,000 sectors in 1s`), the time left in clock form and the finished answer. The search never stops short without saying so.
- The sectors searched are every sector row, every cell of the scatter tables and every birth sector of the standalone phenomena, the one nearest `near` (the sector the person is in) first.

### Changed

- The scramble of a sector's serials is found 49 times faster (0.460 ms to 0.009 ms a sector, the same parameters for every sector, held by a test), which makes a pass over every sector cheap. No name changes.

generation-output: unchanged
