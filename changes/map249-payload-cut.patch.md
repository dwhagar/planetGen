**Headline:** Zooming the Galaxy Map asks the server for about half as many tiles, so it lags far less.

### Changed
- **The Galaxy Map's distance cut is coarser and narrower (MAP.249).** A zoom step used to ask for about 450 tiles (up to 3.2 MB packed, tens of seconds cold on the server). The cut now splits a tile while its edge is more than 3 times its distance (was 2) and keeps a 50 degree half-angle cone (was 70), so a view needs 213 to 271 tiles. The per-tile budgets are unchanged, so each tile is as full as before; the view simply holds fewer tiles, and `GALAXY_VIEW_MAX_STARS` falls from 70,000 to 45,000. The browser twins (`galaxytilecut.js`, `starfade.js`) carry the same numbers.
- **Stars drawn at each zoom level, before and after, on the scaled probe galaxy (the real galaxy was not available, so these are indicative):** 100,000 pc 3,000 to 3,900 (+30%); 30,000 pc 5,126 to 4,813 (-6%); 10,000 pc 19,686 to 9,392 (-52%); 3,000 pc 18,229 to 13,268 (-27%); 1,000 pc 17,511 to 12,075 (-31%); 100 pc 16,154 to 9,608 (-41%); 30 pc 15,969 to 9,325 (-42%). Tiles fell from 414 to 478 down to 213 to 271, packed bytes by 10% to 48% (3.20 to 1.88 MB at 1,000 pc), and the cold server build at 1,000 pc from about 62 s to about 37 s. `docs/design/fly-through-view-distance.md` lists the knobs for tuning the star count back up.

### Added
- Nothing.

TODO items resolved: MAP.249 (payload cut; the cache-warming half, PERF.130 to 132, follows in its own change) / Added: none
