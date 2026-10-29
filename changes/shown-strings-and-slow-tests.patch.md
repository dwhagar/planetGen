**Headline:** Text the program shows no longer carries internal item codes, and the red or too-slow tests of the last two days are repaired.

### Fixed
- **Item codes (GEN.826, GEN.703, DB.48, GEN.441, DB.46) were shown in `--help`, the dataset catalogue, the repair report, the property sheet's note and the configuration page's option table.** They are gone from the text people read (TEST.434).
- **The Galaxy Map's browser tests expected the old map.** The habitable-world layer is on by default (MAP.303) and crowded marks merge (MAP.305), and a planet's habitability is a panel, not a collapsed group (UX.215); the tests say so (TEST.434, TEST.441).
- **The size-estimate tests ran the same 12-sector fill three times, past the test timeout on a busy machine; the phenomena-scatter tests filled a sector next to the core down to the backfill's deepest floor (an hour) and read a log line API.33 had reworded; the mass-backfill same-seed test ran eight whole draws.** The fill is made once per file, the scatter fills skip the deepest floor as the other tests do, and the same-seed test draws over three cells (TEST.422, TEST.441).
- **The albedo-migration tests expected the circular light.** The climate goes by the mean light of the orbit since GEN.727, so one test expects that and the other builds a galaxy that knew only the circle, as a v141 galaxy did (TEST.441).
- **A fresh test machine lacked the white-dwarf data and gave other stars.** The test preflight fetches it (TEST.441).
- **A rogue planet's ring optical depth declared the unit `tau`, which is no unit.** It is a plain ratio (TEST.425).
- **`galaxy/named_regions.py` took the distance of the anchor sector with `math.dist`.** It uses the plain-formula `detmath.dist`, which gives the same answer on every Python and library (TEST.388).

TODO items resolved: TEST.388 (named_regions), TEST.422, TEST.434 (parts 1 and 2), TEST.441 (the albedo and mass-backfill tests, and the data at setup).
