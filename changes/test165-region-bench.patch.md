**Headline:** The region query is timed on a table of three million real star system rows, and the design note records the figures.

### Changed

- TEST.165: `scripts/bench_region_queries.py` grows a scratch database of `star_systems` rows (real generated rows, real sector addresses, real cell keys) and times the cell key statement, the box over indexed x, y, z and the sector walk at 1, 5, 10, 25 and 50 pc, warm and cold, checking that the three return the same systems. The results are in `docs/design/galaxy/spherical-region-queries.md` section 8.1: the key statement is the fastest at every radius (by only 1.3 times over the box at 50 pc), and it is 3 to 10 times slower than the research's 42 ms at 50 pc (427 ms counting 62,336 systems). No stored value changes.
