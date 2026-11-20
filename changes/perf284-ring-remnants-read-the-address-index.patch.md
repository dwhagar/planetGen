**Headline:** The special phenomena stage no longer takes hours on a big galaxy, and it says how many ring sectors it has done and how long is left (PERF.284, OPS.289).

### Fixed

- The ring remnants piece of "Write the special phenomena" reads the light of the cells each chunk reaches by sector address, and the scatter has taken the address index off the phenomenon table for its inserts, so every chunk of 200 sectors scanned the whole table (2 to 3 s each, mysqld at one full core, on Lane 15's run: 2h 33m and counting). The piece now builds that one index first when it is missing; the index stage after it finds it present. On a stand-in table of 3 million rows a 500-cell read took 1.5 s without the index and 0.02 s with it, and the index took 5 s to build.

### Added

- The ring remnants piece says `Neutron stars and black holes of the backfill rings: 12,000 of 91,500 done (13%), 40 a second, 5m 3s elapsed, about 33m 20s left.` at the log cadence and as the unit's second bar (the real count of sectors, OPS.289), where it showed no count.
