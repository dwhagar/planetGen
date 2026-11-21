**Headline:** Every rocky world with an air of 500 Pa or more now keeps the record of its air history in a new column, `air_history` (the history is not read back yet)

### Added

- GEN.1745 slice 4b: the life pass (`life._history_of`) makes the record of `physics/air_record.py` at its end, from the finished body, the air as the life found it and the air its oxygen's rise left (the build record keeps both, with the age of the rise), and the body keeps it as `air_history` (JSON text, `MEDIUMTEXT`, null for a giant, an airless world or an air under 500 Pa). The column is on `planets`, `moons` and `rogue_planets`, with a migration, and is stored and loaded with the rest of the body; an edit that keeps the stored life keeps the record. Nothing reads it back yet, so every other number of a world is the one it was (checked on three seeds, with and without the record).
- Cost: the build of a system takes 11.6 s without the record and 38.8 s with it (+234 percent), because the history solves the air in full at 12 ages of the star's light and heat for each of the 20 percent of bodies that get a record (53 of 263 in the sample; median record 5.9 KB, 90th percentile 20 KB, largest 44 KB).

generation-output: changed
generation-affects: fill
