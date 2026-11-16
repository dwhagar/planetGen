**Headline:** the end of the phenomena scatter is now three commands of its own, so each can be seen, stopped and carried on by itself and the generated galaxy is the same.

### Changed

- OPS.239, part A: the `finish-phenomena` command is gone. `write-special-phenomena` (stage `phenomena-special`) draws the special phenomena and writes them with their classes in one transaction (the classes were added by the record piece before; a write that is repeated adds them once). `build-phenomenon-indexes` (stage `phenomena-index`) builds the table's indexes, a piece to an index, and stamps the hypervelocity stars' epoch. `record-phenomenon-scatter` (stage `phenomena-record`) records the scatter's seed, mass cut and epoch. `needs.py` declares the chain (`phenomena-shells`, `phenomena-special`, `phenomena-indexed`, `phenomena`), and a galaxy run counts two stages more.
