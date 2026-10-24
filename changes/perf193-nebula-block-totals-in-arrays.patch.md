**Headline:** Building the nebula blocks of the Galaxy Map needs about a seventh of the memory.

### Changed
- **The nebula block totals are arrays, and a level's rows are made as they are written (PERF.193).** `fine_cover` and `roll_up` kept the cubic parsecs of nebula volume in each block in Python dictionaries (several million entries at the finest level, each rolled up through five levels), and `level_rows` then made a dictionary for every block of the level at once. They now keep sorted arrays of the packed block keys and their totals, add a level into the next with one array pass, and `LevelRows.blocks` makes a block's row only when it is read, so the writer holds a batch of rows instead of the level. A test galaxy of 100,000 nebulae (7.7 million block rows) peaked at 6.0 GB before and 0.8 GB after, with the same rows, and took 127 s instead of 232 s. The numbers written are unchanged.

TODO items resolved: PERF.193. Added: none.
