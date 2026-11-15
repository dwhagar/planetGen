**Headline:** The sky draw of a New galaxy run no longer waits in line: it runs as 16 units side by side, so one stage takes minutes of waiting instead of hours.

### Fixed

- OPS.253: the sky tiers of a bright-star backfill were drawn one observer cell at a time behind one lock, so a stage with several cells made no progress for hours while the rest of the workers waited. The draw is now cut into 16 units, each drawing its own share of the rows of every sphere (every cell of a sphere is in exactly one share), so nothing waits for anything; the spheres are recorded in the ledger when the stage ends. The stars, the cells' levels and the ledger come out the same as the one-at-a-time draw (census test), and a unit that stops goes on from its last batch. `sky_tier_draws` gains `part_index` and `part_count` (migration 0260).
