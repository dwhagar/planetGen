**Headline:** The backfill of a sector's neighbors is proved safe to start as soon as that sector is saved, while other sectors are still being made.

### Changed
- OPS.242: `independence.json` lists a third by-address edge, `scatter-backfill`: the unit of one source sector (its step targets, each at the mass of its step from that source alone) may start when the source and every sector of its list that the run saves are saved. `tests/backfill_edge_support.py` runs it on a toy galaxy in plan order and by address (random orders on one thread, and on 15 threads) and compares the stored rows by content; a sector reached from two sources is the same in either order, and a unit released before the saves stores other bytes (the control). No generation output changes.
- The mass limit and the sky tier ledger the draw reads are written by earlier stages, so they stay stage edges (Lane 13).
