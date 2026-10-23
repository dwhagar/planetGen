**Headline:** A sector holds the same stars and systems whichever order a galaxy's sectors are generated in.

### Changed
- **The first fill takes a sector to the backfill's deepest floor (GEN.633).** A `galaxy` run backfills dimmer bright stars around the sectors it generates, but that backfill skipped a sector that was already filled, so a sector generated before the backfill reached it held fewer systems than one the backfill reached first. The first fill of a sector now draws its bright stars down to the backfill's deepest floor itself, so no later pass adds stars and a sector comes out the same in either order (GEN.57). Sectors in a capped cell (GEN.281) are left as they were.

generation-output: changed
