**Headline:** The block light totals are now the last stage of a galaxy run, and their bar moves through all three parts.

### Changed
- **Block lights last (PERF.251).** `galaxy_stages` puts "Make the block light totals" after the settling and the run blocks instead of after the nebula field. The totals add up only stars and phenomena that are not built yet, so they were redone for nothing as the backfill, population and settling built them; made last they are made once from what is left. The stage list, `run_galaxy`'s followers, `commands`' stage numbers and the design docs (`redis-and-the-work-queue.md`, `command-inventory.md`) follow the new order.
- **The uncharted placeholder rows come after the backfill.** The totals stage writes the uncharted sector rows, which the bright-star backfill used to count as sectors the run had just made; with the totals last the backfill opens on the run's real sectors only.
- **The bar moves all the way.** "Making the block light totals" counts 100 steps in three parts (the units 70, the levels above the units 15, the uncharted rows the rest) and says the share done in its log line, so it no longer sits at one value through the longest part.
- **Tests.** `test_stages.py` has the new stage numbers and a test that the block lights stage is last; `test_block_lights.py` checks that the bar moves in all three parts.
