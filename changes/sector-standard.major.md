### Changed
- **One sector standard: 4 parsecs.** Every sector, in the galaxy or
  standalone, is now 4 pc (about 13.05 ly) on a side instead of 11.5 ly. At
  the default Milky Way shape the galaxy then reaches about 50,000 ly, the
  real disk's radius. See `docs/design/galaxy-coordinate-system.md`,
  "Sector size".
- **The galaxy is a stack of aligned layers.** Ring `i` now holds
  `round(2π(i + ½))` slots (3, 9, 16, 22, …) instead of a multiple of 4,
  the same on every layer, so sectors line up in vertical columns. The
  skeleton is stored per layer: each layer, from the top of the galaxy to
  the bottom, runs from ring 0 out to the last ring that still expects a
  star per sector (`galaxy_layer`, replacing `galaxy_ring_band`).
- **Generation never lands outside the galaxy.** `generate.py plan` also
  stores each ring's column bound (`galaxy_column`: the highest and lowest
  layer it reaches), and every `generate.py galaxy` mode checks its address
  against the stored outline before generating anything. Explicit
  `--density` or `--num-systems` no longer skip that check; an address
  outside is refused with the reason, and a neighborhood near the edge
  leaves out the sectors past it. `generate.py galaxy` now refuses to run
  before `generate.py plan`.
- **Random starts are drawn from the real outline.** A random start picks a
  uniformly random sector inside the planned galaxy instead of from a fixed
  15,000 pc disk 2,000 pc tall, so the starting neighborhood can no longer
  land outside the galaxy. `--max-ring` defaults to the galaxy's own edge.
- The Galaxy Map's prisms use the generator's own one-star-per-sector
  threshold, so fully zoomed in their outline is exactly the galaxy's
  layers.
- `generate.py plan` no longer takes `--edge-ly` or `--empty-streak-to-stop`;
  the edge is always the standard, and the build takes a few milliseconds.

### Removed
- **Upgrading deletes every galaxy-placed sector again.** Schema v33's
  migration deletes each placed sector with its systems and phenomena, since
  nearly every address moves, and rebuilds the skeleton from the stored
  shape at 4 pc. Sectors that were never placed in the galaxy are kept.
  Regenerate the galaxy afterwards (`generate.py galaxy`), and take a
  backup first if you want the old data.
