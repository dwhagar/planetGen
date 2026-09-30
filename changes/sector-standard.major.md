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
