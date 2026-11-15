**Headline:** The Galaxy Map at galaxy scale no longer shows an even lattice of round brown blobs: a group of nebulae shows only where its block really holds gas, as strongly as the block is full of it.

### Fixed

- MAP.429: the nebula group sprites (`galaxynebulae.js`, MAP.176) were drawn for every block with any nebula, sized to its bounding sphere and as strong as the logarithm of the number of clouds, so the many tiny class N clouds made a sprite of every disk block. A node's `nebula_regions` now carry each block's stored `cover` (the share of it inside a nebula, `nebula_block_stats.cover`, the same share the Nebula cover colour reads), a group under the floor of 3 percent is not drawn, and the strength rises with the cover to its most at 10 percent. The count no longer sets the strength. The stored nodes are rebuilt once (`TILE_FORMAT` 13).

generation-output: unchanged
generation-affects: none
