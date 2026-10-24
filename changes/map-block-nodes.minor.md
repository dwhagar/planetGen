**Headline:** The Galaxy Map's data can now be read as sector-block nodes, with the light of every block stored, the first step to one continuous map that has no zoom stages.

### Added
- **Sector-block nodes (MAP.260).** `GET /api/galaxy/nodes` and the site's `/galaxy/nodes` answer by block of the ladder (a sector, or a block of 3 to 177,147 shells): the brightest stars, generated stars, point phenomena and clouds inside the block, the sums of what is stored about it (its generated sectors and the light of the stars and phenomena scattered in it), and for a level-3 block or a sector the sectors placed in it. A block lists every object of its parent that lies inside it. The site caches each node on disk beside the tiles; a change to a sector drops the nodes of its chain, and the browser is told which, as it is for tiles.
- **Warm-up builds the nodes.** `planetgen warm` has a new family, `galaxy-nodes` (command `warm-galaxy-nodes`), in the `built-levels` phase: every populated block of level 243 and above.

### Changed
- **Block light totals (MAP.259).** New table `galaxy_block_lights` (schema v147, migration 0147): for each block of the eleven levels above a sector, the count, summed luminosity, centre of light and mean temperature of the pre-placed bright stars and the count, luminosity and centre of the scattered phenomena that nothing was built from. A scatter, a clear and the map-luminosity fill make it at their end; the migration fills it from the stored scatter, which can take a while on a big galaxy.
- **The block totals follow a `galaxy` run.** A run now brings the sectors it saved into `galaxy_block_stats` and `galaxy_block_lights` once at its end. Before, only a deletion refreshed `galaxy_block_stats`, so a block's stage and window figures stayed behind what a run generated until the table was made again.
- **`planetgen fill-block-totals`** makes both tables again (it can be run again and stopped).
