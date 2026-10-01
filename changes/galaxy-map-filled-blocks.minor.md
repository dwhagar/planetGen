### Changed

- The Galaxy Map shows generated and not-yet-generated sectors through its blocks alone; the marker dots are gone. Unfilled space is see-through (more so where it's sparse), and a block grows more solid and warmer the more of its sectors are generated, fully solid once all are, so generated sectors can be found by zooming at every scale. Zoomed in to single sectors, a generated sector takes its real density's color and links to its page.
- Clicking a block shows its ring, layer and slot ranges, its exact sector count and how many are generated; a click prefers the nearest block holding generated sectors along the line of sight, and centers on them, so a double-click zooms toward them.
- `/api/galaxy/tiles` tiles carry a `filled` summary counting every placed sector in the tile (per sector, or per cell for large or crowded tiles), which the map sums into its blocks.
