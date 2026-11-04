**Headline:** The Galaxy Map has the density its nebulae will be drawn from, shaped by the galaxy's own gas and star-forming density: bulbous outlines with soft, wispy edges, a dense core, and more filaments the closer you zoom.

### Added
- `galaxy.nebula_density.density_at(galaxy_seed, shape, points_pc, field, detail)` (MAP.365, the shape side): a value in 0..1 at any point in the galaxy, for `"gas"` (dark and molecular clouds) or `"young"` (emission and reflection nebulae). It is a pure function of the position, the seed and the stored shape, so every node and sector reads the same value at the same point and the shapes meet without a seam; the object ID is never read. A point is inside a nebula at `FORM_LEVEL` (0.35) or more.
- `"gas"` is anchored on the seeded cloud field, leans with the model gas (the arms and the thin layer), and a superbubble clears it. `"young"` follows the young-star density. A coarse warp makes the outline bulbous at the coarsest level; each further octave of gradient noise (`wavelength_pc(detail)`, 64 pc halving up to `MAX_DETAIL` 8) adds filaments at the edge and leaves a dense core, and changes a value by at most `octave_bound`, so the coarse blobs stay where they were as the map zooms in.
- Wiring register line for the module until the node response and the drawing (MAP.365 parts b to e, Core 2) read it.

generation-output: unchanged
TODO items resolved: MAP.365 (part a, the shape field; the node response, drawing, tint and blend are Core 2's)
