**Headline:** The bright-star and mass passes now reach the whole grid, out to 262 kpc, so the halo has its bright stars beyond the galaxy's outline.

### Added
- **The scatter reaches past the outline (GEN.322).** A shell's bands beyond the outline are drawn in groups of 64 bands, each bounded in one piece, and the faint shells far out (from about 50 kpc) are drawn together as one task, so a pass to the reach takes about a minute on one process, not hours. The stars inside the outline are exactly the rows they were. At the defaults it adds about 120,000 stars (halo giants down to the halo floor of 500 solar luminosities) beyond the outline. A galaxy scattered before this has none: `planetgen plan --redo-scatters mass luminosity` draws them.
- `tuning.SCATTER_COVERS_REACH` switches it on (the tests turn it off outside the reach tests). `band_majorants` is split into `band_majorants` and `extent_majorants`, which takes any region's nearest point and outer radius, so GEN.134's star-formation factor is bounded over each group's radius range. On the default galaxy it places about 122,000 stars beyond the outline.
- A galaxy scattered before this has no stars beyond the outline: `planetgen plan --redo-scatters mass luminosity` draws them.

generation-output: changed

TODO items resolved: GEN.322. Added: none.
