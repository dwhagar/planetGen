**Headline:** The bright-star backfill draws the stars of a sparse block of cells with one count for the block instead of a draw in every cell.

### Changed
- **The backfill picks its draw path per block (GEN.309).** A block (a band's 256 cells) whose largest per-cell rate in a luminosity or mass band is under 0.05 stars now draws that band's stars block first: one Poisson count for the block, each candidate kept with the chance its cell's own rate asks, and each kept star on a stream of its own (`bright_stars.block_candidates`, built in GEN.308). Denser blocks draw cell by cell as before. The choice depends on the block and the band alone, so a sector gets the same stars whichever other sectors a run reaches, in whatever order, split between workers or restricted to a region. The stars a seed gives in the sparse blocks (the outer disk, the halo, the far edge) are different draws of the same distribution, so existing galaxies keep what they hold and new ones differ there.
generation-output: changed

TODO items resolved: GEN.309 (the part that picks the path; GEN.40 owns weeding sectors out before the backfill). Added: none.
