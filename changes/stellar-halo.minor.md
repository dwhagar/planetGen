**Headline:** The galaxy has a measured stellar halo instead of a flat density floor.

### Changed
- The density model has a stellar halo (GEN.234): a broken power law on a flattened spheroid, 1.5e9 solar masses for the Milky Way (1.25% of its systems), inner slope 2.5, outer slope 4.5, break at 25 kpc, flat inside 500 pc. Its radii scale with the disk and its amplitude is solved like `k_norm`. It is 0.19% of the local density at the Sun and falls to 3e-9 at 262 kpc, so it adds no cell to the disk's outline.
- The flat floor (`MIN_RELATIVE_DENSITY`, 0.1% of the local density everywhere) is gone, with no setting: it was half the real halo at the Sun and 200 to 340,000 times too much at 50 to 262 kpc.
- Two new stellar populations, `halo` and `cluster` (10 to 13 Gyr), join the young, intermediate, old and bulge ones in the density split, the object-first scatter's bounds and the bright-star fractions. The halo's stars are the halo term's; the cluster's density is 0 until the cluster components arrive (GEN.237). Migration 0094 lets `bright_stars.population` hold them.
- The default galaxy's edge moves in by two shells (3,764 to 3,762) because the Sun's density of 1.0 now counts the halo.
- The Galaxy Map's density mirror carries the halo term. The normalisation to the local density (about 2.7 times higher) is kept as `tuning.HALO_LOCAL_DENSITY_FRACTION`.

TODO items resolved: GEN.234. Added: none.
