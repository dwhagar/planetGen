**Headline:** When a Galaxy Map tile cannot list every star or point, it now thins the band it cuts evenly across space instead of keeping the brightest few wherever they cluster.

### Changed
- Galaxy Map: the pre-placed bright stars, the placed point phenomena and the scattered points of a tile are cut by one function (`_even_pick`, MAP.251). The halo and the satellite galaxies still go first (MAP.250); inside a tier the most massive and the brightest bands (half a power of ten of the key) are kept whole, and the band the cut falls in loses the same share of its objects everywhere, chosen by a fixed number every object carries, so a crowded corner no longer takes the whole list and the lists of finer tiles still hold what the coarser ones showed. The database reads four times the limit of each ranking so the cut band has candidates. The distance rule stays in the per-level budgets, and the generated stars (already shared out by sector) and the galaxy-wide sample are as they were.
- Docs: `docs/design/zoom-star-visibility.md` section 11 writes the rule.
- Tests: `test_even_drop.py` (the cut band thins by the same share in a crowded and a sparse region, no brighter band is dropped while a dimmer one is listed, longer lists keep a shorter one's objects, the tier order, order independence) and the brute-force bright-star box test and the MAP.154 weight test now expect whole bands above the cut.

TODO items resolved: MAP.251 / Added: none
