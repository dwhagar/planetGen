**Headline:** A new design note shows how to ask the database for everything inside a sphere with a single cell key and a covering index, measured on real MySQL 8 against the box, the spatial index and walking the sectors.

### Changed
- Added `docs/design/spherical-region-queries.md`: the cell key (Morton, 0.25 pc cubes, 262 kpc box, a generated column), the sphere-to-key-runs cover, the `VALUES` join statement, counts, nearest-first and the dense core, index maintenance cost, what was not measured, and eight proposed TODO items. No code changed.

TODO items resolved: none / Added: none (the proposed items are in the note's section 9)
