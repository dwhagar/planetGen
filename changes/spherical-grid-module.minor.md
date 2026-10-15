### Added
- `planetgen/galaxy/spherical.py`, the pure-geometry module for the spherical sector grid: shells (4 pc deep, shell 0 the 4 pc core ball), bands, wedges (multiples of 6) and polar caps (6 wedges at the rim, chord radius 3.82 pc), with `sector_address_at`, `sector_position_pc`, exact cell volumes, face neighbours, the sector-local frame and a cap-aware `SectorCell` sampler. Shells 1 to 3 are a hand-checked table (shell 1's caps are cut at 37.5 degrees). Nothing uses it yet; the cylinder grid stays in place until the later steps of the Spherical Sector Coordinates work.

TODO items resolved: GEN.224 / Added: none
