**Headline:** The two 20 percent water rows of the ocean-world table in the design notes are recomputed from the built code, and a test now holds every row of that table to it.

### Fixed
- **Ocean-world table (DOC.79).** The x = 0.2 rows of design part C 2.4 (1 and 5 Earth masses) said 1635 km at 26.0 GPa and 2519 km at 114 GPa, which did not follow from the table's own radii. They now read 1294 km at 19.2 GPa and 2024 km at 84.7 GPa, with the gravity and the surface-pressure estimate rechecked, and the hydrosphere note no longer calls them a non-target. `tests/test_water_layer.py` reads the document table and compares each row with `physics/water_layer.py`. Documents and a test only: no generated value changes.
