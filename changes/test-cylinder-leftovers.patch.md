**Headline:** The cylinder grid's last module and tests are gone; the tests save sectors at shell, band and wedge addresses.

### Changed
- `planetgen/galaxy/geometry.py` (the old cylinder grid, unused except for a copy of the binary-companion offset that `frame.py` also holds) and `tests/test_galaxy_geometry.py` are removed; the two importers use `galaxy.frame`.
- Tests that saved sectors with `ring_index`, `layer_index` and `ring_slot_index` keys (which the store silently ignored, leaving the sector without an address) now pass `shell_index`, `band_index` and `wedge_index`; so does the query-budget tier's clone of the sectors table. The planned-slots fuzz test checks the shell-grid position of each slot.

TODO items resolved: none (TEST.156 narrows to its last file checks) / Added: none
