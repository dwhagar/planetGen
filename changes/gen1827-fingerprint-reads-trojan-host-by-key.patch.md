**Headline:** The fingerprint reads a trojan planet's giant by its key, so two orders of charting the same sectors read as one galaxy again.

### Fixed
- GEN.1827: `tests/test_region_order_independence.py` failed on main (`{'9 -1 53': {'planet': 'structure'}}`). The cause was the fingerprint, not the generator: `planets.lagrange_host_id` (GEN.1151) holds the id of the giant a trojan planet shares its orbit with and is no declared foreign key, so the fingerprint read the raw id, which is the id the row got in that run. `fingerprint.ROW_ID_COLUMNS` names such a column and the table it points into, and `_Reader.row_value` writes the key of that row (as it does for `parent_id`). A new test lists every `*_id` column that is no foreign key and requires each named one to be in the map.

generation-output: unchanged
