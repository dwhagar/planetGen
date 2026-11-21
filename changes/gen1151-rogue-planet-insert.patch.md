**Headline:** A rogue planet can be saved again: the trojan planet link columns are not its columns.

### Fixed
- GEN.1151: `lagrange_point` and `lagrange_host_id` are columns of `planets` only (a rogue planet has no giant to share an orbit with), so they join `ROGUE_HOST_ONLY_COLUMNS`; the rogue planet insert no longer names them (`Unknown column 'lagrange_point'` failed every sector unit of the pinned galaxy). A test saves a rogue planet and a moon of it.

generation-output: unchanged
