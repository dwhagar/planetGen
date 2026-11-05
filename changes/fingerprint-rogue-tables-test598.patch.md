**Headline:** The galaxy fingerprint now reads the child rows of a rogue planet, so a rogue's air, layers and spectrum count when two builds are compared.

### Changed

- `db/fingerprint.py` listed neither the four child tables a rogue planet got with GEN.1111 (`rogue_planet_reflection_spectrum`, `rogue_planet_interior_layers`, `rogue_planet_gases`, `rogue_planet_life_class_scores`) nor `designations` (GEN.1138) and `chart_ledger` (GEN.1127), so `test_every_table_is_compared_or_left_out_on_purpose` was red. The four rogue tables are content like those of a planet: they are compared (keyed through `rogue_planet_id`) and sit in the `planet` section. `designations` (retired names and IDs, tied to row ids) and `chart_ledger` (when a sector was charted and by which run) are bookkeeping like the name registries, so they are left out.
- The pinned galaxy's sector digests change for the sectors that hold a rogue planet (12 in the golden galaxy), so `src/tests/golden/galaxy.json` must be rewritten in the same merge.

generation-output: unchanged
