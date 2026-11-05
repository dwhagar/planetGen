**Headline:** The galaxy fingerprint now reads the child rows of a rogue planet, so a rogue's air, layers and spectrum count when two builds are compared.

### Changed

- `db/fingerprint.py` listed neither the four child tables a rogue planet got with GEN.1111 (`rogue_planet_reflection_spectrum`, `rogue_planet_interior_layers`, `rogue_planet_gases`, `rogue_planet_life_class_scores`) nor `designations` (GEN.1138) and `chart_ledger` (GEN.1127), so `test_every_table_is_compared_or_left_out_on_purpose` was red. The four rogue tables are content like those of a planet: they are compared (keyed through `rogue_planet_id`) and sit in the `planet` section. `designations` (retired names and IDs, tied to row ids) and `chart_ledger` (when a sector was charted and by which run) are bookkeeping like the name registries, so they are left out.
- The pinned galaxy's sector digests change for the sectors that hold a rogue planet (12 in the golden galaxy) and the leaves of the sectors that hold a rogue change, so `src/tests/golden/galaxy.json` must be rewritten in the same merge.

- The fingerprint tree's leaf for a sector (`sectors.leaf_digest`) now holds the rogue planets filed in it (TEST.602). A rogue lands in a sector other than its donor's, often one nobody has charted, and was left out of every leaf: the 24 sectors made only for a landing had no leaf at all, and a charted sector's leaf did not change when a rogue landed in it. The leaf is now the sector's content part (the new column `sectors.content_leaf`, migration 0212) joined to a digest of its rogues with their child rows and moons, without their run-time IDs; a sector made only for a landing has a leaf of its rogues alone, and a sector with none keeps the content part as its leaf. Every writer of a rogue (a donor saved or saved again, a top-up, a sector emptied that held rogues born in it) makes the leaves of the sectors it touched again.

generation-output: unchanged
