**Headline:** A test now proves that changing any row of a rogue planet's four child tables changes the digest of the sector it is filed in.

### Changed

- `test_a_changed_row_of_a_rogue_planets_child_tables_changes_the_digest_of_the_sector_it_is_filed_in` (TEST.647) changes one row of each of `rogue_planet_reflection_spectrum`, `rogue_planet_interior_layers`, `rogue_planet_gases` and `rogue_planet_life_class_scores` and sees both `fingerprint.sector_digest` and the leaf's `fingerprint_tree.resident_digest` change, and come back when the row does. The tables were already in the fingerprint (TEST.598) and the leaf (TEST.602) when train 12 merged, so no code, epoch or golden changes.

generation-output: unchanged
