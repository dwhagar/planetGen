**Headline:** The second layered-fixes cleanup: eleven more places where the same code was written twice (child-row inserts, the two compact-remnant inserts, the four class and body-type facets, the orbit-tree writers, the checksum marks, the batched reads of the fingerprint and checksum readers, seven API routes, the galaxy's current shape and seed, the settings tree writer) are one function each, with nothing the generator, the API or the pages produce changed.

### Changed
- **One child-row insert in the store.** `_insert_deposit_types` and `_insert_isotopes` are one `_insert_child_rows(conn, table, id_column, owner_id, rows, fields)`; `body_gases` and `body_class_scores` read through one `_body_child_map`.
- **One compact-remnant insert.** `insert_black_hole` and `insert_neutron_star` share `_insert_compact_remnant` (row, galactic motion, spin, name registration); each passes only its table, fields and kind.
- **One facet builder.** The four class and body-type facet functions in `db/query.py` share `_class_facet` and `_body_type_facet`.
- **One row insert for the orbit trees.** `orbit_trees.write_nodes` and `write_dark_members` insert through `_insert_rows`; each keeps its literal DELETE so the storage checksum's writer ledger still sees it mark the sector.
- **One batched mark in the storage checksum.** `mark_sectors` and `mark_systems` share `_mark_each`; `fingerprint.select_in` is the one batched `IN (...)` reader, and the checksum's reader calls it.
- **API routes written once.** The six routes about one kept body (belt page, belt member, irregular moons and moon, field members and member) share `_kept_body`; the job, page-job and page-job-cancel routes share `_job_answer`; the system and sector facilities routes share `_facilities_of`. Same status codes and messages.
- **The galaxy's current shape and seed are read in one place.** `shape_draw.current_shape_and_seed` serves `potential.current` and `enclosed_mass.current_enclosed_mass_kg`.
- **One settings tree writer.** `util/settings.set_path` (was `_set_path`) is what `settings_writer._put` points at.

generation-output: unchanged
generation-affects: none
