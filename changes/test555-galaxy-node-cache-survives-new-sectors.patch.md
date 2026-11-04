**Headline:** Adding sectors no longer makes the Galaxy Map throw away every node it has cached; it now fetches only the sectors that changed.

### Changed
- **Fix (TEST.555 part 2).** `db.query.galaxy_shape_inputs` read the whole `galaxy_shape` row, including `density_ratio_avg` and `density_ratio_samples`, the decaying average that every sector fill moves (PERF.11). Those two columns feed `tile_base`, so each new sector changed it, and `/api/galaxy/changes` then answered "everything" with an empty history instead of the changed sectors. A browser that had loaded the map reloaded the whole galaxy after any CLI `sectors` run. Only the fill estimate in `db.stats` reads the two columns, and no node, tile or density cloud depends on them, so they are left out of the shape inputs. `test_galaxy_map_and_nodes_show_sectors_the_cli_placed` passes again.
generation-output: unchanged
