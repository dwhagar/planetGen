**Headline:** The third layered-fixes cleanup: the same code written two to six times across the queue, the admin tools, the web pages, the maps, the galaxy model, the physics and the store is one function each (the process-alive check alone had six copies, two of them without the zombie fix the others carried), with nothing the generator, the API or the pages produce changed.

### Fixed
- **One answer to "is this process alive".** `util/process.pid_alive` replaces the six copies in the queue (`workers.process_exists`, `runrecord._pid_alive`, `stopall._alive`, `memorygate.alive`, `unitstuck.alive`) and the update guard (`updateguard.pid_alive`); the names stay as aliases for the tests that patch them. Two copies treated a zombie (a worker that ended but was not reaped) as alive and one treated a permission error as dead; every caller now gets the complete check.
- **The nebula field writer splits long statements.** `store.insert_nebula_field` writes through `_insert_literal_rows`, so its statements are capped by `LITERAL_STATEMENT_CHARS` like every other literal bulk insert; the rows written are the same.
- **Two stale imports and one dead variable** are gone (`store`, `run_plan`, `query`, `system`).

### Changed
- **Queue plumbing written once.** `wrapper._await_plan` and `_await_conclusion` share `_await`; `unitstuck.records` and `unitwatch.read` share `runrecord.unit_json_hash`; `updateguard`'s two Redis sources share `_RedisSource`.
- **Admin account scripts share their frame.** `cli/admin_account_command.py` holds the argument parsing, control-database connection and admin lookup that `forget_devices`, `reset_two_factor` and `accounts` each carried; each script keeps only its action and sentence.
- **The web and the maps.** `tilecache._write_json` and `_write_bytes` share `_write_atomically`; the three map data routes (`galaxy_stage`, `galaxy_window`, `galaxy_near`) share `_tile_json`; `warmup`'s two node warmers share `_warm_node_keys`; `jobs`' three log readers share `_known_job_dir`; the system map's three planet and moon marker loops share `_body_marker_svgs`; the belt and field member pages share `bodyview.member_page_ref`.
- **The galaxy model.** `gravity_cache.prune` and `clear` share `_remove`; `drill.block_key_parts` serves `parse_drill_key` and `window.parse_window_key`; `spherical`'s three sphere walks share `_check_sphere`; `sky_band`'s two cube-face lookups share `_face_coordinates`; `query`'s scattered-point and scattered-nebula box reads share `_scatter_address_clause` and `_scatter_box_params`.
- **The physics.** `giant_interior`'s two table interpolations share `_bracket`; `king.TableMethods` holds the three profile-table methods `KingTable` and `component_profiles.ProfileTable` each had; `position`'s two sector-edge checks are `_sector_edge`; `co2_far_wing`'s two far-wing cross sections share `_layer_cross_sections` (each keeps only its tail).
- **The store and the plan.** `_reserve_id_block` and `_reserve_counter` share `_reserve_range` (one connection, one upsert, one commit); `run_plan`'s `clean_nebula_field`, `redraw_nebula_field` and `clear_scatter` use `_refusals_as_errors` instead of three copies of its body.
- **Settings.** `settings_writer`'s plan and update plan share `_load_document` and `_note_field_errors`; the `system` and `sector` subcommands share `_add_flavor_arguments`.

generation-output: unchanged
generation-affects: none
