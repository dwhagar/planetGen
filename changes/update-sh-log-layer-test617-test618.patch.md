**Headline:** The update script prints its settling wait through the log layer, and the tile cache measuring script names a path that exists.

### Changed

- `update.sh` printed "The update steps are finished: waiting ... seconds before reloading the site." with a bare `echo`, the one line of its own that did not go through `scripts/logrecord.sh`; it now uses `log_info` like the others, so it is timestamped and recorded (`test_update_sh_prints_its_own_lines_through_the_layer`, TEST.617).
- The docstring of `scripts/measure_tile_cache.py` said the sector-block node lives in `db/nodes.py`; it now gives the path `src/planetgen/db/nodes.py` (`test_scripts_only_name_repo_paths_that_exist`, TEST.618). The script itself is unchanged.

generation-output: unchanged
