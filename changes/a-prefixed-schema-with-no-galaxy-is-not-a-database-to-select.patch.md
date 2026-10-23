**Headline:** `?db=` naming a schema that shares the prefix but holds no galaxy now answers 404 instead of 500.

### Fixed
- **A stray `planetgen_*` schema made every route answer 500 when selected with `?db=` (API.32).** `resolve_database` accepted any schema of the prefix that was not the configured control schema, so a server holding a schema such as the real `planetgen_control` (when the configured control name is another) let `?db=planetgen_control` through, and the first missing table (`sectors`, `galaxy_containers`) was a 500. It now checks that the schema holds the galaxy's `sectors` table and otherwise answers the same 404 as an unknown name. `test_fuzz_web_routes.py` and `test_api_bad_input_is_a_clean_4xx` now pass on a machine that holds a `planetgen_control` database, which they failed on (39 failures in one gate).

TODO items resolved: API.32. Added: none.
