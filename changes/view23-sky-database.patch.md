**Headline:** The sky catalogue can now be built straight from the stars the database stores.

### Changed

- `galaxy/sky.py` (VIEW.23, the database entry point): `sky_from_database(conn, observer, at_unix, limit)` runs the sky over the sphere query of VIEW.22 (`db.sky_sphere.stars_in_sphere`) in shells of distance, and its report says how far the stored stars are complete to the limit (`built_to_pc`) and the largest level of the cells reached (`complete_above_sol`), the "built to magnitude X within Y pc" of the plan. `sphere_source` can now collect each shell's cell levels. Nothing reads it yet: the sky route (VIEW.31) does. generation-output: unchanged
