### Changed
- The hand-written `_migrate_vN_to_vM` steps (v8 to v61), their old-schema test fixtures and the tests of each step are removed; Alembic carries every migration from the v61 baseline. A database older than v61 is refused with `SchemaTooOldError` (run `update.sh` from an earlier checkout first). `PLANETGEN_MIGRATIONS_DIR` points the migrations elsewhere for tests (DB.11 complete).
