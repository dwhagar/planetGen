**Headline:** Migrations have helpers for slow changes on big tables, so the server refuses instead of silently copying.

### Added
- `planetgen/db/migrations/helpers.py`: `add_column_instant`, `drop_column_instant`, `add_index_online`, `drop_index_online`, `batched_update` (pages by key, a commit per page, progress on the migration's bar), the `exists_*` checks, `lock_wait_timeout`, and a check that refuses to start behind a long-open transaction. Each names its `ALGORITHM`, and each skips work already done so a stopped revision can be run again.
