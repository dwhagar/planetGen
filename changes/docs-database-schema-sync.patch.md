**Headline:** The database format guide now matches the v91 schema and the Alembic migrations.

### Changed
- `docs/database-schema.md` says the schema is at v91 and moves by Alembic revisions, adds the tables it was missing (`phenomenon_scatter`, `phenomenon_scatter_classes`, `sector_system_counts`, `object_timing_stats` and `object_timing_samples`, `id_blocks`, `id_counters`, `nebula_shape_balls`, `generation_run_arguments`, `alembic_version`) and the columns it lacked, describes the control schema's 18 tables at v14, and trims the planned-changes list to what is still open.
