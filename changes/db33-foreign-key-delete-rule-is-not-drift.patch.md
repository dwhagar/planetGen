**Headline:** The database check holds the rule that a foreign key's default delete rule is not drift on any supported server.

### Changed
- `check-db` and the model comparison already pass on a healthy MySQL 8.0.46 database (the models declare `ON DELETE RESTRICT` where schema.sql does, and Alembic reads a key MySQL 8.4 reports with no rule as the same). DB.33 is therefore closed by a regeneration that landed earlier; this adds the tests that keep it closed.
- `tests/test_db_check.py` rebuilds `fk_facilities_type` with no rule, with `ON DELETE RESTRICT` and with `ON DELETE NO ACTION`, and none is drift; with `ON DELETE CASCADE` it is.

TODO items resolved: DB.33. Added: none.
