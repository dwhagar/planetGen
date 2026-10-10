### Changed
- CI: the `linux-update` job's "needs migrating" and "failed migration" steps rebuild the database from the v61 baseline fixture (the oldest schema `update.sh` upgrades from) instead of faking v48, which the code has refused since the Alembic cleanup. The old steps failed with "database is at schema v48, older than v61".
