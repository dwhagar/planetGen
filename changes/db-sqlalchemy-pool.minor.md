### Changed
- The database connection pool is now SQLAlchemy's `QueuePool` (5 kept, up to 10 open, a dead connection replaced when checked out) instead of DBUtils; `dbutils` is no longer a dependency. Id blocks, batching, statement timeouts, UTC sessions and the never-fail-over rule behave as before (DB.11, first step).
