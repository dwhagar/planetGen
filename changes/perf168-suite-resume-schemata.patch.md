**Headline:** The benchmark suite's `--resume` looks for its scratch database in a way the SQL portability test accepts.

### Changed
- `planetgen benchmark --suite --resume` finds the scratch database with a query on `information_schema.SCHEMATA` instead of `SHOW DATABASES`, whose `DATABASES` is a reserved word on MySQL 8.0 and 8.4. Nothing else about the resume changes.
