**Headline:** Pages that name a database (`?db=`) no longer scan every schema of the MySQL server on each request: the Galaxy Map's near-sector read took 20 to 26 seconds on a server with 147 schemas and takes about one second now.

### Changed

- `db.store.resolve_database` checked that a requested schema exists by calling `list_databases`, which reads the size and modification time of every schema with the prefix (one `information_schema.tables` scan each) and throws the numbers away. It now asks for that one name with one query on `information_schema.schemata` (`_is_listed`), with the listing's rules unchanged: the prefix, not the control schema, and the name exactly as written.
- Why it mattered: a Galaxy Map view fetches a dozen stages, each resolving its schema, so one `/galaxy/near` call made 59 resolutions and 8,383 queries. On a test server that had gathered 147 schemas that was 20 to 26 seconds a call, and the whole-galaxy wheel-in browser test (TEST.696) failed because the sector scene arrived long after the camera had flown past the star.
- The `list_databases` page of the admin screens is unchanged (it does want the sizes).

generation-output: unchanged
