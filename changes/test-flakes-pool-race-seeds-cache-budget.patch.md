**Headline:** Several threads asking for a database connection at once no longer leave extra connections open on the server.

### Fixed
- **A first burst of connections built one connection pool per thread (TEST.218).** `_get_pool` now builds a database's pool under a lock, so threads that are the first to ask share one pool instead of each opening a connection that the last pool stored would drop without closing.

### Changed
- **Tests: the object index budget test (TEST.273) keeps its records under its own Redis keys,** so another worker's database cannot take bytes off the tile cache budget between two reads.
- **Tests: the prevalence percentage tests (TEST.233) use a seed that is the same in every process,** so a failure repeats and the assert message prints the seed.
- **Tests: the control connections test (TEST.218) ends its pool and asserts it left no connection open.**
