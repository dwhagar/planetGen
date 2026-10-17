**Headline:** The benchmark records the buffer pool, table sizes and worker start-up cost, and the test suite is quieter and uses the shell, band and wedge words.

### Changed
- **The PERF.31 benchmark records its setting (PERF.47).** A run now stores the InnoDB buffer pool size, the size of each table and how long a worker takes to start, and prints them in the setup lines, so two runs on different machines can be compared.
- **A multi-row literal INSERT is built in one place (PERF.48).** The bright-stars and phenomenon-scatter inserts share one builder that splits a statement at a size limit; a numeric-only row no longer goes through the driver's escaping.
- **Browser-test failures are no longer buried in the request log (TEST.187).** The test server's request log is set to WARNING.
- **The marching-cubes call no longer floods the log with a NumPy 2.5 deprecation (TEST.188).** The warning comes from inside scikit-image, so the nebula mesh builder filters it around the call; remove the filter when scikit-image stops raising it.
- **Tests that named ring, layer and slot now name shell, band and wedge (TEST.156).** The migration tests that describe the old schema keep the old words.
