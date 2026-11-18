**Headline:** The band-table cache's lock poll interval carries its tag, so the constants census passes again.

### Fixed

- PERF.281: `band_table_cache.LOCK_POLL_SECONDS` says in its docstring that it is machinery, not physics, and defined where it stands, which `test_constants_completeness` (TEST.715) asks of a numeric constant in the generator's own modules. No value or behavior changes.
