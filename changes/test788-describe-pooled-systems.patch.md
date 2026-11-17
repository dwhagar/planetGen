**Headline:** `test_describe` makes the systems it reads once per process (and the compact-remnant sweeps in parallel) instead of generating the same systems again in each test; every test checks the same systems and keeps its assertions.

### Changed
- The systems of the file's seven variants (the first eight of each variant's seeded run, which `_systems(n)` already read as its first n) are made once per process, by forked workers at the same time, and shared by the tests that read them. The documents of a system do not change when it is described again, after other systems, or after a `pickle` round trip (checked on 14 systems); no test changes a system it is given.
- The three compact-remnant tests (6, 12 and 60 seeded systems) make their systems in forked workers, each on its own seed; an assertion inside a worker reaches the test, and the set of pulsar types is gathered as before.

generation-output: unchanged
