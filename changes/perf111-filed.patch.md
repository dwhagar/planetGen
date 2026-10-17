**Headline:** A page cache that never trims its no-database entries is filed, with the reproducer checked on main.

### Changed

- PERF.111 is filed (bug, found by Core 2 in the MAP.217 gate): `ResponseCache.put` without a database skips the eviction loop, so the cache passes `max_entries`; the reproducer gives 6 entries for a bound of 3 on main.
- TEST.213 gets a sighting with a second symptom (the system's scene never came in, under load 10); it passed in 124.68 s on an idle box.

TODO items resolved: none / Added: PERF.111
