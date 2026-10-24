**Headline:** The page cache no longer grows past its entry limit for pages kept without a database.

### Fixed
- **`ResponseCache.put` without a database never trimmed the cache (PERF.111).** Entries of the no-database kind were stored and never evicted, so the cache could hold more than `max_entries`. It now trims like the keyed branch; a test fills 20 entries into a limit of 5.
