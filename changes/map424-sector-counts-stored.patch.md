**Headline:** the Galaxy Map's sector list no longer times out on a big galaxy: each sector's system count is a stored count made in the background.

### Changed

- MAP.424: `GET /api/galaxy/sectors` counted every system in the galaxy on every request, which passed the 10 second web limit once the galaxy held millions of systems (the Galaxy Map gave 504 and the Warm galaxy-map step of the cache warm-up stopped with QUERY_TIMEOUT). Each sector's count is now read from `sector_system_counts`, the stored table the Sectors list already reads, which one background thread refreshes through the stored-count cache; a sector made since the last refresh has `system_count` null (not counted yet) and the map draws it as zero.
