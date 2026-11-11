**Headline:** the Galaxy Map's sector list no longer times out on a big galaxy: each sector's system count is a stored count made in the background.

### Changed

- MAP.424: `GET /api/galaxy/sectors` counted every system in the galaxy on every request, which passed the 10 second web limit once the galaxy held millions of systems (the Galaxy Map gave 504 and the Warm galaxy-map step of the cache warm-up stopped with QUERY_TIMEOUT). The counts now come from the same stored-count cache the table pages use (`countcache`): the answer made under the galaxy's current content stamp is served as it is, an older one is served while one background thread makes the next, and with no answer at all the sectors come back at once with `system_count` null (not counted yet) until the first count is made.
