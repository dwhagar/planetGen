**Headline:** The site's caches are rebuilt by themselves after a generation run, an orbital update or an edit, without waiting for the next update.

### Changed
- **The warm-up starts after every change (PERF.141).** Anything that moves a data version tag (a generation run, `planetgen orbits`, an edit through the API, a change found in the galaxy) now asks for `planetgen warm`, which starts in the background as the web server's user, waits `cache.warm_settle_seconds` (default 60) so a burst of edits is one warm-up, and rebuilds only what the changes touched. A change that arrives during a run gets another round; nothing starts during an update (it warms at its end) or with `cache.warm_on_change` off (`PLANETGEN_WARM_ON_CHANGE=0`). Two new settings, `cache.warm_on_change` and `cache.warm_settle_seconds`.

TODO items resolved / Added: PERF.141 (slice 2: the automatic start; the admin page follows) / Added: none
