**Headline:** The Generate page answers a repeated size-and-time estimate from a cache, and costly cached answers are refreshed just before they expire.

### Added
- **The Generate page's estimate is cached (PERF.109).** The estimate (`--estimate-only`, up to two minutes for a large shell) is kept in Redis for five minutes under the job kind, its command line and the database, and served only while the generator epoch, the galaxy's data version and the recorded speed statistics are unchanged. Two admins asking for the same estimate at once run it once.
- **XFetch for costly entries that expire by age (PERF.110).** An answer in the page cache that took half a second or more to fetch (the Search page's facets, for one) is fetched again by one reader shortly before it expires, so a burst of visitors does not all miss in the same second. The admin Stats page keeps its database figures for 30 seconds when they take half a second or more to count, and one request makes them again just before they go stale; a galaxy that counts quickly always shows live figures.

### Changed
- `singleflight.claim` and `singleflight.wait` take the length of the claim and of the wait, for builds that take longer than a tile.
