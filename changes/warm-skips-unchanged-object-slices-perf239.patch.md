**Headline:** A warm-up with nothing changed no longer visits every object page again: its object phase takes seconds, not minutes.

### Changed

- **The warm-up skips object slices that are still stored (PERF.239).** A run with nothing changed built nothing but still took 4 minutes 15 seconds on a 14-sector test galaxy, 214 seconds of it the objects phase, because the "already stored" check was two requests through the Flask app for every system, sector and phenomenon. While a slice of an object family is warmed, the page cache now records the Redis key of every answer its pages read or stored, the counter of each tag they read and the incarnation; the next run reads the counters and checks the keys in two round trips and does not visit a slice whose counters, incarnation and answers all stand. A tag that moved (any change that touches the slice), an answer that left Redis, a page that failed when the slice was first warmed, or a Redis that does not answer visits the slice as before, so the rule for when an answer is stale is the cache's own. The objects phase of an unchanged run on the test galaxy went from 214 s to 4.5 s (10,014 objects, one process). The slice record lives half as long as its answers (`page_cache.shared_seconds`), so an answer nearing its expiry is looked at again in time.

generation-output: unchanged
