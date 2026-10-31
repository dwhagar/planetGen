**Headline:** While the caches are rebuilt, requests that would build one are answered 503 instead of killing the API daemon.

### Changed
- A request that needs a cache that takes minutes to rebuild (the Galaxy Map's nodes, the shape, the lists, the search, and the galaxy's cluster field wherever a route reaches it) is now refused with a 503, `Retry-After` and a plain message (HTML page or JSON) while an update, a generation run, a hand-started cache warm-up or maintenance mode is on (`admin/cachebuild.py`, `web/rebuild_gate.py`). Cheap requests and anything the cache already holds still work. Before, such a request built the cache in the web daemon, ran past mod_wsgi's 60-second `request-timeout`, and the daemon was killed with every other request in flight ("Truncated or oversized response headers"; the shape took 56 to 65 s to rebuild on the test server because a running generation keeps changing the rows it is made from).
- Routes that read only the stored grid numbers (the sector edge, the shell limits, the disk scale length) no longer build the cluster field: `store.get_galaxy_shape(conn, light=True)`.
- A stale marker is ignored: a dead process, or one older than 12 hours.
