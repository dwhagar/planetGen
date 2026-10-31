**Headline:** `update.sh` builds the site's caches before it reloads the site and checks `/api/health`.

### Changed
- The cache warm-up runs inside the update (it no longer starts in the background after it), so the first visitor finds the pages built and the health probe does not meet a site that is still warming.
- The health probe waits 10 seconds after the reload (`HEALTH_SETTLE_SECONDS`).
