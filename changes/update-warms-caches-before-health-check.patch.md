**Headline:** `update.sh` starts building the site's caches early and reloads the site only once the opening views are built.

### Changed
- The cache warm-up starts in the background as soon as the code, database, permissions and Redis are ready, and the update waits (at most 20 minutes) for its opening phase before it reloads Apache. The other phases go on after the update ends.
- The health probe waits 10 seconds after the reload (`HEALTH_SETTLE_SECONDS`).
