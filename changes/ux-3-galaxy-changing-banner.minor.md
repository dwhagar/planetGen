**Headline:** Every page now warns visitors while a background run is changing the galaxy, with the time left to the next whole hour.

### Added
- **A banner on every page while the galaxy is being changed (UX.3).** Any visitor, signed in or not, sees "The galaxy is being changed by a background job" with an ETA in whole hours, rounded up and held steady (it falls as the estimate falls and rises only when even the low end of the estimate is longer). It reads the Redis run record of PERF.96, so a run started from the command line shows too; before a rate is known it says "a guess from earlier runs" or "working out how long it will take". A run with no sign of life, or on another database, shows nothing. Pages draw it when they are made and `/galaxy-status` refreshes it every 30 seconds.
