**Headline:** A Redis that stops answering now shows up in the log as one warning, not as a silently slower site.

### Changed
- **A silent Redis fallback now says so (OPS.165).** The page cache, the single-flight claims, the estimate cache, the galaxy status, the count cache, the object index and the queue's page-job and warm-up helpers used to catch a Redis error and carry on with no line at all. Each now calls `log.degraded`, which writes a debug line every time and a warning (so it reaches error.log) at most once in five minutes for the same place. generation-output: unchanged
