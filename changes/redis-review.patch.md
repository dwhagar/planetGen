### Added

- Research: how Redis is used for generation, and how to use it better.

  Added a design note that reviews how planetGen uses Redis today, sets that against Redis and RQ best practice, and recommends how to rebuild the generation task process around Redis (one run record, each stage's units as queue items, progress in Redis, a properly configured Redis server). No code changes.
