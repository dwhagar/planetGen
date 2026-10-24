**Headline:** The cache warm-up no longer fails on a galaxy with more chunks than one stage can queue as units.

### Fixed
- **A warm-up stage failed with "a stage queues at most 4000 units, not 4051".** The warm stages cut their keys into
  fixed chunks, so a galaxy with more than 4000 chunks was refused by the queue. A stage now makes its units a whole
  multiple of the chunk when the chunks would pass the cap (`units.coarsen`), and a unit builds its keys 64 nodes per
  request, stopping part way when the cache fills. The habitability and map-luminosity fills cut their id ranges the
  same way.
