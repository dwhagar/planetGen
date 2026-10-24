**Headline:** Saving the sector paths of a large run is spread over the workers.

### Changed
- **The settle stage is cut into units of 100 sectors (PERF.115).** A run saved the paths of every sector it made, and of the sectors around them, one after the other in a single process (182 ms a sector at the default galaxy size). The sectors are found once and saved by the workers in chunks, each path as before: the same rows are stored, and on 100 sectors the stage took 11 s on three workers instead of 18 s.
