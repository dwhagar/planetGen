**Headline:** The table sweep now tests belts and asteroid fields at a physical disk temperature, and the two taxon E rows that no belt reaches are marked unreachable on purpose.

### Changed
- **The table sweep draws a disk temperature (TEST.545).** The sweep had no range for `temperature_k`, so the taxa table of a belt or asteroid field was tried at 0 K and raised "the temperature must be positive". It draws from 10 to 2,500 K now. The taxon E rows are never reached because the taxon kernel gives that taxon no weight at any temperature (its source is not found yet); the sweep notes say so.

generation-output: unchanged

TODO items resolved: TEST.545.
