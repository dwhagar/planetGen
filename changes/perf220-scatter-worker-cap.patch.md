**Headline:** The bright-star and massive-star scatter stages no longer start more workers than the machine's memory holds, so a big run is not killed for running out of memory.

### Changed

- The scatter stages (massive stars, bright stars, the neighborhood backfill and the dimmer band) now use at most 60 percent of the host's memory divided by 600 MB per worker (15 workers on a 16 GB host), whatever the run's worker count, and ask the workers an earlier stage left to finish their unit and exit when there are more. On 2026-10-31 36 workers took 16,285 MB of a 15,992 MB host in the bright-star stage and the run was killed. Other stages are unchanged. What a stage generates does not depend on the number of workers (each unit has its own seed), so the stars are the same; only the time can be longer on a big machine. generation-output: unchanged
