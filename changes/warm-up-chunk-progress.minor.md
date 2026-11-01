**Headline:** The warm-up now says which chunk of how many it is on, and shows a timed progress bar on a terminal.

### Changed
- `planetgen warm --progress` (update.sh stages 1 and 2): a stage cut into chunks, such as `galaxy-deep-nodes`, shows `chunk N of M` with the total known up front, and the overall bar advances with the chunks instead of standing still until the stage ends.
- With a terminal the rich bars show elapsed time, time left and seconds per chunk. In a log (`/tmp/planetgen-warm-map.log`, cron) each chunk leaves one line with its number, the nodes it held (the 64 is per chunk), the built and already stored counts so far, elapsed and time left, one line per chunk up to 20 chunks and every few chunks beyond that, then one summary line for the stage.
