**Headline:** The Galaxy Map spends less time on the node cut's bookkeeping, so the map's busy time on a dive, flight and climb falls by about a fifth (PERF.200).

### Changed

- A block's key is built by joining numbers into one string instead of an array and a join, the cut reads each block's key from its kept geometry (kept in one table per sector edge, no `key@edge` string per look-up) instead of building it again, and the cone test no longer makes an array per block.
- The page's node memory re-inserts a hit only when the node was last used over 500 uses ago, and gives a node used since a second chance at the front before it is evicted, instead of moving every hit; a draw looks each node of the cut up once instead of six times.
- Measured with the MAP.249 benchmark on this box: the named node-cut and memory functions took 4,229 ms of the profile before and 2,583 ms after (39 percent less), the page's whole busy time 9,720 ms before and 7,569 ms after (22 percent less); node requests (104) and bytes (2.81 MB) are the same.
