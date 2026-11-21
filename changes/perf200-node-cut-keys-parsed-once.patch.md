**Headline:** The Galaxy Map's node cut parses a node key once and keeps the keys of a block's children, so the cut's own functions take about a fifth less time.

### Changed

- `parseNodeKey` keeps the block it made for each node key (the page parsed the same few thousand keys on every cut), `settledExpanded` keeps the keys of each block's children instead of parsing and formatting them again for every expanded block of every cut, and `blockOfKey` reads the number before the first dot without making an array (PERF.200 slice 2). The keys themselves are unchanged, so the stored nodes, the URLs and the order of the cut are the same.
- On the MAP.249 benchmark, the node cut and node memory functions went from 3,163 ms to 2,475 ms in the profile (22 percent less); the page's stars written (3,719) and stars (3,812) are the same, the requests 8 against 9 and the bytes 315,343 against 333,029 follow how the fetches fell in the run.

generation-output: unchanged
