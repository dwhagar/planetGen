# Galaxy Map performance (MAP.249)

The Galaxy Map's node path is held to a budget by `tests/test_web_browser_map_perf.py`. This page says what it
measures, the numbers it found, and where the time goes.

## The benchmark

`large_site` (`tests/map_site_support.py`) serves a synthetic galaxy of 46,000 stars: an exponential disc (scale length
2,800 pc, 350 pc thick), 40 clusters of 300 stars and a core cluster of 4,000 stars (scale 25 pc) at the dive's target,
beside the fixture's own few hundred. Every node lists its brightest stars up to the real budget (`NODE_STARS` 96, a
sector `NODE_SECTOR_STARS` 300) and carries the glow of all the stars in its block, so each level of the dive loads what a
real galaxy's would, with the table scan replaced by a block index so the fixture server is no part of the figures.

The test opens the Galaxy Map, centres on the core cluster and then, step by step:

1. zooms in with the mouse wheel (400 delta a step) from 18,600 pc down to the closest view, 2 pc, through the 81, 27, 9,
   3 and 1 levels;
2. at the 3-block level, 60 pc from the target, jumps the centre six times round a circle of 150 pc and back;
3. climbs out with Escape, one level a press, to the whole galaxy.

After each step it waits for the view to settle (no flight, no node request in the air, no change to the star buffer for
400 ms) and records what the step cost: the time to settle, the node requests and bytes, the stars written and moved in the
buffer and the syncs that did it, the frames' median and longest interval and the longest task on the page's thread.
Setting `PLANETGEN_MAP_PERF_REPORT=<file>` writes every step's figures as JSON; `PLANETGEN_MAP_PERF_PROFILE=<file>` runs the
dive under the browser's sampling profiler and writes the 30 functions with the most self time.

## The limits

A step fails the test when it passes any of these. Counts and bytes are held close (the dive measured less than a third
of each); milliseconds are under software WebGL (SwiftShader) on a lane machine, so they are held at five times what was
measured or more, to catch a stall and not to rank two builds.

| Figure, per step | Measured at most | Limit |
| --- | --- | --- |
| Node requests | 8 | 12 |
| Node bytes | 261 kB | 1 MB |
| Time to settle (including the 400 ms quiet) | 3.2 s | 15 s |
| Longest task on the page's thread | 154 ms | 1.5 s |
| Longest frame interval | 217 ms | 2 s |
| Stars written or moved, against the buffer after the step | 1.3 times | 2 times |

The test also fails when the dive never passes through the 81, 27, 9, 3 and 1 levels. A new Boss sighting of lag is filed as
its own bug, not folded in here.

## What the run found (40 steps, 46,000-star galaxy)

| Phase | Steps | Slowest settle | Most requests | Bytes in all |
| --- | --- | --- | --- | --- |
| Open and centre | 1 | 3.0 s | 8 | 261 kB |
| Dive in | 26 | 3.2 s | 8 | 1.8 MB |
| Flight across a level | 7 | 1.2 s | 3 | 278 kB |
| Climb out | 7 | 3.1 s | 8 | 515 kB |

- **A step that crosses a coarse level takes 3 s to settle**, and most of that is the nodes arriving one after another: up
  to eight syncs of the star buffer in one step, each as another node lands. Below the 27 level a step settles in about
  a second, and in 0.4 s (the quiet wait alone) from the 3 level down.
- **Stars written per step stay inside the buffer.** The climb is the heaviest (1.3 times the buffer), as each level
  it climbs to replaces nearly every star on screen; the dive and the flight rewrite only the chunks that came or went.
- **The page's thread is mostly idle**: the longest task was 154 ms and the longest frame 217 ms (the frames are
  software-rendered, with a median of 33 to 67 ms).

## Where the main thread's time goes

Under the sampling profiler (the busy time of the dive, flight and climb, the idle time left out), about half the page's
own time was the node path's bookkeeping, not drawing:

| Self time (ms) | Function |
| --- | --- |
| 733 | `consider` in `nodeCut` (`galaxynodecut.js`), the cut's per-block test |
| 706 | `formatDrillKey` (`galaxyprisms.js`), the block's key built with `Array.join` for every block of every cut |
| 693 | the garbage collector |
| 498 and 330 | `Map.get` and `Map.put`, the node cache (`memoryGet` re-inserts on every hit: 386) |
| 415 | the cut's sort comparator (`galaxynodecut.js:254`) |
| 294 | `geometry` (`galaxynodecut.js`), a block's box worked out again at every cut |
| 178 | `parseNodeKey` |

That is 4.1 s of 7.9 s in the listed functions. Causes, in the order they would pay back: the cut builds the key string
of every block it considers (and parses it again to settle it) where the block could carry it; `geometry` is worked out
again for the same blocks at each cut; `memoryGet` moves an entry on each hit even when it is already the newest. None of
this is a fault in what the user sees (the longest task is 154 ms), so the fix is a separate item once a sighting needs it.
