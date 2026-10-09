# Galaxy Map wire format: what the browser downloads and what could be smaller (MAP.147)

Research only, 2026-10-09 (Research Lane 3). No code, TODO.md or PR changes. Boss decides.
Scripts and raw numbers: `research/map-wire-format/bench/`.

## 1. Answer in brief

**Section 5.6's 52 bytes a star is the GPU buffer, not the download.** On the wire each star is a JSON object of about **270 bytes** (45 bytes after gzip), five times what 5.6 assumes. At the 70,000-star cap that is about 19 MB of JSON, 3.2 MB gzipped, not the 3.6 MB of 5.6 or the "about 8 MB" in the `GALAXY_VIEW_MAX_STARS` docstring. Quantising the GPU attributes (52 to 20 bytes) would not touch any of that; shrinking the tile format will.

Recommendation, in two steps that need no new design decision to start:

1. **Now (small, no client contract change): trim the tile JSON on the server.** Drop what the client never reads and round the rest. Measured on 400 real tiles: **gzip bytes fall 48%** (3.52 MB to 1.83 MB), raw bytes 52%, JS keeps working unchanged. Also: stop prefetching 18 times more than the view needs, serve prebuilt bytes instead of re-parsing and re-serialising JSON on every hit, and precompress with brotli at cache-write time (another 10 to 29%).
2. **With MAP.146 and the nested tile lists (zoom-star-visibility item 2), one cache-stamp bump: a packed binary tile** (quantised planes, 12.9 to 16.7 bytes a star). Measured **gzip bytes fall 78%** against today (3.52 MB to 0.79 MB; 4.3x), raw bytes 96%. Store tiles in IndexedDB, not localStorage (section 4.4 shows localStorage is full after one view).

Answer to the open question ("build the quantised format straight from 5.6?"): **no**. 5.6 aims at the wrong buffer. Step 1 first (days, reversible), step 2 once MAP.146's keys and the nested lists are decided, because any format change bumps the cache stamp and the three should land together.

Ranking (best value first):

| # | Option | Wire change (gzip, 400 real tiles) | Client change | Cost | Verdict |
|---|---|---|---|---|---|
| 1 | Server trim: drop `placed`/`planned`/`filled` and unused star fields, letter-only `star_type`, rounding (V1) | -48% (3.52 to 1.83 MB) | none | about 1 day | Do now |
| 2 | Fix prefetch (zoom-in only, after idle, not on save-data), serve prebuilt bytes, brotli precompressed | cuts the 533 KB opening-view prefetch; warm hit 22-75 ms to about 2 ms; brotli -10% (binary) to -29% (JSON) | prefetch policy only | about 1-2 days | Do now |
| 3 | Packed binary tile, quantised planes (V3/V4) | **-76% / -78%** (3.52 to 0.83 / 0.79 MB) | decoder, storage, tests | 3-5 days | Do with MAP.146 + nested lists |
| 4 | Nested lists: a finer tile omits stars a coarser tile already sent (Lane 1's item 2) | 31-43% of bright records repeat at levels 3-6, 100% at 7-12 | key by rank, carry parents | in Lane 1's item | Fold in with 3 |
| 5 | IndexedDB (or per-tile immutable HTTP caching) instead of localStorage | no change to first load; revisit hit rate up about 10x | storage layer | 1-2 days | With 3 |
| 6 | Columnar JSON (V2) | -60% | rebuild objects in `tileStars` | about 2 days | Skip: gets 60% of the way to binary for most of the work |
| 7 | Quantise the GPU attributes (5.6 as written) | none | `setStars` | 1-2 days | Defer: saves GPU memory and re-upload, not download |
| 8 | Request batching | none: already 1-3 requests per action | none | none | Nothing to do |

## 2. What the browser downloads today

Measured with real Chromium (software GL) against the real Flask app on a local database: 412,924 bright stars (floor 10,000 L☉, see section 7), 222,739 phenomena, 17,176 generated stars in one cluster at 7,963 pc from the core, 1,369 sectors. gzip level 6 was applied to every response in the harness, as Apache's `mod_deflate` does in `examples/apache/planetgen.conf.example`. Throttling was set with the Chrome DevTools protocol: "DSL" = 10 Mbit/s, 40 ms; "slow 4G" = 1.6 Mbit/s, 150 ms.

### 2.1 The format

`GET /galaxy/tiles?tiles=level/ix/iy/iz,...` returns `{stamp, generation, tiles: {key: {placed, planned, filled, clouds, stars, generated, points}}}` as JSON with `Cache-Control: no-store`. The opening view's tile (level 0) is embedded in the `/galaxy` HTML. A star is `{id, x, y, z, luminosity_sol, temperature_k, radius_sol, star_type ("B4V Blue-White Main Sequence Star"), population, yerkes_class, ring_index, layer_index, ring_slot_index, system_id}`, plus `name` for generated stars. Floats arrive rounded to 4 decimals (`_round_floats` in `tilecache.py`), which is still 5.5 digits more than the screen can show.

### 2.2 Bytes by zoom (one view = the tiles the JS asks for at that camera radius)

| Camera orbit (pc) | Tile level | Tiles | Stars+points | Raw JSON | gzip | brotli-11 | Unused sections (placed/planned/filled) | Cold server | Warm server |
|---|---|---|---|---|---|---|---|---|---|
| 60,000 | 0 | 1 | 400 | 166 KB | 29 KB | 21 KB | 32% | 181 ms | 5 ms |
| 15,000 | 1 | 8 | 3,200 | 998 KB | 171 KB | 119 KB | 9% | 558 ms | 25 ms |
| 4,000 | 3 | 8 | 3,219 | 1,055 KB | 192 KB | 134 KB | 13% | 1,038 ms | 23 ms |
| 1,000 | 5 | 8 | 2,108 | 738 KB | 137 KB | 95 KB | 19% | 2,008 ms | 20 ms |
| 500 | 6 | 14 | 2,652 | 890 KB | 160 KB | 110 KB | 16% | 1,142 ms | 25 ms |
| 250 | 7 | 12 | 1,277 | 487 KB | 97 KB | 67 KB | 29% | 436 ms | 14 ms |
| 120 | 8 | 12 | 1,004 | 404 KB | 84 KB | 59 KB | 35% | 274 ms | 12 ms |
| 60 | 9 | 16 | 4,625 | 1,422 KB | 242 KB | 156 KB | 18% | 518 ms | 31 ms |
| 30 | 10 | 24 | 6,157 | 1,897 KB | 320 KB | 202 KB | 18% | 368 ms | 75 ms |
| 15 | 11 | 16 | 7,215 | 2,126 KB | 353 KB | 202 KB | 15% | 377 ms | 54 ms |
| 8 | 12 | 14 | 4,369 | 1,281 KB | 195 KB | 132 KB | 15% | 99 ms | 27 ms |
| 4 | 12 | 8 | 3,534 | 999 KB | 155 KB | 105 KB | 11% | 29 ms | 27 ms |

A view is 1 to 24 tiles, 0.17 to 2.1 MB raw, 29 to 353 KB gzipped. Per star (65,403 stars in 233 tiles): **270 B raw, 45 B gzip, 41 B brotli-5**. Server time: a cold tile is built from the database in 0.1 to 2.5 s; a warm (disk-cached) request takes 5 to 75 ms, almost all of it parsing the stored JSON and serialising it again (the cache keeps parsed dicts, not bytes).

Cap check: 70,000 stars x 270 B = **18.9 MB raw, 3.2 MB gzip** at today's format. The worst case is reachable on paper (27 tiles x (400 bright + 850 generated) + 8 detail tiles x (400 + 4,000) is about 69,000); my local data does not reach it (section 7).

### 2.3 Sessions (cold browser profile)

| Scenario | Tile requests | Tiles | Tile bytes (raw / gzip) | Static (requests / gzip) | First star: local / DSL / slow 4G |
|---|---|---|---|---|---|
| Open `/galaxy` | 3 | 28 | 3.05 MB / 533 KB | 173 / 536 KB | 1.29 s / 3.04 s / 7.86 s |
| Open a sector link (`?sector=`), 10,409 stars | 3 | 38 | 4.93 MB / 839 KB | 173 / 543 KB | 1.02 s / 3.02 s / 8.11 s |
| Reload, same profile | 0 | 0 | (localStorage hit) | 109 / 58 KB (dev-server revalidation) | 1.1 s |

- **The first star does not wait on a tile.** The opening tile is inside the 35 KB gzipped HTML (29 KB of it). What the first star waits for is the 173 static files (543 KB gzip, 1.8 MB decoded; JS modules and three.js). On slow 4G that is the 7.9 to 8.1 s. A smaller tile format will not move this number; bundling or preloading the modules would. Not part of MAP.147; a candidate item. (In production Apache serves `/static/` with `immutable` year-long caching per the example config, so a repeat visit makes no static requests; the 109 reloads above are the Flask dev server revalidating.)
- **The opening view downloads about 18 times more tile data than it shows.** Of the 3 requests, 1 is the visible tile (embedded); the other 2 are the prefetch of one zoom step in and out (`PREFETCH_ZOOM = 3`): 28 tiles, 533 KB gzip, 3 MB raw. On a phone that is data spent before the user does anything. It starts only after the real fetch ends, so it does not slow the first view, but it does compete with the user's next action.
- **Drill-down steps**: crossing a level costs 1 to 3 requests. Typical: 8 tiles 206 KB gzip (1.3 MB raw) for the new level, then 12 tiles 165 KB gzip for the prefetch. Measured duration of the 206 KB request: 338 ms local, 332 ms DSL, **1,196 ms slow 4G**. That matches bytes / throughput (206 KB / 200 KB/s = 1.0 s) so the model below holds.
- **Stars on screen change**: 10,409 to 10,999 to 6,245 and so on, each change a full rebuild of the star geometry: 3,200 stars = 166,400 bytes (52 B x 3,200) uploaded; 10,999 stars = 572 KB. Every tile arrival re-uploads everything (`renderFromCache` builds a signature of every star key and calls `setStars` on any change). Cheap on a real GPU; can't be measured here (software GL).
- **Stage (drill-down block) responses are small**: 0.25 to 7.4 KB (0.4 to 1.5 KB gzipped). Nothing to gain.

## 3. What the client uses

`galaxymap3d.js` reads only `tile.clouds`, `tile.stars`, `tile.generated` and `tile.points` (`tileStars`, `renderFromCache`). **`placed`, `planned` and `filled` are never read** by any file under `html/static` (grep of the whole directory; the block view uses `/galaxy/stage` instead). They are 9 to 35% of every tile (section 2.2, "Unused sections") and are written to the disk cache too. Of each star's fields the JS uses `x y z luminosity_sol temperature_k radius_sol star_type (first letter only, for the class filter) id system_id` (system_id only for `null` = "uncharted"); `ring_index`, `layer_index`, `ring_slot_index`, `population`, `yerkes_class` and a generated star's `name` are never read (stars cannot be picked, MAP.101). Treat this as true of main at 8.0.783; the build re-greps before removing.

## 4. The options, measured

All numbers are on the same 400 real cached tiles (65,403 stars), each variant built from the real records. Prototype encoders live in `bench/protolib.py`; nothing was added to the repo.

| Variant | What it is | Per star raw / gzip | 400 tiles raw / gzip / brotli-11 | At the 70,000 cap (gzip) |
|---|---|---|---|---|
| V0 | Today | 270 / 45 B | 20.2 MB / 3.52 MB / 2.49 MB | 3.2 MB |
| V0b | Today minus `placed`/`planned`/`filled` | n/a | 17.7 MB / 2.99 MB / 2.12 MB | n/a |
| V1 | Same dict shape; only the 9 used fields; `star_type` one letter; x/y/z to 3 decimals, luminosity 4 and radius 3 significant digits, temperature to 10 K; no `name` | 148 / 28 B | 9.7 MB / 1.83 MB / 1.30 MB | 1.9 MB |
| V2 | V1 values as columns (`{x:[...], y:[...]}`) | 60 / 21 B | 4.0 MB / 1.42 MB / 1.10 MB | 1.5 MB |
| V3 | Binary planes: x/y/z as uint16 inside the tile's cube; log10 luminosity as uint16 (0.0002 dex); log temperature and radius as uint8; class + `has system` flag as uint8; `id` uint32; `system_id` as uint32 | 16.7 / 12.3 B | 1.12 MB / 0.83 MB / 0.72 MB | 0.86 MB |
| V4 | V3 sorted along a Morton curve, key deltas and id deltas as varints | 12.9 / 11.7 B | 0.87 MB / 0.79 MB / 0.72 MB | 0.82 MB |

Per view (gzip, the table from section 2.2 re-encoded):

| Orbit (pc) | Stars | V0 | V0b | V1 | V2 | V3 | V4 | V4 / V0 |
|---|---|---|---|---|---|---|---|---|
| 60,000 | 400 | 29 KB | 20 KB | 12 KB | 10 KB | 5 KB | 5 KB | 16% |
| 15,000 | 3,200 | 170 KB | 154 KB | 93 KB | 72 KB | 39 KB | 35 KB | 21% |
| 4,000 | 3,219 | 191 KB | 155 KB | 95 KB | 73 KB | 40 KB | 38 KB | 20% |
| 1,000 | 2,108 | 136 KB | 99 KB | 62 KB | 48 KB | 27 KB | 26 KB | 19% |
| 500 | 2,652 | 159 KB | 122 KB | 76 KB | 59 KB | 34 KB | 33 KB | 21% |
| 250 | 1,277 | 96 KB | 60 KB | 36 KB | 28 KB | 17 KB | 16 KB | 16% |
| 120 | 1,004 | 84 KB | 47 KB | 28 KB | 22 KB | 14 KB | 12 KB | 15% |
| 60 | 4,625 | 242 KB | 186 KB | 111 KB | 83 KB | 54 KB | 51 KB | 21% |
| 30 | 6,157 | 320 KB | 251 KB | 151 KB | 114 KB | 74 KB | 69 KB | 22% |
| 15 | 7,215 | 352 KB | 293 KB | 175 KB | 132 KB | 87 KB | 82 KB | 23% |
| 8 | 4,369 | 195 KB | 167 KB | 100 KB | 75 KB | 49 KB | 47 KB | 24% |
| 4 | 3,534 | 155 KB | 135 KB | 81 KB | 61 KB | 40 KB | 38 KB | 24% |

### 4.1 Compression on top

Brotli-5 saves about 10% over gzip-6 on JSON (3.52 to 3.16 MB), brotli-11 29% (2.49 MB), zstd-19 22%. On the binary form the gap shrinks to 9% (0.79 to 0.72 MB) because there is little left to find. The production Apache example enables `mod_deflate` (gzip) only; `mod_brotli` is an extra module and the example does not load it. The tile cache already writes each tile once to disk, so the build can precompress at brotli-11 at write time and serve the `.br` file with zero request-time CPU. That is the cheap route to the 29%.

### 4.2 Cost of the binary form

- Server: encoding 65,403 stars in Python took 0.29 s as V3 (4.5 µs a star, the same as `json.dumps` at 5.2 µs); the Morton sort in pure Python takes 22 µs a star (1.46 s), which numpy removes. Either way it runs once per tile at cache-write time next to a 0.1 to 2.5 s database build.
- Client: `JSON.parse` of the 17.7 MB of today's records takes 86 ms (node 22; trimmed V1 57 ms); decoding the V3 planes into the typed arrays `setStars` needs takes 12 ms. Parsing was never the browser's cost (in-page `Response.json` took 5 to 47 ms per step); the rebuild of all star geometry on each tile arrival is. Not measured on a phone CPU.
- Precision: uint16 inside a tile cube is edge/65,535 (1 pc at level 0, 0.25 mpc at level 12), under 1/100 of a pixel at the zoom each level is used. Luminosity and radius colour the star from `logShare` over fixed ranges (luminosity -4 to 6 dex, radius -1 to 3 dex); a uint8 step is 0.4% of the range. The 16-bit luminosity keeps the dimmest-star filter exact to 0.0002 dex.
- Rank order: Lane 1's smooth-fade design needs each star's rank in the tile list. V3 keeps list order (rank implicit, 16.7 B). V4 re-sorts, so it would add a 2-byte rank plane: 14.9 B, still smaller than V3.
- A "lean" record without `id`/`system_id` is 11 B (6 position + 2 luminosity + temperature + radius + flags). Nested lists remove cross-level duplicates, so the client no longer needs `id` to de-duplicate; `uncharted` becomes one flag bit. **Do not put the new 80-bit object ID (10 B) on the tile wire**: that would add 60% to a V3 record for nothing, since stars are not pickable.

### 4.3 Reuse between zoom steps

Of the bright-star records in finer tiles, **31 to 43% at levels 3 to 6 and 100% at levels 7 to 12** were already sent in a coarser ancestor tile (measured on the cached tiles that have an ancestor cached; generated stars repeat 34 to 80%). That is the "nested lists" change Lane 1 proposes (zoom-star-visibility item 2), and it also helps this item: it shaves the repeated records from every fine-level step. It does not need the binary format, and it changes the stamp, so land it with the format.

### 4.4 Browser cache

The tile cache is `localStorage`, capped by Chromium at **5,242,880 characters** (measured: 524 x 10,000-character entries then `QuotaExceededError`). A one-view session already stores 22 tiles / 2.2 M characters (the prefetch included); a sector link fills it (**40 tiles, 5.08 M characters**, largest tile 265 K). When a write fails `storeTile` purges other generations, then every stored tile of the database. Visiting 7 sector links in one profile left 40, 18, 12, 44, 19, 48, 65 tiles stored (5.08 M, 1.41 M, 0.64 M, 4.07 M, 1.11 M, 4.53 M, 5.07 M characters): the cache is wiped and refilled again and again, so revisits mostly miss. Each deep link downloaded 82 to 839 KB gzip. Binary tiles in base64 in `localStorage` would fit about 10 times more tiles; IndexedDB has no practical cap and stores the buffer directly, with no synchronous `JSON.parse` on the main thread. The alternative is immutable per-tile URLs (`/galaxy/tiles/<generation>/<key>`, `Cache-Control: public, immutable`) so the browser's own disk cache and Apache serve them; that gives up the stamp history's per-tile invalidation unless the stamp enters the URL, so it belongs with PERF.38/PERF.41, not here. (`Cache-Control: no-store` today means the browser cache is never used.)

### 4.5 Modelled time on a phone-class link

Throughput x bytes + RTT, validated by the measured 1,196 ms above (model: 1.03 s + RTT):

| Step | Today (gzip) | V1 | V4 | Slow 4G today / V1 / V4 | DSL today / V4 |
|---|---|---|---|---|---|
| Level change in a sector (8 tiles) | 206 KB | about 107 KB | about 43 KB | 1.2 s / 0.7 s / 0.35 s | 0.33 s / 0.08 s |
| Opening view's prefetch (28 tiles) | 533 KB | about 280 KB | about 110 KB | 2.7 s / 1.4 s / 0.55 s | 0.4 s / 0.09 s |
| Worst view at the 70,000 cap | 3.2 MB | 1.9 MB | 0.82 MB | 16 s / 9.5 s / 4.1 s | 2.6 s / 0.66 s |

(V1/V4 columns use the per-view ratios of section 4; not measured end to end because the client does not read them yet.)

## 5. Section 5.6, the 13 floats

`setStars` builds 13 floats a star (position 3, colour 3, size, core, glow, bright, born, clipped, uncharted): 52 B, 3.6 MB at 70,000 stars. Quantising colour, scalars and flags to bytes (about 20 to 24 B) would cut GPU memory and the re-upload on each tile arrival by about 60%, but not the download, and nothing here measured it on a GPU. The bigger client cost is that every tile arrival rebuilds and re-uploads all stars; a client that appends the new tile's stars (and drops the old ones) would remove it and is independent of the wire. Both are optional and small next to the 270 to 13 bytes on the wire. Above about 1e6 points 5.6 recommends level of detail by tile; the tile budgets (400 bright, 150 to 850 generated, 4,000 detail) already do that, so that part is satisfied.

## 6. Related items

- **MAP.146 (click-centred drill-down).** Research Lane 2's study (`research/drilldown-region-sizes/report.md`) finds the **cube tiles are already free-floating** and need no key change; the block-keyed cache is the *stage* cache (small, 0.4 to 1.5 KB gzipped, nothing to gain). So MAP.146 does not force a tile-key change. The format change and the nested lists are the only things that bump the tile stamp, and they should share one bump.
- **Zoom-star-visibility (PR #878, Research Lane 1).** Rank fade needs rank per star and nested lists; both fit V3 (list order = rank) and section 4.3. Item 3 (point objects from level 8 instead of 10) raises tile size at levels 8 and 9 by the points section (about 4 KB a view at level 10 today); small.
- **PERF.38 / PERF.41** (single-flight tile builds, what stamps the cache): serving stored bytes and precompressed files belongs with them.
- **Object IDs.** A star needs only a tile-local key on the wire; the 80-bit ID is for the details page.
- **GEN.64 position IDs as names.** Phenomenon points carry the 19-hex name (`0D01F3FFDF704001560`), 22 B a point with JSON overhead; at most 4 KB a view. Leave as JSON (the binary tile keeps `clouds` and `points` as JSON).

## 7. What could not be measured

- **A phone, or any real GPU.** Chromium here draws with software GL (SwiftShader) on a 4-core container, so frame time, first-paint, long tasks and rebuild cost are not representative; throttling is emulated (DSL 10 Mbit/40 ms, slow 4G 1.6 Mbit/150 ms) and does not model packet loss, radio wake-up or TCP slow start.
- **Production-scale data.** The bright-star floor was 10,000 L☉ (412,924 rows) rather than 1,000 L☉ (26.9 M rows), and only 17,176 generated stars in one cluster exist. So 400-star bright caps saturate only at coarse levels here, no tile reaches the 4,000-star detail cap, and no view reached 70,000 stars (the most was 10,999). Per-star bytes are properties of the record, so the extrapolation to the cap holds; counts and cold-build times (database) do not.
- **Apache.** gzip was applied by a level-6 WSGI wrapper standing in for `mod_deflate`; brotli sizes are from the Python library, not `mod_brotli`; HTTP/2 not tested (the example config does not enable it). Flask's dev server was the origin.
- **The binary client.** The decoder was timed in node on real V3 blobs, not wired into the page; time-to-first-star and step times with it are modelled from bytes, not measured. IndexedDB speed was not measured.
- **Safari/Firefox**, mobile quota rules for `localStorage`/IndexedDB, and concurrent users.
- **GPU upload and rebuild cost** on a real device (only the bytes uploaded were counted: 52 B a star per rebuild).
- **MAP.146 and the nested lists** are not built, so their effect on tile counts is not measured; 4.3 uses today's tiles.

## 8. Method

`bench/measure_tiles.py` requests the real `/galaxy/tiles` route in-process for 15 camera radii (replicating `neededTiles`), recording raw/gzip/brotli/zstd, the section split and cold/warm time. `bench/measure_browser.py` drives Chromium with Playwright: CDP network log (encoded and decoded bytes, request counts), a `bufferData` hook, `Response.json` timing, long-task observer, `localStorage` size, and `performance` resource timing; scenarios `home`, `sector` (a `?sector=` deep link plus wheel zoom), `drill` (clicks down the block stages), `--net dsl|slow4g`. `bench/prototype.py`, `view_totals.py`, `star_only.py` size the encodings on the cached tiles; `dump_corpus.py` + `decode_bench.js` time encode and decode; `ls_quota.py`, `ls_thrash.py` measure the `localStorage` quota and refill; `overlap.py` counts repeats between levels. Results are the `*.json` files beside them.
