# Cache warm-up and invalidation (PERF.104)

Informs: PERF.104, PERF.2 (the page cache), PERF.64 (stored counts), MAP.134 (the Galaxy Map warm-up), PERF.86 (Redis rules), OPS.51 (maintenance)
Status: built (Core 2, 2026-10-16). Nothing here changes a table; no migration.

## The problem

After `update.sh` the first visit to a page builds it from the database, and today the warm-up covers one thing, the Galaxy Map's opening tiles. The page cache (`web/lib/pagecache.py`) is a dict inside each WSGI process, so no outside script can fill it, and it drops a database's whole cache whenever the galaxy's content stamp moves. A warm-up therefore needs a store every process shares, and "changing one sector rebuilds only the caches that read it" needs the cache to know what each answer read.

## Design

**One shared store, Redis (PERF.86 rules).** The page cache keeps its in-process copy for speed and adds Redis behind it. Redis is the only place a warm-up can write for every worker to read. It runs with AOF `everysec` and `noeviction`, so every key carries a TTL and every value is capped (one body at most a quarter of `page_cache.max_mb`). A cold Redis is an empty cache, never wrong: a read that cannot reach Redis is a miss and the answer is not stored.

**Every entry names what it read (tags), and every change bumps tags.** Redis holds one hash per database, `planetgen:ver:<db>`, of tag name to counter, and a field `_id` (an incarnation token made when the hash is first written, so a flushed Redis cannot revive an old in-process copy). An entry stores the counters of its tags at the time it was fetched; it is served only when every counter still matches. A stale copy is therefore never served as current, in this process or another.

| Tag | Bumped when |
|---|---|
| `epoch` | the plan, the galaxy shape, the naming key, a release, a bulk delete, or a write the table below cannot place. Every entry carries it. |
| `sector:<id>` | that sector is edited, regenerated, has systems added or removed, or a system in it is edited (ids as printed, API.23) |
| `system:<id>` | that system is edited or deleted |
| `phenomenon:<type>:<id>` | that phenomenon is edited or deleted |
| `sector-set`, `system-set`, `phenomena-set` | an object of that kind is added or deleted (who is next to whom, counts) |
| `sector-data`, `system-data`, `phenomena-data` | any object of that kind is added, deleted or edited (what a list shows) |
| `system-any` | a star, planet, moon or facility is edited and the route does not name its system |
| `population` | the population pass ran, or a species, polity or territory was edited |
| `positions` | the orbital update moved objects |

Every bump also adds one to `_n`, which `begin`/`put` compare so an answer fetched while data changed is not stored.

**What an answer reads** (`web/lib/pagetags.py`, the one place that says; the ids come from the answer's body):

| Answer | Tags |
|---|---|
| Sector page | `epoch`, `sector:<its id>`, `sector:<each neighbour's id>`, `sector-set`, `phenomena-data`, `population`, `positions` |
| System page, text, sections, scene | `epoch`, `system:<id>`, `sector:<its sector>`, `system-any`, `phenomena-data`, `population`, `positions` (the two set tags when the body names no sector) |
| Phenomenon page | `epoch`, `phenomenon:<type>:<id>`, `sector-data`, `phenomena-set` |
| Sector list | `epoch`, `sector-data`, `system-set`, `phenomena-set`, `population`, `positions` |
| System lists | `epoch`, `sector-data`, `system-data`, `system-any`, `population`, `positions` |
| Phenomena list | `epoch`, `sector-data`, `phenomena-data` |
| Polities, species | `epoch`, `population`, `system-data` |
| Anything not listed (search, near, nav, galaxy) | `epoch` plus every tag above except the per-object ones |

So renaming one sector bumps `sector:<id>` and `sector-data`: its page, its neighbours' pages and its systems' pages are rebuilt, and the lists; no other sector's page and no other system's. Renaming one system bumps `system:<id>`, `system-data` and `sector:<its sector>`: its page, its sector's page and its siblings' pages are rebuilt, no other sector.

**Who bumps.** (1) A write through the API (`web/__init__.py`) maps the route to a change kind (`pagetags.change_for_write`) and bumps its tags before the response goes out; for a system the sector is looked up before the write runs; a route it cannot place bumps `epoch`. (2) The command line records its change when it finishes (`cli/generate.py`: `plan` is `everything`, `population` is `population`, every other writing command `generated`; `cli/orbits.py`: `orbits`). (3) The backstop for everything else, as before: the content stamp (`GET /api/galaxy/changes`, checked every `page_cache.stamp_seconds`). The first process to see a new state claims it in Redis (`dataversion.claim_state`) and bumps `sector:<id>` for each sector named (`galaxy_changes` now lists them, placed or not) and `sector-data`, plus every set and data tag when something was added; a change it cannot place (`full`) bumps `epoch`. The state the last claim recorded is what the next check asks about, so a change made while no web process ran is found by the first check after. While a fill makes the galaxy `busy` (PERF.38) nothing is bumped from the second busy check on and entries age out, the one documented case where an answer can be up to `max_age_seconds` old. The age limit stays the last backstop.

**Counts.** `db/countcache.py` already keys its answers by the stamp and refreshes them in the background; the warm-up asks for each table page's counts once, which stores them.

## The warm-up

Superseded in part (Boss, 2026-10-18, 5:43 PM EDT, PERF.139): the flat family order below is to be replaced by budgeted warm-up phases in place of the level cap.

`planetgen warm` (PERF.127; it was `python -m planetgen.cli.warm_map`), run after `update.sh` as the web server's user. It covers these families, each its own command on the work queue (`warm-galaxy-map`, `warm-galaxy-tiles`, `warm-lists`, `warm-sectors`, `warm-systems`, `warm-phenomena`, `warm-search`; `queue/commands.py`) that warms that family and nothing else and can run alone (`planetgen command-run`). `planetgen warm` builds the list of them and the wrapper runs them in order as a maintenance run (`runrecord.MAINTENANCE`: its own place in Redis, so it does not block a generation run and the visitors' banner never shows it):

| Family | What it warms |
|---|---|
| `galaxy-map` | the opening tiles and stage (existing `warm_opening_view`), the galaxy shape, sector list, phenomena, superbubbles, shell specs, bright-star status |
| `galaxy-tiles` | the Galaxy Map's tiles (PERF.130): every tile of the first `tile_cache.warm_levels` levels (default 8, 4,713 tiles on the scaled probe galaxy), then the tiles of a zoom run from the opening view in to the centre, down to level `tile_cache.warm_zoom_levels` (default 14), 7,112 tiles in all; in units of 128 tiles on the queue, coarsest level first |
| `lists` | the first page of the sector, system, uncharted-system, phenomena, polity and species lists, and each list's counts |
| `sectors` | every sector page (the sector detail), nearest the galaxy centre first, up to `page_cache.warm_max_objects` |
| `systems` | every system page, same order and cap |
| `phenomena` | every phenomenon page, same cap |
| `search` | the search page's empty and default answers and its index counts |

**`galaxy-tiles` (PERF.130).** The only family cut into many units. Its expansion lists the tiles (`warmup.tiles_to_warm`) and cuts them into chunks of `TILE_CHUNK` (128, a tile request's most); a unit builds one chunk through `tilecache.fetch_tiles`, the call the map's own request makes, so what is stored is what a visitor is served and a tile that is already stored and current is not built again (a second run builds nothing). Two parts, in this order: (1) every tile of levels 0 to `warm_levels` meeting the galaxy's cut box, coarsest level first and nearest the origin first within a level, so a stopped run has built what the most views use; (2) the zoom run: the union of the cuts (`viewport.fetch_cut_tiles`) for the camera on the opening camera's line through the origin (the page's starting angles round the central black hole), from the opening distance in to the finest level, one fade band (1.15) closer at a time, keeping tiles down to `warm_zoom_levels`. The count is the galaxy's own: each level has eight times the tiles of the one above inside the box, so `warm_levels` 9 is 37,481 tiles (about an hour of one worker) and a run past 150,000 tiles is refused with its count. Tiles deeper than the run are built on first request (PERF.132 builds the ones around a visit ahead of time). Measured on the scaled probe galaxy (no generated sectors; Boss's galaxy will cost more at depth): a tile builds in 0.05 s (finest) to 0.5 s (coarsest) and is 0.1 to 21 KB packed, so the defaults (7,112 tiles) took 12 minutes 17 seconds in one process on that galaxy, 40 MB on disk (under `tile_cache.max_mb`, 200; to be replaced by `cache.max_mb`, default 1024, a settings-page slider and a config key, Boss, 2026-10-18, 5:43 PM EDT, PERF.139), and a second run took 1.0 second and built none.

*Settled-view time* (the server's time to hand a fresh visitor every tile of one view, asked in request-sized chunks one after the other, for the camera on the zoom run's line at each distance, on that galaxy): without warming 5.1 / 9.0 / 18.6 / 53.5 / 12.5 / 8.4 / 3.0 / 1.2 s at 100,000 / 30,000 / 10,000 / 3,000 / 1,000 / 300 / 100 / 30 pc, 111 s over the eight views; after the warm-up 0.0 / 0.0 / 1.7 / 1.2 / 0.1 / 0.1 / 0.0 / 0.1 s, 3.2 s over the eight. What remains is tiles the run's one-step grid does not hold (a distance between two grid steps asks for a few tiles of a neighbouring cut); PERF.132 builds those around a visit. The tile cache is keyed by the tile stamp (PERF.133, below), which leaves the release number out, so a release that changes neither the tile format nor the generator epoch keeps every stored tile and the `planetgen warm` that `update.sh` starts finds them stored ("0 built").

**Building around a visit (PERF.132).** `web/lib/tileprefetch.py`. When a `/galaxy/tiles` response has gone, the route hands the request's tile keys, and the camera and target the page sent (`?at=cx,cy,cz,tx,ty,tz`, parsecs, six plain numbers or ignored), to `tileprefetch.schedule`, which queues what the visitor is likely to ask for next and does not have: (1) the tiles of the cut for the camera one, two, three and four fade bands closer along its line to the target (`ZOOM_STEPS`; the same `fetch_cut_tiles` the map asks with), nearest step first, then (2) the neighbours of the requested tiles at their own level and the tiles of the next level inside them, for a visitor who pans. They go to the bulk queue (never the interactive one, so an admin's work is not behind them) as jobs of 16 tiles that a worker builds into the same disk cache through `tilecache.fetch_tiles`; one worker is started if none serves the queue. Rules, each a `tile_cache` setting or a fixed constant: at most `prefetch_per_visit` (128) tiles per request, `prefetch_per_minute` (600) in all, nothing while `prefetch_max_queued` (8) jobs wait, nothing for ten minutes after `prune` deleted files (the cache is full), nothing during an update, nothing for a tile stored or already claimed (a Redis key per tile and generation, 15 minutes), so a visit to a built area queues nothing and a second visitor to the same area causes no build. It fails open: no Redis, no cache directory, a bad `at` just means no background build. It does not queue stages or windows: those are not tiles of the octree and the map asks for them one block at a time.

**A tile is built again only when what it holds changed (PERF.133, Boss 2026-10-18 6:05 PM EDT).** The stored tiles are keyed by `query.galaxy_tile_stamp`, a hash of the galaxy shape, the bright-star seed, the naming key, `query.TILE_FORMAT`, the generator epoch and the change counters; the release number is not in it (it stays in `galaxy_content_stamp`, which the pages and counts follow). So (1) nothing changed: a second warm-up, or the `planetgen warm` after a release, builds no tile and says `N built, M already stored`; (2) a sector edited, added or given a system: `GET /api/galaxy/changes` names its tile at every level, only those files are deleted, and the ones that were stored are queued at once for a rebuild on the bulk queue (`tilecache.take_rebuilds`, `tileprefetch.rebuild`, coarsest first, at most 2,048 a change, off with the prefetch, during an update and while the cache is full), so the next visit finds them built; (3) a deleted sector, a new shape, tile format or epoch, more than 1,000 changed sectors, or a backfill: a new generation, rebuilt by the next warm-up. A change to what a tile holds or how it is packed must raise `query.TILE_FORMAT`; `tests/test_tile_format.py` pins the bytes of one fixed tile and fails until the number and its digest are updated. The stored tiles are never dropped for age; what the budget does when the cache is full is PERF.139 and PERF.133's second part.

Each family ends with one line, `<family>: <n> warmed (<m> already cached, <k> left out by the cap) in <s> s`, printed by the command, written to the activity log (the admin's log) and returned as the unit's result. `planetgen warm` exits non-zero if any command fails, and says which (the run stops at the failing one, as every wrapper list does); it has no silent fallback (a Redis that does not answer, another warm-up running, or an update running are errors, printed). A cap that leaves objects out says so; those fill on first visit. `planetgen warm lists search` runs just those; there is no `--no-rq` and no `--only`, since the order and the choice live in the wrapper's list and Redis is required. A family warms by visiting the pages through the site's own Flask app (`warmup.build_app`: open during the sector-switch maintenance, not rate limited), so what it caches is exactly what the page asks for; a page that does not answer 200 fails the family. The activity log has `GEN warm.family` or `GEN warm.failed`. `planetgen warm` waits for a running update to finish before it queues (the queue takes no job during one, OPS.53), so `update.sh` can start it in the background.

It runs inside the maintenance flag the update set (OPS.51). A list of every family ends with the command `end-map-rebuild-maintenance`, which ends the sector-switch maintenance after every family has succeeded: it is the list, not a family, that knows the maps are then rebuilt. It registers with the update guard (OPS.53) so a later update lists and cancels it.

## Done when

After `update.sh`, the first visit to the Galaxy Map, a sector, a system, a phenomena or system list and a search answers from the cache. Changing one sector rebuilds only the caches that read it. A test per kind of change proves a stale copy is never served as current.

## Not in this slice

Warming past the cap (a galaxy of millions of sectors warms its nearest `warm_max_objects`, the rest fill on first visit); resolving a star, planet or moon edit to its system (those routes bump every system page, `system-any`, instead); a per-page progress bar on the runs page (the node shows its state and the one result line).

## Without Redis

The shared tier needs Redis (PERF.86 makes it required for generation anyway). With none, the page cache stores nothing (a copy it could not validate is never served), `planetgen warm` refuses with a printed error (Redis is required for every run, PERF.127 removed the inline `--no-rq` mode).
