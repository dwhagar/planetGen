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

`planetgen warm` (PERF.127; it was `python -m planetgen.cli.warm_map`), run after `update.sh` as the web server's user. It covers these families, each its own command on the work queue (`warm-galaxy-map`, `warm-lists`, `warm-sectors`, `warm-systems`, `warm-phenomena`, `warm-search`; `queue/commands.py`) that warms that family and nothing else and can run alone (`planetgen command-run`). `planetgen warm` builds the list of them and the wrapper runs them in order as a maintenance run (`runrecord.MAINTENANCE`: its own place in Redis, so it does not block a generation run and the visitors' banner never shows it):

| Family | What it warms |
|---|---|
| `galaxy-map` | the opening tiles and stage (existing `warm_opening_view`), the galaxy shape, sector list, phenomena, superbubbles, shell specs, bright-star status |
| `lists` | the first page of the sector, system, uncharted-system, phenomena, polity and species lists, and each list's counts |
| `sectors` | every sector page (the sector detail), nearest the galaxy centre first, up to `page_cache.warm_max_objects` |
| `systems` | every system page, same order and cap |
| `phenomena` | every phenomenon page, same cap |
| `search` | the search page's empty and default answers and its index counts |

Each family ends with one line, `<family>: <n> warmed (<m> already cached, <k> left out by the cap) in <s> s`, printed by the command, written to the activity log (the admin's log) and returned as the unit's result. `planetgen warm` exits non-zero if any command fails, and says which (the run stops at the failing one, as every wrapper list does); it has no silent fallback (a Redis that does not answer, another warm-up running, or an update running are errors, printed). A cap that leaves objects out says so; those fill on first visit. `planetgen warm lists search` runs just those; there is no `--no-rq` and no `--only`, since the order and the choice live in the wrapper's list and Redis is required. A family warms by visiting the pages through the site's own Flask app (`warmup.build_app`: open during the sector-switch maintenance, not rate limited), so what it caches is exactly what the page asks for; a page that does not answer 200 fails the family. The activity log has `GEN warm.family` or `GEN warm.failed`. `planetgen warm` waits for a running update to finish before it queues (the queue takes no job during one, OPS.53), so `update.sh` can start it in the background.

It runs inside the maintenance flag the update set (OPS.51). A list of every family ends with the command `end-map-rebuild-maintenance`, which ends the sector-switch maintenance after every family has succeeded: it is the list, not a family, that knows the maps are then rebuilt. It registers with the update guard (OPS.53) so a later update lists and cancels it.

## Done when

After `update.sh`, the first visit to the Galaxy Map, a sector, a system, a phenomena or system list and a search answers from the cache. Changing one sector rebuilds only the caches that read it. A test per kind of change proves a stale copy is never served as current.

## Not in this slice

Warming past the cap (a galaxy of millions of sectors warms its nearest `warm_max_objects`, the rest fill on first visit); resolving a star, planet or moon edit to its system (those routes bump every system page, `system-any`, instead); a per-page progress bar on the runs page (the node shows its state and the one result line).

## Without Redis

The shared tier needs Redis (PERF.86 makes it required for generation anyway). With none, the page cache stores nothing (a copy it could not validate is never served), `planetgen warm` refuses with a printed error (Redis is required for every run, PERF.127 removed the inline `--no-rq` mode).
