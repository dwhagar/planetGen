# planetGen API

A JSON API over the planetGen database (`src/planetgen/db/schema.sql`),
built with [Flask](https://flask.palletsprojects.com/). This is `TODO.md`'s
Phase 5 backend — backed by the same MySQL database every other tool in this
project uses (see [`database-schema.md`](database-schema.md) for the MySQL
port). It now has a frontend: the interim `../src/html/` browser
([`html-interface.md`](html-interface.md)) is a thin client over this API
rather than a direct database consumer — see that doc's own note on the
switch, and "Deploying" below for how both are served by one app.

Every read endpoint is fully implemented and unauthenticated (read-only,
no account needed), except the admin stats endpoints. Every write endpoint
(create/modify/delete a sector, system or facility, rename a body,
generate a sector's neighborhood) requires an authenticated admin whose credentials aren't still the
seeded first login — see "Authentication" and "Write endpoints" below.

## Why Flask

Comparison against FastAPI/Django REST Framework: the persistence layer
(`planetgen/db/store.py`) is deliberately plain SQL over a small `pymysql`
wrapper with no ORM, this API is read-heavy with no concurrency pressure yet,
and it needs to deploy onto a plain Apache2/VPS setup (or nginx, Caddy
or macOS; see [`deployment/`](../deployment/README.md)).
Flask has no opinion
about the data layer (route handlers call straight into `planetgen.db.query`'s and
`planetgen.db.store`'s existing functions), deploys as a plain WSGI app
(`mod_wsgi` or gunicorn), and
lives at `../src/planetgen/api/`, mounted at `/api/` by the same Flask app that
serves the HTML pages (`../src/planetgen/web/`, see "Deploying"
below). Those pages are this API's own frontend, calling it in-process. FastAPI's
headline advantages (async, auto-generated OpenAPI docs) still don't pay
for themselves: this API is read-heavy and low-concurrency regardless of
framework, and the pages are plain server-rendered HTML with no use for
generated API-contract docs the way a JS single-page app
would. Worth revisiting if a richer JS frontend is ever built against
this API instead.

## API version and object IDs

Every `/api` answer carries the API version in the `X-PlanetGen-API-Version` header, and `/api/health` reports it as `api_version` (`src/planetgen/api/version.py`, API.22). A breaking change bumps it, the route snapshot `src/tests/fixtures/api_routes.json` moves with it, and it gets a line here and a row in [api-compatibility.md](api-compatibility.md) (API.4: each version, the release that first carried it, the routes and payload formats it changed, and which clients speak it), which is built from `docs/reference/api-versions.json`.

- **11** (API.33): `GET /api/admin/backups` and `GET /api/admin/backups/<name>`: an archive and each of its parts have `size_bytes` where they had `size`.
- **10** (API.37): `GET /api/systems/<id>`: a planet and a moon no longer have `life_chemical`. They have `life_class` (the number of the class that lives on them), `life_class_variant` (its letter), `life_class_best` and `life_class_best_score` (the class they support best and its score, 0 to 9), `life_animal_grade` (0 none to 3 large and active) and `life_complexity`, and `life_class_id` (the class and variant as one id such as `1a`, or `null`) and `life_class_scores` (`{class id: score 1 to 9}` for every class the body supports, empty for none). A species has `life_class` where it had `life_chemical`.
- **9** (MAP.287): `GET /api/systems/<id>/scene`: a star, planet and moon no longer has `color`; it has `look`, the body look record (`bodylook.py`), which holds its `color` and everything the painters need. A comet keeps `color`.
- **8** (MAP.262): `GET /api/galaxy/tiles` is gone: the map's data is read by sector-block node (`GET /api/galaxy/nodes`), and `GET /api/galaxy/changes` no longer has `tiles` (it has `nodes`).
- **7** (API.26): every sector address in a request or an answer is `shell`/`band`/`wedge` (`shell_index`/`band_index`/`wedge_index` in answers), where it was `ring`/`layer`/`slot` (`ring_index`/`layer_index`/`ring_slot_index`). The switch to shells, bands and wedges was made on 2026-10-15 without raising the version, so 4, 5 and 6 were taken by later changes and this one carries it. No route that takes the old names is kept.
- **6** (PERF.99): `GET /api/admin/work` and `GET /api/admin/work/<id>` return the runs in Redis (`items` are runs with their stage, units and time left, not job tree roots; `<id>` is a run id and answers the run with its stages and log), `POST /api/admin/work/<id>/control` takes `cancel` or `clear` (`pause` and `resume` of a node are gone), `POST /api/admin/work/<id>/delete` and `clear-finished` delete runs, `GET /api/admin/work/<id>/log` is new, and `POST /api/admin/work/lease/clear` is gone (clear the run instead).
- **5** (DB.51): `GET /api/systems/<uid>/sections` is gone: the description is a document built from stored data, not a table of paragraphs.
- **4** (ADM.50): `GET /api/admin/galaxy-settings` and `GET /api/admin/galaxy-settings/<name>` are gone with the galaxy settings files.
- **3** (API.23, stage 2): stars, planets, moons, asteroid belts, comets and facilities are named by their object ID too, in the same way: `id`, `star_id`, `planet_id`, `moon_id`, `belt_id`, `comet_id`, `host_id` and `homeworld_planet_id` fields, `<kind>:<ID>` references (including the `around`/`parent` fields of a system scene), the `ref` of every section and link of `/api/systems/<id>/description`, and the URLs (`/api/planets/<id>`, `/api/moons/<id>`, `/api/facilities/<id>`, ...). A facility's `host_id` (on `POST /api/facilities` and `GET /api/facilities/orbit`) is the host's ID: a sector's for `host_type` `space`. A row number is refused.
- **2** (API.23, stage 1): sectors, systems and phenomena are named by their object ID, printed as `0008000000-4000000-000`, and not by their database row number. IDs appear in every `id`, `sector_id` and `ref` field of those answers, in the URLs (`/api/systems/<id>`, `/api/sectors/<id>`, `/api/phenomena/<type>/<id>`, `/api/nebulae/<id>/...`) and in `from`/`to`/`course` references such as `system:0008000000-4000000-000`. A row number is refused (404). A sector's ID is the hex designation of its grid address (shell, band and wedge bit-packed, e.g. `FE81000A2B`), or `(1 << 40) + id` when it was made by hand with no grid address.
- **1**: the first numbered version.

The sector grid changed to shells, bands and wedges on 2026-10-15 (the spherical switch). Every sector address in a request or an answer is now `shell`/`band`/`wedge` (`shell_index`/`band_index`/`wedge_index` in answers) where it was `ring`/`layer`/`slot` (`ring_index`/`layer_index`/`ring_slot_index`), and sector IDs were renumbered from the new designations. Version 6 records it.

## Units of the numbers

Every number the API takes or returns says what it measures. A response field named after a stored column carries its unit as the name's suffix (`radius_km`, `age_gy`, `position_x_mpc`, `period_years`), the same unit the column registry declares (`db/column_units.py`, DB.76). A request body field either has the suffix too (`distance_km`, `phase_deg`, `offset_ly`, `radius_ly`, `edge_ly`, `expires_days`) or a description starting "Unit:" in its schema: a count (`num_orbits`, the worker `count`), an id (`type_id`, `affiliation_polity_id`), an index (`sort_order`), `position` in light-years from the sector's centre, and the facility type field values, whose unit is the one the field itself declares (`unit`). `tests/test_api_field_units.py` fails for a request field that is a number and names no unit. For answers, `tests/test_api_response_units.py` reads the real answers of the sector, system, search and galaxy routes and fails for a number whose key (or the key of the map holding it) neither ends in a unit suffix, nor is a column the registry gives a unit that says what it is (a fraction, a score of 9, an optical depth, a pH), nor is listed in the test with the thing it counts or names (a count, an id, a grid address, an eccentricity, a spin-axis part).

## The User-Agent and the version header

Every request the software makes to another site sends one `User-Agent` in the RFC 9110 form of product tokens with a comment, `PlanetGen/<release> (+<project URL>) PlanetGen-API/<API version>`, for example `PlanetGen/8.0.2017 (+https://github.com/dwhagar/planetGen) PlanetGen-API/7` (API.31). `planetgen.util.user_agent.user_agent()` builds it and is the only place that spells it; the API client the web pages use, the setup-time dataset downloads (`planetgen.datasets.fetch`, [gap-filling-plan.md](../design/research/gap-filling-plan.md) section 6.5) and the build script for the element tables all send it. A page's own JavaScript cannot set `User-Agent`, so `static/apiversion.js` (loaded first by every page) adds the `X-PlanetGen-API-Version` header, with the same API version, to each `fetch` and `XMLHttpRequest` the page makes to this site's `/api/`.

## The API call log

Every request to `/api/...` leaves one `API call` line in the activity log (`planetgen.log`, category `API`, API.15): the time, the account that made it (the key's owner, the signed-in user, or `god` for a command-line tool), how it came in (`via=key`, `session`, `anonymous` or `console`), the method, route, the HTTP response code, the time taken and the sizes. The log uses HTTP response codes, for parity with the web server's own logs. The same id is sent back in the `X-Request-Id` header (a client may send its own, 8 to 64 letters, digits, dots, dashes or underscores), so a client can quote it. Keys, tokens, passwords and request bodies are never written. `api_log.internal` set to `off` leaves out the calls the web pages make in-process. The fields are listed in [config.md](config.md).

## Endpoints

All under `/api/`, all JSON in, JSON out. Every endpoint below except
`/databases` accepts an optional `db=<name>` query
parameter selecting *which* MySQL schema on the configured server to
read from (validated against the same prefix-filtered list
`/api/databases` itself returns — an unrecognized name is a `404`, same
as an unrecognized sector/system id); omitted, it falls back to
`MYSQL_CONFIG`'s own configured default database (`config.py`). This is
what lets one API process serve several databases — see
`planetgen.db.store.list_databases`/`resolve_database`. (The HTML pages
always show the one configured database.) `/api/health` also honors `db=` — passing it checks
connectivity to that specific schema rather than the default one.

### Read

- `GET /api/health` — liveness/readiness check: confirms the process is up
  and the configured (or `db=`-selected) database can actually be opened.
  Returns `{"status": "ok", "schema_version": <n>, "schema_current":
  <bool>, "generator_version": <n>}` (`generator_version` is the one counter
  of the generator rules, OPS.28's generator epoch, bumped only when output
  changes for the same seed; an added field, not a breaking change) (plus `maintenance`: `{"enabled": true, "cause", "reason"}` while the site is closed, since `/api/health` stays open to everyone then so that `update.sh` can check its reload before it opens the site again, OPS.75; plus a `detail` naming the fix when the schema is behind the
  code's `SCHEMA_VERSION` or was never initialized; still a `200`), or
  `{"status": "error", "detail": "database unavailable"}` with a `503` if
  the database can't be reached (the real
  reason, which names the MySQL user and host, goes only to the server's
  log). Rate-limited per client IP (`ratelimit.pages.health`, default 60 per
  minute), not by the default limits; see "Rate limiting".
- `GET /api/databases` — every MySQL schema on the configured server
  whose name starts with this deployment's prefix (literally: `_` and `%`
  in the prefix are not wildcards), except the control schema
  (`control_database`, which holds admin logins and API keys and is never
  listed or selectable with `db=`) (`planetgen.db.store.list_databases`),
  each with `name`, `size_bytes`, `modified_at`, and a quick-glance
  `sector_count`/`system_count` (`null` for a matching schema missing this
  project's own tables, e.g. mid-migration, rather than failing the whole
  listing).
- `GET /api/sectors?limit=<n>&offset=<n>` — every sector, paginated (see
  "Pagination" below), each with `id`, `name`, `edge_mpc`, `edge_ly`,
  `system_count`, and its galaxy placement (`center_x_pc`/`center_y_pc`/
  `center_z_pc`/`galactic_radius_pc`/`galactic_radius_ly`/`shell_index`/
  `band_index`/`wedge_index` (the spherical grid address), all
  `null` together for an unplaced sector, plus a `placed` bool), nearest
  the galactic core first, then unplaced sectors by name
  (`queryDb.list_sectors`/`count_sectors`). The Sectors table (UX.41) also
  takes `sort` (`name`, `systems`, `density`, `position`, `distance`) with
  `order=asc|desc` (unplaced sectors sort last for `density`, `distance` and
  `position`; ties follow the sort's direction by sector ID), the repeatable filter
  `quadrant` (`I`-`IV`, or `unplaced`) and `facets=1`, which adds
  `facets: {"quadrant": [{"value", "count"}]}` (the menu ignores its own
  filter). `total` counts the sectors that pass the filter; a bad `sort` or
  `order` is a 400. `after=<sector ID>` (PERF.137, PERF.208) starts the page after that
  sector in every sort and ignores `offset`, so a page costs the same at any depth: the
  placed sectors are read by the distance index and the unplaced ones after them by name; the
  systems, density and position sorts read the stored count, density and Quadrant
  (`sector_system_counts`) by index, with the sectors that have no galaxy position last, and
  sectors that tie follow the sort's direction by id; a sector made since the last refresh of
  those figures follows the others. An `after` that is not a sector ID is a 400.
- `GET /api/sectors/<id>` — one sector's full display detail: the same
  fields as the listing above, plus `systems` (every system placed in
  it, nearest the sector's center first: `id`, `name`, `quadrant`, `location`, `is_binary`, `binary_type`,
  `position_x_mpc`/`position_y_mpc`/`position_z_mpc`, `center_distance_ly`,
  `runaway_class`/`runaway_speed_kms`, `inside` (the nebula or supernova
  remnant around it, or `null`), `nearest` (its stored nearest systems,
  `{id, name, distance_ly}`), and `stars`, each
  with `role` (`single`, `primary`, `secondary`, or `member` for a star beyond the pair)/`name`/`star_type`/`temperature_k`/`radius_km`/`luminosity_w`),
  `star_count`, `interstellar_debris_count` (an estimate of the
  interstellar comets and planetesimals drifting through it, far too many
  to store as rows), `neighbors` (every grid cell sharing a face with a
  placed sector: its address, `designation`, `direction_pc`, `exists`
  and `sector_id`/`sector_name`, empty for an unplaced sector), and
  `phenomena` (every placed standalone phenomenon generated as part of
  this sector, plus every cloud from elsewhere whose sphere could
  plausibly reach into this sector's cell, nearest its center first; a
  point-like object from another sector is left out — `id`,
  `type` (`"nebula"`/`"asteroid_field"`/`"black_hole"`/`"neutron_star"`/
  `"supernova_remnant"`/`"rogue_planet"`/`"interstellar_comet"`/`"quasar"`), `name`,
  `descriptor`, `radius_ly` (always 0 for a black hole/neutron star/rogue
  planet/interstellar comet/quasar — point-like at this scale), `distance_ly`,
  `offset_x_ly`/`offset_y_ly`/`offset_z_ly`, its center
  relative to this sector's own, `age_years` (a supernova remnant's age, MAP.286; `null` for every other type), and `home` (`true` when generated as
  part of this sector, `false` for a neighbor's cloud) — `queryDb.phenomena_near_sector`, empty
  for an unplaced sector; see `schema.sql`'s "v18"/"v21"/"v28" header
  notes),
  and `stats` (PERF.11, GEN.44: the sector's
  `sector_stats` row — `bright_level_sol`, `relative_density`,
  `expected_systems`, `actual_systems`, `actual_stars`,
  `mean_temperature_k`, `mean_luminosity_sol`, `mean_age_gy` and
  `total_luminosity_sol` (DB.14: what the Galaxy Map colors from) — or `null` for a sector off
  the grid) (`queryDb.sector_detail`). Distinct from
  `planetgen.db.store.load_sector(...).to_dict()`'s *generation* object
  graph (config/provenance, no database ids) — this is the flat,
  ids-and-display-fields shape the sector page's (`/sector/<id>`) Contents table
  and Sector Map actually need.
- `GET /api/systems?star_type=<prefix>&sector_id=<id|none>&limit=<n>&offset=<n>` —
  filtered, paginated system listing (`queryDb.list_systems`/
  `count_systems`), each with `id`, `name`, `sector_id`, `sector_name`,
  `quadrant` (its octant in the sector), `is_binary`, and
  `star_summary` (the single star's `star_type`, or a binary's
  `binary_type`). `star_type` is a literal prefix (`%` and `_` match only
  themselves). `sector_id=none` matches only standalone systems
  (`sector_id IS NULL`, the `/systems` page's Standalone Systems table) — distinct
  from omitting `sector_id` entirely (no sector filter at all). The Systems
  tables (UX.41) also take `sort` (`name`, `sector`, `octant`, `binary`) with
  `order=asc|desc` (systems with no sector or octant sort last), the filters
  `binary=yes|no`, `placement=sector|standalone` and `octant` (repeatable),
  and `facets=1`, which adds `facets: {"placement", "binary", "octant"}`
  (each `[{"value", "count"}]`; every menu's counts apply the other filters
  but not its own). `total` counts the systems that pass the filters; where
  no stored count exists yet it stops at 10,000 and `total_capped` is `true`
  (the exact figure follows once the background count finishes). `after=<system
  ID>` starts the page after that system in the name, sector, octant and binary sorts and ignores
  `offset`, so page 40,000 costs what page 1 does (PERF.75, PERF.206; an unknown ID is a 400);
  `from=<text>` (PERF.75, the Systems table's jump-by-value box) begins the
  list at the first name from that text on under the name sort (`M` starts
  at the first name from M; under `order=desc` the names that begin with M
  come first and the list runs down from them) or at the first sector from
  it on under the sector sort (systems with no sector still come last),
  `total` counts the systems from there, and `offset` and `after` count from
  there; the other sorts ignore it;
  a bad `sort`, `order`, `binary` or `placement` is a 400.
- `GET /api/uncharted-systems?sort=<luminosity|temperature|type|sector|age>&order=asc|desc&limit=<n>&offset=<n>&after=<star id>&from=<value>` —
  (`after` starts the page after that star in the sort and ignores `offset`, PERF.75; equal sort values are ordered by id in the
  sort's own direction, so the luminosity and sector sorts walk their index with no sort; a non-number is a 400;
  `from=<value>`, PERF.135, begins the list at the first star at or after a luminosity in solar luminosities, a temperature in kelvin, an age in Gyr, a
  shell number or the start of a star type, in the sort's direction, `total` counts from there and `offset` and `after` count from there; a value
  the sort cannot read is a 400, and one past the last star answers the last page, `past_end` true and `total` its length)
  the stars the brightness scatter placed that have no system yet (UX.87:
  `bright_stars` rows with a NULL `star_system_id`, so both the mass-limit
  and luminosity passes count). Each item has `id`, `designation`,
  `shell_index`, `band_index`, `wedge_index`, galaxy-frame `x`, `y`, `z`
  and in-sector `local_x`, `local_y`, `local_z` (parsecs), and the star's
  own data (`star_type`, `population`, `mass_solar`, `radius_solar`,
  `temperature_k`, `luminosity_sol`, `age_gy`, `galactic_radius_pc`).
  Brightest first unless `sort`/`order` say otherwise; a bad `sort` or
  `order` is a 400. Each item also has `member_count` (0 for a single star) and
  `system_mass_solar` (the mass the point carries, `null` for a single star): the stored
  members of a star that has company (UX.114).
- `GET /api/object-lists/<kind>` — one keyset page of the objects of a kind from the object index (UX.141): `kind` is
  `system`, `star`, `planet`, `moon`, `belt`, `comet`, `facility` or a phenomenon type. `sort` (`name`, the default,
  `class`, `mass` or `radius`; the last three list the objects that have one) with `order` (`asc` or `desc`), `limit`
  (50 by default, 200 at most), `after=<object ID>` (the `next` of the previous page; there is no `offset`), `name`
  (names starting with it, case not significant), `class` (exactly this class), `sector_id` and `from` (begin a
  `sort=name` list at the first name from it on; `past_end` is true when none lies at or past it and the last page is
  answered). Answers `items`
  (`kind`, `id` the printed object ID, `source_id`, `name`, `sector_id`, `system_id`, `parent_id`, `object_class`,
  `mass_kg`, `radius_km`, `placed`; a planet's or moon's also `phi4` and the six `tier_<domain>` PHI-6 ratings, 0 Blue to 3 Red, `null` for a body that was not scored, UX.100; a star's also `temperature_k`, `luminosity_w` and `age_gy`, a planet's or moon's `distance_km` (its orbit size), `eccentricity` and `surface_temperature_k`, as stored, `null` where the row holds none, UX.187), `next` (null at the last page), `total` and `total_capped`, `limit`. The index is
  what `planetgen warm object-index` last wrote, so a kind it has not been built for lists nothing. 400 for an unknown
  kind, sort or order, a `limit` or `sector_id` that is not a number, a `from` under another sort, or an `after` that is not an object ID
  (`objectlist.list_objects`).
- `GET /api/uncharted-systems/<bright star id>` — one scattered star (UX.114): the fields of a list item,
  `rotation_period_hours` (GEN.846: its spin period; `null` on a row made before schema v185), `system_id` (the system built around it, `null` while it is uncharted), `hosted_nebula` (the nebula the star lights, the one
  its sector builds: `name`, `nebula_class`, `nebula_type`, `class_name`, `radius_ly`, `age_myr`; `null` for a star that
  lights none or is built, GEN.530) and, for a star with company that is not built, `members` (`member_no`, `kind`, `state`, `mass_solar`, `parent_node_no`, `flags`, `radius_au`: GEN.952, `null` on a row made before schema v184) and `nodes` (`node_no`,
  `parent_node_no`, `kind` of `close`, `wide`, `compact` or `group`, `mass_solar`, `separation_au`,
  `period_years`, `eccentricity`, `interaction`, `flags`) and `compact_pairs` (UX.186: one per accreting neutron star or black hole
  with a companion or ablating pulsar with its companion: `node_no`, `kind` of `high-mass-xray-binary`, `be-xray-binary`,
  `low-mass-xray-binary`, `black-widow` or `redback`, `compact` and `companion` (member numbers), `microquasar` and `x_ray_luminosity_w`: GEN.952, the watts an accreting pair radiates, worked out on read from the stored masses, radii and orbit, `null` for a spider pair and for a tree stored before v184). 404 for an id that is no row (`queryDb.uncharted_system_detail`).
- `GET /api/systems/<id>` — one system's full display detail: `id`,
  `name`, `sector_id`, `quadrant`, `location`, `is_binary`, `binary_type`,
  `binary_configuration` (`"close"`/`"wide"`/`null`),
  `binary_mutual_position_x/y/z_km` (the secondary's position relative to
  the primary — `null` for a single star; the System Map's own real
  binary-star placement is derived from this plus each star's `mass_kg`),
  `runaway_class` (`"runaway"`/`"hypervelocity"`/`null`) and
  `runaway_speed_kms`, `stars` (each with `parent_node_no` and `slot`, the node of the orbit tree it hangs in and which of its two places it holds), `orbits` (MAP.213: the
  system's tree of orbits, one entry per node with `node_no`, `parent_node_no`, `slot`, `kind`, `mass_solar`, `separation_km`, `eccentricity`, `inclination_deg`, `period_years`, `interaction` and `flags`; empty for a single star), `planets` (each with
  its own nested `moons`; every planet and moon also carries `habitable` (its overall PHI-6 colour is Blue or Green, GEN.643),
  `habitable_class` (its class can carry life), `life_class_id` and `life_class_scores` (API.37: the life class that lives on it, such as `1a`, and `{class id: score}` for each class it supports; `life_class_id` is `null` for a world that suits a class but whose life has not started, a star younger than the first step of its schedule, GEN.1001, where the stored `life_class` stays), `life_stage` — the furthest evolutionary milestone its timeline reached,
  or `null` — and `inhabited`; a scored rocky one carries `equipment` too — the survival equipment a human needs, GEN.655: `required` and `recommended` item names, a short `code` such as `B4p F100 S2 H3 [R2]`, `required_line` and `recommended_line` (UX.163; a line with nothing in it is empty), the `items` (`slot`, `item`, `severity`, `domain`), the `line` as text, `shelter_column_g_cm2`, `outside_hours` (the hours a year the dose allows), `outside_hours_day` (the hours of one local day a person can be outside, GEN.768) with `outside_hours_limit` (what sets them: `daylight`, `dose`, `heat`, `cold`, `air`, `uv` or `vehicle`) and `regions` (GEN.656: the regional temperatures in kelvin by short name — `polar`, `north`, `equator`, `south` and the day and night `swing` for a world that turns on its own, `day_side`, `terminator` and `night_side` for a tidally locked one; GEN.848 adds the swing of every other place by the same short name with `_swing`: `polar_swing`, `north_swing`, `south_swing`, `mean_swing` and, for a locked world, `day_side_swing`, `terminator_swing`, `night_side_swing` (0, it has no day); `null` when it has none) — and the raw `eq_*` levels; the old `equipment_tier` is gone; `null` for a gas giant, an unscored body and a row not yet re-scored), `belts` (each with its `composition`,
  `{component, concentration}` largest share first), `comets`, `sector_siblings`
  (`{id, name}` for every system in the same sector), and
  `nearest_neighbors` (`{id, name, distance_ly}` for the up-to-3 closest
  systems, nearest first, from the stored `nearest_systems` rows (schema
  v41), which are searched across sector boundaries; a system with no
  stored rows falls back to its own sector's closest systems, computed
  from current positions and names. The names inside `location` are
  frozen at generation time and can be stale), `inside` (the nebula or
  supernova remnant around the system, `{type, id, name, class,
  density_cm3, temperature_k}`, or `null`), and
  `heliopause_au` (the heliopause squeezed by that cloud, the edge of the
  system's navigation frame) with `heliopause_open_space_au`
  (`queryDb.system_detail`) — same "flat display
  shape, not the generation object graph" relationship to
  `planetgen.db.store.load_star_system(...).to_dict()` as `/api/sectors/<id>`
  above.
- `GET /api/systems/<id>/description` — the system's description (API.25,
  `planetgen/describe/`, `docs/design/architecture/system-description.md`): `{"kind":
  "system", "ref", "title", "sections": [...]}`. A section has a `kind`
  (`system`, `star`, `pair`, `planet`, `moon`, `belt`, `comet`), a `ref`
  (`planet:<ID>`, or `null` for a pair), a `title`, `blocks` and `children`.
  A block is `facts` (rows of `key`, `label`, `value`), `ratings` (the PHI-6
  `overall`, the `gear` a visitor needs and a row per category), `composition`
  or `prose` (sentences). Text is a list of strings and typed items: a
  quantity `{"qty": {"kind", "si", ...}}` (`si` in the kind's base unit, never
  formatted), a link `{"ref", "name"}`, `{"strong"}`, a chemical formula
  `{"chem"}` or a colour `{"swatch": {"hex": "#99bfff", "word": "pale blue"}}`.
  A planet or moon with a solid surface and air has an `atmosphere_colour`
  fact whose value is a swatch (the clear-sky colour of its air under its
  star, [atmosphere-colour.md](atmosphere-colour.md); absent for an airless
  body and for a giant). Built from the database rows on each request.
- `GET /api/systems/<id>/text?format=markdown|text` — the same
  description written out on request, `{"id", "format", "content"}`
  (`planetgen/db/render.py`; no page text is stored since schema v29).
  `format` defaults to `markdown`; any other value is a `400`.
- `GET /api/planets/<id>/life-readings` and `GET /api/moons/<id>/life-readings` — the three
  words for each kind of life on a world (UX.251, `physics.life_readings`), worked out on
  request from the stored body and its star and stored nowhere. `{"id", "life_class_id",
  "readings"}`. `readings` is empty for a gas giant and for a body whose kind scores are not
  stored; otherwise it has one entry for each of microbes, flora, fauna and sentient life, in
  that order: `kind` (`microbial`, `flora`, `fauna` or `sentient`), `label`, `score` (the stored
  0 to 9), `abundance` (Absent, Rare, Occasional, Frequent, Abundant or Dominant), `trend`
  (Growing, Stable or Declining), `security` (Secure, Marginal, Precarious or Dying out),
  `years_left_gyr` (Gyr until the living class stops scoring on the world as its star
  brightens; `null` when it lasts or is not known) and `text` ("Microbes: abundant, stable,
  secure."). The trend and security of a kind that is Absent are `null`. An unknown id is a `404`.
- `GET /api/planets/<id>/life-air` and `GET /api/moons/<id>/life-air` — the data of
  the Life and air chart on a planet's or moon's page (UX.237,
  `physics.life_air.air_history`): the air through the age of its star as far as the
  life class that lives on it changed it, worked out on request from the stored air and
  stored nowhere. `{"id", "life_class_id", "life_air"}`. `life_air` is `null` for a world
  with no living class, a giant, a class that changed nothing, or an air the life did not
  change; otherwise `age_gyr` (the star's age), `producers_gyr` (when the producers began),
  `rose_gyr` (when the oxygen rose, `null` for a class that makes none), `gases` (`{"formula",
  "name", "steps": [[age_gyr, partial_pressure_kpa], ...]}` for oxygen, carbon dioxide and methane,
  each value holding until the next age; a gas that is never above 0.1 Pa is left out),
  `stages` (`{"key", "age_gyr"}` for each milestone the star is old enough for) and
  `speculative` (the ledger's `{"tier", "tweak"}` for a class whose numbers are tuned, else
  `null`). An unknown planet or moon is a `404`. Additive: the API version does not change.
- `GET /api/systems/<id>/scene` — the 3D system view's data (MAP.69):
  `system`, `epoch`/`epoch_unix` (when the stored phases were last
  advanced, `null` if never), `stars`, `planets` (each with `moons`),
  `tree` (the tree of orbits the stars hang in, GEN.311: `nodes`, one per
  `orbit_nodes` row with its own epoch, `members`, one per star with its
  `parent_node_no`, `slot` and `mass_solar`, and `origin`, the star at the
  origin of a wide pair or `null`), `belts` and `comets`. Every object has a `ref` (`star:3`, `planet:12`,
  `moon:40`, `belt:5`, `comet:8`), `radius_km`, `look` (a star, planet or moon) or `color` (a comet), `position_km`
  relative to what it orbits, and an `orbit` (`around`, `distance_km`,
  `period_years`, `inclination_deg`, `ascending_node_deg`, `phase_deg`; a
  comet's holds its Kepler elements under `kepler`). A planet or moon also carries `habitable`, `inhabited`,
  `body_type` (`t` terrestrial, `g` gas giant), `phi4` (its overall habitability score, `null` unscored) and `tier`
  (that score's colour, `Blue` to `Red`, or `null`), which the Galaxy Map draws its symbol from (MAP.321), and the four life scores `life_microbial`, `life_flora`, `life_fauna`, `life_sentient` (0 none to 9 abundant, `null` unscored) with `life_stage` (kept for a class that can carry life only) and `planet_class`, which its symbols for life, sentient life, a civilization and a planet class read (MAP.343).
- `GET /api/near?from=<ref>&distance=<pc>` (or `point=x,y,z` in galaxy-frame
  parsecs) — everything generated within `distance` parsecs (at most 50) of
  a place (NAV.43, `planetgen.db.near.objects_within`): systems, stars,
  planets, moons, belts, comets, facilities and the standalone phenomena,
  nearest first, each with `ref`, `kind`, `id`, `name`, `parent` and
  `distance_pc`. Optional `kinds` (comma-separated), `limit` (50, at most
  200) and `offset`. The answer also carries `total`, `by_kind`, and
  `not_generated` (`{kind: count}`, with `not_generated_total`): the
  pre-placed bright stars (`bright_star`) and scatter objects in the sphere
  that are not built yet, counted and never listed; the search never
  generates and never reads the sectors (NAV.60: it reads each table by its
  cell key). `nearest=K` (NAV.61) asks for the nearest K rows instead: the
  search tries 2, 4, 8, 16 and 32 pc and stops at the first sphere that
  holds K, so it is quick at the dense core too, and `distance` is optional
  and caps how far it looks. A sphere that holds more than 100,000 systems
  is refused (400) with the text that says to shrink the distance. The
  answer carries `systems_in_range` (the systems in the sphere, counted no
  further than 10,000) and `systems_in_range_capped` (`true` above that),
  which the What's nearby page shows as "10,000+".
- `GET /api/objects/<ref>` — one object by its reference (NAV.7), the one
  form every kind shares: `<kind>:<id>` with kind `sector`, `system`, `star`,
  `planet`, `moon`, `belt`, `comet`, `facility` (UX.183: its parents run down to the body it is built on, or to its
  sector, and its siblings are the facilities of that host) or a phenomenon type (`nebula`,
  `asteroid_field`, `black_hole`, `neutron_star`, `supernova_remnant`,
  `rogue_planet`, `interstellar_comet`, `quasar`); a bare number is a system.
  Returns `ref`, `kind`, `id`, `name`; `parents` (from `galaxy` down to the
  direct parent, each `{ref, kind, name}`); `siblings` (references of the
  other objects of the same kind under the same parent, at most 200, with
  `siblings_truncated`); and `positions`: `galaxy_pc` (parsecs),
  `sector_ly` (light-years from the sector's center) and `system_km`
  (kilometers from the system's origin), each `[x, y, z]` or `null` where
  the object has none (a body inside a system takes the system's galaxy and
  sector positions; a belt is a ring with no `system_km`; a moon of a free planet (DB.94) has `rogue_planet` as its parent after the galaxy, the rogue planet's galaxy point, no `system_km`, and the planet's other moons as siblings). A bad reference
  or a missing row is a `404`. A planet or moon also carries `phi4` and the six `tier_<domain>` ratings (UX.100) and its four life scores `life_microbial`, `life_flora`, `life_fauna` and `life_sentient`, whole numbers 0 (none) to 9 (abundant), `null` for a body that was not scored (UX.152). `keep_out` (NAV.24) is the object's keep-out
  radius for course planning: `{radius_km, basis, note}`, `basis` one of
  `hill` (a planet or moon), `perimeter` (a star or system: the widest
  stored `system_perimeter_km`), `galactic_hill` (a black hole, neutron
  star, quasar or rogue planet, never smaller than the object itself), `radius`,
  or `none` (`radius_km` is `null`: nebulae, supernova remnants and asteroid
  fields pass through with a `note`; belts, comets and sectors have none).
  `queryDb.resolve_object`,
  `planetgen.galaxy.objectref` and `static/objectref.js` hold the logic.
- `GET /api/objects/<ref>/sheet` — the property sheet of a star, planet, moon,
  belt, comet, facility (UX.183, with a `ratings` group) or standalone phenomenon (UX.168, UX.175, `describe/sheet.py`):
  `kind` and `groups` (`identity`, `orbit`, `physical`, `rotation`,
  `composition`, `atmosphere`, `climate`, `water`, `magnetism`, `radiation`,
  `habitability`, `life`, the star's `zones`, `stellar`, `galactic`, `disk`, and
  a phenomenon's `gas`, `compact`, `activity`), each `{key, label, rows}`. A
  row is `{key, label, value, text, raw, stored, absent, source, show_absent}` (`source`: a one-line source for a physics row, how its value is worked out, `describe/sheet_sources.py`, UX.187; `null` for a row with none written) (`raw`: the stored number, `null` for text and for a column with no value) (`absent`: the value is none, 0, false or empty; `show_absent`: a page lists it even then, UX.188, because its absence matters, such as no air or no magnetic field): `value` is the Document's inline
  list (a quantity is `{"qty": {"kind", "si"}}`, shown in the reader's units),
  `text` the same as plain text, and `stored` false, with the text "not
  stored", for a property the row holds no value for. Every stored column of
  the kind's table is either a row or on the registry's `NOT_SHOWN` list with
  its reason. A reference that names nothing, or a kind with no sheet (a
  sector, a system, a phenomenon), is a `404`.
- `GET /api/objects/<ref>/children[?limit=&after=]` — what belongs to an
  object, one keyset page at a time (UX.169, `db/children.py`). A sector holds
  its systems, then the phenomena placed in it; a system its stars, then the
  planets, belts and comets that orbit the whole system (a close pair's
  circumbinary bodies); a star its planets and belts in their order out from
  it, then its comets; a planet, or a free (rogue) planet, its moons. A moon,
  belt, comet or any other phenomenon holds nothing. The answer is `ref`,
  `kind`, `total` (every child it holds), `items` and `next`. An item is
  `{ref, kind, id, name, summary, counts, child_count}`: `summary` is what its
  own row says (a star's `star_type`, `mass_kg`, `luminosity_w`; a planet's or
  moon's `planet_class`, `distance_km`, `period_years` and the PHI-6 `phi4` and
  `tier_<domain>` ratings and the four `life_<name>` scores (UX.152); a system's `is_binary`, `quadrant` and `star_types`;
  a belt's `density` and limits; a comet's orbit), `counts` how many of each
  kind lie anywhere below it (a system's `star`, `planet`, `moon`, `belt`,
  `comet`; a star's or planet's; none for a sector) and `child_count` how many
  it holds directly, so a row can show its expand arrow without another
  request. `limit` is 1 to 200 (default 50); `after` is the previous page's
  `next`, an opaque token (`null` at the end), and a page never uses `OFFSET`.
  A reference that names nothing is a `404`; a bad `limit` or `after` is a `400`.
- `GET /api/nav?from=<ref>&to=<ref>` — course, distance, and an optimal route
  between two objects (systems, bodies, phenomena) (`queryDb.nav_between`) — see "NAV" below.
- `GET /api/galaxy/sectors` — every galaxy-placed sector (non-`null`
  galaxy placement), each with `id`, `name`, `x`/`y`/`z`
  (`center_x/y/z_pc`), `galactic_radius_pc`, `shell_index`, and
  `system_count` (`queryDb.galaxy_placed_sectors`) — the data
  the Galaxy Map page (`/galaxy`, `../src/planetgen/web/galaxy_views.py`) plots. Not paginated: bounded by
  how much of the galaxy has actually been generated (the lazy galaxy-scale generation
  design, finished in the original roadmap's phase 4), not by the addressable galaxy's own
  scale.
- `GET /api/galaxy/phenomena` — every galaxy-placed standalone
  phenomenon of all eight types
  (non-`null` `center_x/y/z_pc`), each with `id`, `type`
  (`"nebula"`/`"asteroid_field"`/`"black_hole"`/`"neutron_star"`/
  `"supernova_remnant"`/`"rogue_planet"`/`"interstellar_comet"`/`"quasar"`), `name`,
  `descriptor` (a nebula's `nebula_type`, a field's `density`, a black
  hole's accretion state, a neutron star's `pulsar_type`, a remnant's
  `morphology`, a rogue planet's body type, a comet's activity, or a
  quasar's radio loudness), `radius_ly`
  (0 for every point-like type: everything but nebulae, asteroid fields
  and supernova remnants), `x`/`y`/`z`
  (`center_x/y/z_pc`), and `galactic_radius_pc`
  (`queryDb.galaxy_placed_phenomena`) — the phenomenon counterpart to
  `/api/galaxy/sectors`, plotted as a small fixed-size dot on the same
  Galaxy Map (a nebula/asteroid field's own real physical extent is
  instead shown as a translucent cloud on the Sector Map of any sector it
  reaches into — a black hole/neutron star stays point-like there too — see
  `/api/sectors/<id>`'s `phenomena` key above). Not paginated, for the
  same reason `/api/galaxy/sectors` isn't.
- `GET /api/galaxy/shape` — `{"shape": ...}`, the galaxy's stored
  density-skeleton shape (`planetgen plan`'s output): every
  `planetgen.galaxy.density.GalaxyShape` field plus `edge_pc`,
  `outer_shell_index`, and `expected_system_count_at_density_1`, with
  `pinned_stellar_mass_solar` (the stellar mass the plan pinned, else `null`)
  and `pins` (DB.81: `{lever name: value}` of the shape levers the plan was
  told, empty when the seed drew everything) (`queryDb.galaxy_density_shape`). `shape` is `null` when the skeleton
  has never been built. The Galaxy Map shades its "expected density"
  cloud from this real model (falling back to a generic illustrative
  gradient when `null`) instead of a placeholder.
  `version_warning` (DB.7) is a sentence when the galaxy holds sectors generated by another PlanetGen version than the running one, else `null`.
  `bright_stars` (`queryDb.bright_star_scatter_status`) is the scatter's floor and seed, and its `limits` (UX.116, `null` before a plan) are the settings stored with the plan that the Generate page shows: `core_min_luminosity_sol` and `core_mass_limit_sol` (the core's own, else the galaxy-wide ones times 1.1), `core_follows` (neither was set), `halo_min_luminosity_sol`, `cell_star_ceiling` and `fill_luminosity_cap_sol`.
- `GET /api/galaxy/superbubbles` (GEN.254) — `{"superbubbles": [{"x", "y", "z", "radius_pc"}, ...]}`, the
  cavities the supernovae of OB associations blow in the gas, in galaxy-frame parsecs, biggest first; drawn from the
  galaxy's seed (nothing is stored) and empty before a plan. The Galaxy Map draws them as faint shells. A sector's
  detail (`GET /api/sectors/<id>`) carries `superbubble`: the cavity its centre lies in (`radius_pc`,
  `ambient_density_cm3`, `temperature_k`, `centre_distance_pc`) or `null`.
  It also carries `region` (UX.231): the charted region its centre lies in, `{name, parent}` ("Vang Zarh Region", "Vang Reach"), or `null` when
  that cell holds no generated sector or no naming key is in force. A system's detail (`GET /api/systems/<id>`) and an object (`GET /api/objects/<ref>`) carry
  `region` the same way, for the system's or object's galaxy position.
  It also carries `unresolved` (GEN.281): `null` for a cell filled whole, else what a cell above the cell star
  ceiling leaves unresolved (`ceiling_systems`, `luminosity_cap_sol`, `expected_systems`, `unresolved_systems`,
  `unresolved_mass_solar`), worked out from its address and never stored. It also carries `scatter_tier` (GEN.333): `null`
  outside the galaxy or before a scatter ran, else the tier of the sector's band (`core`, `baseline` or `halo`) and the
  limits the galaxy scatter drew its stars at there (`min_luminosity_sol`, `mass_limit_sol`, the latter `null` when the
  scatter had no mass pass).
- `GET /api/galaxy/cell?shell=<i>&band=<j>&wedge=<k>` (or `?x=&y=&z=`,
  parsecs, for the cell holding that point) — one cell of the spherical
  sector grid, whether or not anything was generated there
  (`spherical.describe_sector_cell`): `shell_index`/`band_index`/
  `wedge_index`, `designation`, `cartesian_pc`, `cylindrical`
  (`r_pc`, `theta_rad`, `z_pc`), `spherical` (`r_pc`, `theta_rad`,
  `polar_rad` from galactic north), `bounds` (`r_inner_pc`, `r_outer_pc`,
  `colatitude_north_rad`, `colatitude_south_rad`, `theta_start_rad`,
  `theta_end_rad`), `mean_band_arc_pc`, `mean_wedge_arc_pc`, `volume_pc3`,
  `vertices_pc` (8 corners), `edge_pc` (the stored skeleton's shell depth,
  else the standard 4 pc), `sector_id` (the generated sector there, or
  `null`) and `in_galaxy` (whether the cell lies inside the planned
  galaxy's stored outline; `null` before any plan). `400` for a bad or
  incomplete query, a negative shell, an address that is not a cell, or a
  cell too far out to compute.
- `GET /api/galaxy/profile?shell=<i>&band=<j>&wedge=<k>` (or `?x=&y=&z=`, parsecs, for the cell holding that point) — the
  Sector Profile of any cell of the grid, charted or not (API.27), so an uncharted sector can be previewed. Read only:
  it is closed form in the stored skeleton and the address (`galaxy.profile.sector_profile`), draws no random number and
  stores nothing, so two calls give the same answer. The body has `shell_index`/`band_index`/`wedge_index`, `center_pc`,
  `edge_pc`, `volume_pc3`, `region` (`nuclear`, `bulge`, `disk`, `halo`, `cluster` or `satellite`), `density` and
  `expected_systems` (the mean of the cell's Poisson draw), `mix` (population -> density), `populations` (for each of
  `young`, `intermediate`, `old`, `bulge`, `halo`, `cluster`, `satellite`, `nuclear_cluster`, `nuclear_disc`, `satellite_young`, `satellite_intermediate`: `density`,
  `share`, `age_gy` as `[youngest, oldest]`, `fe_h` as `[weight, mean, sigma]` peaks, `dispersion_kms`, `rotation`),
  `gas` (`cold_fraction`, `ambient`, `ambient_density_cm3`, `mass_solar`, and `superbubble`: `null` or the cavity's
  `radius_pc`, `ambient_density_cm3`, `temperature_k`, `center_pc`, `centre_distance_pc`), `dark_matter`
  (`density_solar_pc3`, `mass_solar`), `kind_factors` (kind -> the factor the place applies to that kind's rate per star), `type_ia_share` (the share of the place's supernova remnants that are Type Ia; 1 where no young stars are),
  `sector_id` (the generated sector there, or `null`), `in_galaxy` (as in `/api/galaxy/cell`) and `unresolved` (as on a
  sector: what a cell above the cell star ceiling leaves unresolved, else `null`) and `estimate` (UX.113: what filling this one
  cell would make and take, from this server's stored speeds and sizes: `systems`, `stars`, `bytes`, `seconds`, and
  `measured`, whether the speed is this server's own). `404` before a galaxy is planned,
  `400` for a bad or incomplete query, a negative shell, a cell that does not exist or one past the grid's outermost
  shell.
- `GET /api/galaxy/uncharted?shell=<i>&band=<j>&wedge=<k>` — one cell's
  place and what the scatters left in it, for opening a cell nothing was
  generated in (MAP.162): `shell_index`/`band_index`/`wedge_index`,
  `designation`, `center_pc`, `edge_pc`, `sector_id` (the generated sector
  there, or `null`), `stars` (its waiting bright stars,
  `queryDb.bright_stars_in_sector`), the unbuilt scattered phenomena
  (`queryDb.uncharted_sector_contents`) and `hosted_nebulae`, the nebulae its waiting stars light (`star_id`, the
  star's `x`, `y`, `z` in parsecs and the fields of `hosted_nebula` above, GEN.530). `400` for a missing or
  non-integer address, a negative shell, or a cell outside the grid.
- `GET /api/galaxy/made?since=<unix>[&until=<unix>]` — the galaxy-placed
  sectors created in a time window, for the Galaxy Map's "Show on Galaxy Map"
  after a generate job (ADM.31): `{"total", "items"}` (`queryDb.sectors_made`),
  `items` the first `MADE_SECTOR_LIMIT` by id, each `{id, name, x, y, z,
  shell_index, band_index, wedge_index}`. `400` unless `since` (and `until`
  when given) are numbers.
- `GET /api/galaxy/nodes?keys=<key>,<key>,...` — the single-space Galaxy Map's data, one sector block ("node") at a
  time (MAP.260, `db.nodes.galaxy_nodes`). A key is `m.shell.band.wedge` of the block ladder (`m` one of 177147, 59049,
  19683, 6561, 2187, 729, 243, 81, 27, 9, 3, or 1 for a sector), at most 128 keys. Returns `{"nodes": {"<key>":
  {"aggregate": {...}, "stars": [...], "generated": [...], "points": [...], "clouds": [...], "placed": [...]}},
  "edge_pc": ..., "has_shape": ...}`. `aggregate` is what is stored about the block whole: `generated` and `stats` (the
  sums of its generated sectors, `galaxy_block_stats`) and `lights` (the star and point count, luminosity, centre of
  light and mean temperature of the pre-placed objects nothing was built from, `galaxy_block_lights`). The lists are
  the brightest objects inside the block, nested down the ladder: a block lists every object of its parent that lies
  inside it. `placed` is only on a level-3 block and a sector. The nebula cover is not in a node (it changes without a
  sector changing; read `/api/galaxy/nebula-cover`). The node of a charted region (a 9-block that holds at least one generated
  sector, `region_stats`) also has `name` ("Fardoo Beez Region") and `parent` (the name of its 81-block family, "Fardoo Reach"), from the
  gated phoneme codec under the naming key in force (API.35); no other node has them, and none has them while no key is in force. The same node has `extent` (MAP.319), `{"centre": [x, y, z], "radius"}` in parsecs: the bounding sphere of the region's largest face-connected clump of generated sectors, at least the 12 pc framing sphere, which the map flies to; it needs no key; and `charted` (MAP.314), `{"sectors", "percent"}`, how much of the region is generated. A node of 27, 81 or 243 sectors that holds charted regions has `regions` (MAP.314), a list of `{"key", "name", "parent", "centre", "radius", "sectors", "percent"}` per charted region inside it, most sectors first (names need a key in force); no other node has it. `400` for a malformed key or too many keys.
- `GET /api/galaxy/charted-regions` (MAP.316) — every charted region (a 9-block that holds at least one generated sector, `region_stats`), by share charted then cell: `{"regions": [{key, m, shell, band, wedge, centre: [x, y, z], radius, sector_count, capacity, share, clumps, name?, parent?, inhabited?}], "edge_pc"}`. `key` is the cell's node key (`9.shell.band.wedge`), `centre` and `radius` the extent in parsecs (as the node's `extent`), `share` the sectors generated over the sectors the cell's outline holds (0 to 1, above 0 for every region), `name` and `parent` as on the node (API.35, only with a naming key in force), `population` (MAP.318: the people counted in the region, `db/population_slot.py`; present only on a region that has a count, so absent on every region until the world population simulation, POP.30, stores head counts), `inhabited` (MAP.321) only on a region with an inhabited world: `{systems, worlds, score?, tier?}`, the systems with an `inhabited` interest in the region, their inhabited worlds, and the best PHI-6 score and its tier when one is scored (read from `system_interests` when asked, so a colony just founded is in it). The Galaxy Map draws its galactic-scale region spheres and its Menu's charted-regions list from it, through the web page's `GET /galaxy/regions`.
- `GET /api/galaxy/nebula-cover?at=<block key>[&g=<n>]` — the share of each cell of a Galaxy Map stage (`at` alone,
  or none for the galaxy) or data window (`at` the middle cell and `g` the cells across) that lies inside a nebula
  (`queryDb.galaxy_nebula_cover`, MAP.174), for the cells that hold a generated sector and the ones that do not: the
  stored aggregates of MAP.172 are made from the seed's nebula field and the stored nebulae, so they need no sector.
  Answers `{"at", "child_m", "cells"}` with `cells` the `[shell, band, wedge, cover]` of every cell whose cover is
  above zero (`cover` 0 to 1), sorted. The cells of a sector have none (the aggregates are of the blocks above a
  sector), so their list is empty. `/galaxy/stage` and `/galaxy/window` add it to their payloads.
- `GET /api/galaxy/stage?at=<block key>` and
  `GET /api/galaxy/window?at=<block key>&g=<n>` — the Galaxy Map's drill-down
  stage (`queryDb.galaxy_stage`) and data window (`queryDb.galaxy_window`,
  MAP.151b, MAP.205): how many generated sectors each child block holds, with
  the `sector_stats` sums. A block key is `m.shell.band.wedge` with `m` one of
  243, 81, 27, 9, 3 (a sector is `1.shell.band.wedge`), as
  `design/galaxy/spherical-block-ladder.md` defines the ladder (`galaxy/drill.py`).
  A stage answers `{"at", "edge_pc", "child_m", "children", "sectors"}`: `at`
  is the canonical key (`null` when `at` is left out: the galaxy, which lists the
  level-243 blocks that hold a generated sector), `children` is
  `[{"shell", "band", "wedge", "generated", "stats"}]` for each child block
  one size down that holds a generated sector (`stats`:
  `{"systems", "expected_systems", "stars", "mean_age_gy",
  "luminosity_sol", "mass_kg"}`; `mass_kg` is the total mass of the generated contents, from `sector_stats.total_mass_kg`, which the Galaxy Map ranks sector names by, MAP.322), and, when the children are sectors (`child_m` 1),
  `sectors` lists each generated one: `{"shell", "band", "wedge", "id",
  "name", "system_count", "stats"}` (otherwise `null`). When the children are 9-blocks (`child_m` 9) each child also has `name` and
  `parent` (the region and its family, as on the node, API.35), and the stage of a charted 9-block has them itself. Above the sector the
  figures are read from `galaxy_block_stats`. A window is a ball: `g` (3, 5, 7,
  9 or 11) cells of the key's size across, every block whose center lies within
  `g * m / 2` sector edges of the center block's center. It answers `{"at",
  "g", "child_m", "cells", "children", "sectors"}`: `cells` is every cell of the
  ball (`[[shell, band, wedge], ...]`), `children` is the ones that hold a
  generated sector, as a stage lists them, and `sectors` is `null`; a
  window of sectors (`m` 1) is refused, because its names come from the stages of
  its parents. A malformed key or width is a 400.
- `GET /api/galaxy/locate?q=<part of a name>` — the Galaxy Map address
  bar's name lookup (`queryDb.galaxy_locate`): `{"matches": [{"kind"
  (`"sector"`, `"system"`, `"planet"` or `"moon"`), "ref", "id", "name",
  "parents" (`[{"ref", "kind", "name"}]`, the sector down to the direct
  parent), "sector_id", "sector_name", "shell_index", "band_index", "wedge_index"}]}`, at most
  8, exact names first, then names that start with the term. A planet or
  moon also has `system_id` and `system_name`; its `ref` opens through
  `/object/<ref>`. Asteroid belts have no names, so none is found. Sectors
  with no address and systems outside a sector are left out, since the map
  can't fly to them. Each `/api/search` row likewise carries `ref` and
  `parents`.
- `GET /api/galaxy/stamp` — `{"stamp": "<16 hex characters>", "state":
  "<token>"}` (`queryDb.galaxy_content_stamp`). `stamp` changes whenever
  node contents could: sectors placed, edited or removed (their
  `modified_at`, schema v27), new star systems, a re-planned galaxy shape,
  or a new planetGen release. `state` is what `/api/galaxy/changes` takes.
  `tile_stamp` (in `/api/galaxy/changes`) is the same without the release, for the stored nodes (PERF.133).
- `GET /api/galaxy/changes?since=<state>` — `{"stamp", "tile_stamp", "state", "full",
  "tile_full", "stages", "nodes"}` (`queryDb.galaxy_changes`): the keys of the nodes
  that changed since that `state`, found from `sectors.modified_at`, new
  sector ids and new star-system ids (`nodes`, the keys `queryDb.galaxy_node_keys` gives for each changed
  sector: every block of the ladder over it, from the root down, and the sector's own), and the drill-down stages
  holding those sectors (`stages`, the keys `queryDb.galaxy_stage_keys` gives for each changed
  sector: `galaxy` and every block above the sector up to level 243). `tile_full` is `true` (and both lists empty) when that can't
  be pinned to nodes: a deleted sector, a new shape, node format or generator epoch, more than
  1,000 changed sectors, new bright stars (a backfill), or a missing or
  unreadable `since`; `full` is `true` for all of those and also for a new release alone, which stales the pages and counts cached by `stamp` but no stored node (PERF.133). When `tile_full` is only because so much changed at once
  (more than 1,000 sectors, or new bright stars), the answer also has
  `"busy": true`, and the web layer keeps its cached nodes for up to ten
  minutes instead of refetching them all under that load (PERF.34).
  `../src/planetgen/web/lib/tilecache.py` (the web layer's disk cache) calls it
  about once a minute and deletes only the listed nodes, and passes the
  list on to the map's browser cache.
- `GET /api/phenomena?limit=<n>&offset=<n>&after=<scatter id>|<type>:<ID>` — every exotic phenomenon (`after`, the `scatter_id` of a scattered
  row, starts the page after it in the tail of unbuilt scattered phenomena and ignores `offset`, PERF.75; `after=<type>:<ID>`, the `type` and `id` of a
  built row, starts it after that row in the order of `sort` and then the tail, PERF.205, PERF.211: the first page of every sort (`name`, `type`,
  `descriptor`, `radius`, `placed`, `sector`) and every page after a built row cost the same at any depth, each table read on an index ending in
  `name` and the best rows of all taken, a tie settling by name, type and id in the direction of the sort; a page by `offset` past the first still
  sorts the whole union of built rows; anything else is a 400;
  `from=<scatter id>`, PERF.135, lists that tail from the row with that id on, `total` counting from there; a value that is not a whole number from 1 up
  is a 400, and one past the end answers the last page of the tail, `past_end` true and `total` its length),
  across every sector and regardless of galaxy placement (unlike
  `/api/galaxy/phenomena`, which only returns the galaxy-placed subset) —
  paginated the same way `/api/sectors`/`/api/systems` are (`items`,
  `total`, `limit`, `offset`). Each item has `id`, `type`
  (`"nebula"`/`"asteroid_field"`/`"black_hole"`/`"neutron_star"`/
  `"supernova_remnant"`/`"rogue_planet"`/`"interstellar_comet"`/`"quasar"`), `name`,
  `descriptor`, `radius_ly` (same shape as `/api/galaxy/phenomena`'s own
  items), plus `sector_id`/`sector_name` (both `null` if never linked to a
  sector) and `placed` (bool, whether it has a galaxy position at all) —
  `queryDb.list_phenomena`. Excludes a black hole/neutron star that's
  actually anchored to a normal star system (`star_id` set) — that one's
  already shown on its own system's page, not a standalone phenomenon.
  The data the `/phenomena` page shows. The Phenomena table (UX.41) also
  takes `sort` (`name`, the default, `type`, `descriptor`, `radius`,
  `sector` or `placed`) with `order=asc|desc` (ties fall back to name, so
  pages never overlap), the filters `type` and `descriptor` (each repeatable,
  any of) and `placed=yes|no`, and `facets=1` to add `facets`:
  `{"type": [{"value", "count"}], "descriptor": [{"value", "count"}]}`, the
  options for the two filter menus. Each menu's counts apply every filter
  except its own, so picking a type narrows the descriptor options and the
  other way round. `total` counts only the rows that pass the filters; an
  unknown `sort`, `order` or `placed` is a 400.
- `GET /api/nebulae/<id>/shape[?lod=low|full]` — a nebula's shape as a
  triangle mesh (GEN.75): `{id, radius_ly, center_pc, lod, vertices,
  faces}`, the vertices `[x, y, z]` in units of the nebula's radius from
  its center (multiply by `radius_ly`), the faces triples of vertex
  indexes. `low` (the default) is coarse, for the Galaxy Map. 404 for an
  unknown nebula or level.
- `GET /api/nebulae/<id>/surroundings` — the brightest stars round a
  nebula, for its page's 3D view (MAP.105): `{radius_pc, half_width_pc,
  stars: [{x, y, z, luminosity_sol, temperature_k}]}`, positions in
  parsecs from the nebula's center, the most luminous first (at most 300,
  from a box 2.5 radii each side, at least 30 pc). No stars for a nebula
  never placed. 404 for an unknown nebula.
- `GET /api/uncharted-phenomena/<scatter id>` — one scattered phenomenon by its
  `phenomenon_scatter` row (MAP.211): `type`, `kind`, `subtype`, `built`, `address`
  (`[shell, band, wedge]`), `x`/`y`/`z` and `distance_from_core_pc` (parsecs), `mass_solar`,
  and for a row with a velocity `speed_kms`, `direction` and `flight_years`; a built row has
  `target` (`{"kind", "id", "name"}` of the object it became, `kind` `system` for a
  hypervelocity star) or `null` when none matches. It also has `member_count` and `system_mass_solar`
  (UX.114) and, for an unbuilt row with members, `members` and `nodes` as in `GET /api/uncharted-systems/<id>`.
  404 for an id that is no row (`queryDb.scattered_phenomenon_detail`). The rows of `GET /api/phenomena` that are
  scattered carry `scatter_id`, `member_count` and `system_mass_solar`. The `bright_stars` of
  `GET /api/galaxy/bright-stars` and the `stars` and `scattered` of `GET /api/galaxy/uncharted` carry the same two.
- `GET /api/phenomena/<type>/<id>/contents` — what sits inside a nebula or supernova remnant (UX.144): `counts` by kind
  (`system`, `nebula`, `asteroid_field`, `rogue_planet`, `interstellar_comet`, `black_hole`, `neutron_star`), `total`, and
  up to `limit` (default 100, at most 500) `items` (`{kind, id, name}`, by kind then name). A type that cannot contain
  anything, or an id that is no row, is a 404 (`db.cloudcontents`). A system's own compact remnant is not counted.
- `GET /api/containers[?kind=globular_cluster|superbubble]` (DB.85) — the stored containers, a globular cluster inside its tidal radius or a superbubble inside its cavity, by kind then id: `{"containers": [{id, uid, kind, name, component_id, center_pc, extent_pc, region: {shell, band, wedge}, bubble?}]}`. `uid` is the printed object ID (the container kind of the layout); `bubble` (`[column i, column j, draw]`) is set on a superbubble only. Another kind is a 400 (`db.containers`).
- `GET /api/containers/at?x=&y=&z=[&kind=]` — the containers holding a galaxy-frame point, parsecs, deepest first, each with `depth` (distance from the centre over the extent) and `distance_pc`. A missing or non-numeric coordinate is a 400.
- `GET /api/containers/<object ID>` and `GET /api/containers/<object ID>/contents[?limit=]` — one container (a globular cluster also carries `component`: its `mass_solar`, `age_gy`, `metallicity_fe_h`, `half_mass_radius_pc`, `tidal_radius_pc`, `king_w0` and `core_collapsed`), and what sits inside it (`counts` by kind, `total`, up to `limit` `items` whose `id` is the printed object ID, as the phenomenon contents above). An ID that is no container is a 404. A container is also an object reference (UX.184): `globular_cluster:<ID>` and `superbubble:<ID>` resolve at `GET /api/objects/<ref>` (the galaxy above it, its centre in `positions.galaxy_pc`, no keep-out radius), serve a sheet (`/sheet`: radius of its sphere and birth shell, band and wedge) and hold no children; a cluster's ID named as a superbubble (or the reverse) is a 404.
- `GET /api/phenomena/<type>/<id>` — one phenomenon's full detail (every
  column its own table has, e.g. a nebula's `composition`/
  `formation_cause`, a black hole's `mass_solar`/`spin`/
  `has_accretion_disk`, a supernova remnant's `morphology`/`progenitor_type`/
  `age_years`), plus `type`, `sector_name` and `nearest` (its three
  nearest star systems, `{id, name, distance_ly}`), `previous` and `next` (UX.185: `{id, name}` of the free phenomenon of the
  same type before and after it, or `null`) and, for a supernova remnant with a stored centre, `compact_remnant_name`. An asteroid field also
  has `composition`, its saved composition rows in order as
  `{component, concentration}` (`concentration` is `high`/`moderate`/
  `small`/`trace`), and an interstellar comet `composition`, its
  components in order — `queryDb.
  phenomenon_detail`. `type` is one of `nebula`/`asteroid_field`/
  `black_hole`/`neutron_star`/`supernova_remnant`/`rogue_planet`/
  `interstellar_comet`/`quasar`; an unrecognized type or
  a nonexistent id is
  a 404. The data the `/phenomenon/<type>/<id>` page shows — this
  project's first per-phenomenon info page (previously a phenomenon had no
  detail page of its own, only a hover tooltip on the Sector/Galaxy Map),
  now including a to-scale AU diagram (`../src/planetgen/web/maps/phenomenonmap.py`).
- `GET /api/search?sector_q=&system_q=&star_q=&planet_q=&moon_q=&<facet>=<value>...` —
  the faceted search behind the `/search` page: click-to-filter tags
  (object type; star spectral/luminosity class; planet/moon class, body
  type, supported life chemistry, the equipment a human needs
  (`equipment`/`moon_equipment`, one value per item the body needs, such as `pressure_suit` or `hardened_habitat`; a body matches when it needs any ticked item; GEN.655); asteroid belt density; phenomenon type
  and phenomenon class, the class tag written `<type>:<class>` such as
  `phenomenon_class=nebula:D` — repeat a facet
  name for multiple active values, e.g. `class=M&class=K`), a per-entity
  name search, and a min/max size range per entity (UX.214 also adds the score filters below) —
  `star_min_radius_km`/`star_max_radius_km` (likewise `planet_`/`moon_`),
  in km, either bound optional (a size range is a continuous quantity,
  not a discrete facet value, so it's its own pair of query parameters
  rather than a tag). `stars` (UX.239) keeps the systems with that many stars: a count (`stars=3`, exactly three), a count and a plus (`stars=3+`, three or
  more) or a name (`single`, `binary`, `trinary`, `quaternary`: exactly that many; `quintenary`: five or more); the Systems panel then runs with or without
  a name, every system row carries `star_count`, and anything else is a `400`. Returns `facets` (one `{value, label, count,
  tooltip}` list per facet, built from the distinct values actually
  present), `autocomplete` (`sectors`/`systems`/`stars`/`planets`/`moons`
  name lists), `facet_labels` (`"facet:value"` -> label, for an
  active-filter chip), and `results` (`sectors`/`systems`/`stars`/
  `planets`/`moons`/`belts`/`phenomena` -> `{"rows": [...], "total", "limit",
  "offset", "truncated"}`, or
  `null` for an object type with no active reason to query it — see
  `queryDb.search`'s docstring for the exact inclusion rule; a size range
  alone is reason enough, same as a tag or name term). `stars`/`planets`/
  `moons` result rows each include their own `radius_km`, and a planet or moon row also `equipment_tier`, `phi4` and the six
  `tier_<domain>` ratings (0 Blue to 3 Red; `null` for a body that was not scored), which the page draws as symbols (UX.100), and
  the four life scores as the cluster beside them (UX.152). Each result
  panel is paged on its own: `limit` sets the rows per panel (default
  300, clamped to 500 like "Pagination" below), and `sectors_offset`/
  `systems_offset`/`stars_offset`/`planets_offset`/`moons_offset`/
  `belts_offset`/`phenomena_offset` pick each panel's page; `total` is that panel's full
  match count and `truncated` is true when `rows` isn't all of them. An
  offset past the last match returns the last page (with its real
  `offset`). `panels` (comma-separated panel names, or `none`) runs just
  those panels, the others answering `null`; `none` still answers `facets` and `autocomplete`. The Search page asks for each
  panel on its own, so a panel that hits the statement limit (a `504`) times out alone. A name search of one word
  shorter than the database's full-text minimum (usually 3 letters) matches the start of a name.
  UX.214, score filters (planets by the plain name, moons with `moon_` in front, each repeatable, a body matches when it
  holds for any ticked value of a facet and all facets must match): `habitability=blue|green|yellow|red` (the overall PHI-6
  colour: Blue 1 or more, Green 0.75 to under 1, Yellow 0.4 to under 0.75, Red under 0.4); `domain=<domain>:<colour>`
  (domains `pressure`, `temperature`, `air`, `water`, `radiation`, `suitability`; `suitability:green`); `life_level=<kind>:<1|5>`
  (kind `microbial`, `flora`, `fauna` or `sentient`, that life score 1 or more or 5 or more, of 9); `living=alive|lifeless|food|civilization`
  (`food` is flora, fauna or a civilization); `supports=<class id>:<1|5>` (MAP.345: the body supports that life class, `1a`, `2`, `6`, with a class score of 1 or more or 5 or more of 9, read from the class score table's index; the class that actually lives on the body is the `life=<class id>` facet). Unscored bodies (gas giants, bodies saved before the scores) match none of them.
  `planet_sort`/`moon_sort` order a panel by `phi4`, `suitability`, a domain name, `microbial`, `flora`, `fauna` or `sentient`:
  worst first, `-` in front for best first (`-phi4`); an unscored body comes last either way. The facets answer their options with counts.
  UX.226, range filters and sorts on stored physical values (`src/planetgen/db/search_values.py`): for the star, planet and moon panels
  (`star_`, `planet_`, `moon_` in front), `<entity>_min_<key>` and `<entity>_max_<key>` keep the bodies whose value lies in the range
  (either bound optional, both inclusive; a body with no stored value is out of any range on that key), `<entity>_<key>_unit` names the
  unit the bounds are given in (any unit of the same kind in the unit registry, GEN.568; without it the key's own unit below), and
  `<entity>_sort=<key>` or `-<key>` orders the panel smallest first or largest first, a body with no value last. Several ranges all
  have to hold, and a range alone is reason enough to run its panel. The keys: `gravity` (g0, surface gravity, planets and moons),
  `surface_temperature` (K), `surface_pressure` (bar), `surface_dose` (mSv/yr, the surface radiation dose), `magnetic_moment` (A*m2),
  `ocean_depth` (km), `orbit_distance` (au for a planet, km for a moon), `day_length` (h, the solar day), `rotation_period` (h),
  `year_length` (yr), `mass` (Mearth for a planet or moon, Msun for a star), `moon_count` (a planet's number of moons), and for stars
  `luminosity` (Lsun), `age` (Gyr) and `temperature` (K); and `gas_<formula>` (such as `gas_N2`, `gas_CO2`, `gas_O2`, one per gas the
  generator knows) for the partial pressure of that gas in the surface air (bar; a body without the gas has 0). Every key reads the
  stored column, so a stored value moves the answer and nothing is computed from other columns. A row carries `values`,
  `{<key>: {"value", "unit", "symbol"}}` for each key the request used, in the unit the bounds were given in. The answer also gives
  `values` (the whole key table: label, the panels it applies to, the unit and symbol of each) and `active_ranges` (the ranges that were applied,
  with their units). A bound that is not a finite number, a minimum above the maximum or a unit that is not of the key's kind is a `400`.
  Example, planets of 0.8 to 1.2 g where the surface dose is under 2 mSv a year:
  `GET /api/search?planet_min_gravity=0.8&planet_max_gravity=1.2&planet_max_surface_dose=2`; the same pressure floor in kilopascals:
  `planet_min_surface_pressure=50&planet_surface_pressure_unit=kPa`.
- `GET /api/facilities/<id>` — one facility (see "Facilities" below).
- `GET /api/systems/<id>/facilities` — `{"items": [...]}`, every facility
  on a star, planet, moon or asteroid belt in that system.
- `GET /api/sectors/<id>/facilities` — `{"items": [...]}`, the
  stand-alone facilities parked in that sector and those on its asteroid
  fields.
- `GET /api/facilities/orbit?host_type=star|planet|moon&host_id=<id>[&distance_km=<km>]`
  — the orbit an orbital facility would get there, without saving:
  `{"distance_km", "period_years", "orbital_speed_kms", "min_distance_km",
  "max_distance_km"}`; the last two are the orbits the host allows, from
  just above its surface to the edge of its sphere of influence (a
  planet's or moon's Hill sphere, a star's heliosphere). Without
  `distance_km` it orbits at 3 host radii (or the sphere's edge, if
  nearer). `400` for a distance inside the host or outside its sphere,
  `404` for an unknown host. `POST /api/facilities` refuses the same
  distances; an asteroid facility in a belt takes no distance and gets a
  random spot in the belt, with the circular orbit around its star from
  there, which `planetgen.cli.orbits` advances like an orbital facility's.
- `GET /api/galaxy/bright-stars?shell=<i>&band=<j>&wedge=<k>[&all=1]` —
  `{"items": [...]}`: the bright stars the plan pre-placed in one sector
  cell, brightest first (`queryDb.bright_stars_in_sector`), only those not
  yet built into a system unless `all=1`. Each has `id`, `x`/`y`/`z` (pc),
  `luminosity_sol`, `temperature_k`, `star_type`, `yerkes_class`, the
  address and `system_id`. Empty when no scatter has run. `400` for a
  missing or non-integer address.
- `GET /api/population` — `{"generated", "species", "polities",
  "territories"}` booleans: whether a population pass has run and whether
  any species, polity or owned system exists (all `false` on a database
  without them). Pages that show population data hide themselves when
  the matching flag is false (`population.population_status`).
- `GET /api/species?spacefaring=true|false&limit=<n>&offset=<n>&after=<species id>` — the
  dominant species of every life world, by name, paginated (`after`, the `id` of the last species of the page before, starts the page after it in
  the same order and ignores `offset`, PERF.205; see
  docs/design/habitability/population-and-politics.md). Each item has its homeworld
  (`homeworld_planet_id`, `homeworld_name`, `star_system_id`,
  `system_name`), `life_class` (the id of its life class), `life_stage` (`multicellularity` or
  `technological_civilization`), `build`/`climate`/`size`,
  `civilization_age_years` and `era` (`null` without a civilization),
  `spacefaring`, and its `polity_id`/`polity_name` (`null` unless
  spacefaring). `400` for any other `spacefaring` value. The Species table
  (UX.41) also takes `sort` (`name`, `homeworld`, `system`, `era`,
  `spacefaring`, `polity`) with `order=asc|desc`, the repeatable filter
  `era`, and `facets=1` for `facets: {"spacefaring": [{"value": "yes"|"no",
  "count"}], "era": [...]}` (each menu ignores its own filter); `total`
  counts what passes the filters. Tech (POP.9): each item also has the six
  indices `energy_index`, `materials_index`, `information_index`,
  `medical_index`, `propulsion_index` and `defense_index` (whole numbers 0
  to 7), `tech_level` (0 to 7, two decimals), `tech_source` (`generated` or
  `admin`) and, derived, `tech_band` (the band's name), `limiting_domain`
  (the lowest of energy, materials and information) and `lopsided` (the
  strongest and weakest indices are 3 or more apart); all `null` until the
  species' tech is drawn (the next population pass). `sort=tech` orders by
  level, species with none last. `tech_min`/`tech_max` and
  `<domain>_min`/`<domain>_max` (`energy`, `materials`, `information`,
  `medical`, `propulsion`, `defense`) keep a range, ends included; a value
  that is not a number from 0 to 7 is `400`.
- `GET /api/species/<id>` — one species; `404` if unknown. It carries
  `limits` (POP.12): one row per environment limit of the human table
  (`name`, `unit`, `meaning`, `human`, `value`, `changed`), the species' own
  value in `value` and `changed` true where it differs from the human one.
  The limits are made from the species' homeworld until they are stored
  (DB.30).
- `GET /api/planets/<id>/habitable-species`, `GET /api/moons/<id>/habitable-species`
  — the species that could live on that body (POP.12), scored against each
  species' own limits through the six habitability domains: `unaided` (no
  domain is Red) and `equipped` (could with the kit its tech level can make,
  Water not Red, and not already unaided), each entry with `id`, `name`,
  `homeworld_planet_id`, `tech_level`, `phi4`, `equipment` (the required items its limits
  ask for, as the equipment line), `failing` (the Red domains) and `domains` (`{score, tier}` per
  domain), best first. `scored` is false for a gas giant or a body saved
  without a score (both lists empty). `404` if unknown.
- `GET /api/systems/<id>/habitable-species` — the same for every rocky
  planet and moon of a system that some species could live on: `bodies` of
  `{kind, id, name, unaided, equipped}`.
- `GET /api/species/<id>/habitable-worlds/systems/<system>` and
  `.../sectors/<sector>` — the other way round: the rocky planets and moons of
  one system, or of every system of one sector (at most 2,000 bodies; `truncated`
  says if there were more), that this species could live on, as `unaided` and
  `equipped` lists of `{kind, id, name, star_system_id}` and the same
  verdict fields.
- `GET /api/planets/<id>/species` — the species whose homeworld that
  planet is; `404` when it has none.
- `GET /api/polities?limit=<n>&offset=<n>` — every polity (one per
  spacefaring species), by name: `name`, `government`, `color`
  (`#rrggbb`), `reach_ly`, its species, `era`, capital and
  `system_count`. The Polities table (UX.41) also takes `sort` (`name`,
  `species`, `government`, `capital`, `systems`, `reach`) with
  `order=asc|desc`, the repeatable filters `government` and `era` (the
  species'), and `facets=1` for `facets: {"government": [...], "era":
  [...]}`; `total` counts what passes the filters.
- `GET /api/polities/<id>?limit=<n>&offset=<n>` — one polity plus a page
  of the systems it owns (`id`, `name`, `distance_ly` from the capital),
  nearest first, or by `sort=name|distance` with `order=asc|desc`; `404`
  if unknown.
- `GET /api/systems/<id>/owner` — `{"owner": {polity_id, polity_name,
  color, distance_ly}}`, or `{"owner": null}` when no polity holds it;
  `404` for an unknown system.
- `GET /api/territories` — `points`: up to 20,000 owned systems with
  galaxy-frame positions in parsecs (`x`, `y`, `z`) and their polity's
  `color`, nearest their capitals first; `polities`: each polity's `id`,
  `capital_pc` and `reach_ly`. For a 3D territory overlay.

### Write (admin auth required — see "Authentication" and "Write endpoints")

- `POST /api/sectors` — create a sector.
- `PATCH /api/sectors/<id>` — modify a sector.
- `DELETE /api/sectors/<id>` — remove a sector.
- `POST /api/sectors/<id>/generate-neighborhood` — generate every
  not-yet-generated sector within `radius_ly` (optional JSON body field,
  default 12 pc, at most the generator's own radius cap) of this
  galaxy-placed sector. The run is queued on Redis (see "Queued jobs"
  below) and the answer is `202` with `job_id` and `status_url`; its
  result is the counts `generated`, `already_existed`, `skipped`,
  `candidates` and `outside_galaxy`, plus the size and time `estimate`
  (see [`cli.md`](cli.md#size-and-time-estimates)).
  `"estimate_only": true` answers `200` at once with the counts and
  `estimate` and queues nothing. `"disregard_missing_data": true` (GEN.599)
  starts the run although a downloaded data file is missing, with each
  missing dataset's labelled fallback; without it the answer is `409` naming
  every file to download with its URL (see [`cli.md`](cli.md#datasets)). Before anything is queued: `404` for an
  unknown or unplaced sector, `409` when the galaxy has never been planned
  (`planetgen plan`) or a downloaded data file is missing, `507` when the database
  disk can't hold the run, `413` when the radius holds more than 2 million
  sectors (a web request lists no more). `503` when no Redis server answers. Every new
  sector also gets the bright stars (100 L_sun and up) within 100 ly of it
  (GEN.23). A large radius (100 ly) covers thousands of candidate slots,
  so this can run for a long time. The sector page asks for its
  `estimate`, then starts the real run as a Generate page job.
- `POST /api/galaxy/blocks/<key>/generate` — generate every
  not-yet-generated sector the galaxy's stored outline holds inside one
  block of the Galaxy Map's ladder (GEN.310; `<key>` is `m.shell.band.wedge`
  as in the map's stage links, `m` 243, 81, 27, 9 or 3), shell by shell
  outward, with the bright-star backfill and the settle after them. Queued
  and answered like generate-neighborhood: `202` with `job_id` and
  `status_url`, the same counts (`candidates` are the block's addresses
  inside the outline) and `estimate`; `"estimate_only": true` and
  `"disregard_missing_data": true` work the same. `400` for a key that names no
  block, `409` for a galaxy never planned or a missing data file, `413` for a
  block of more than 2 million sectors, `507` when the disk can't hold the
  run, `503` when no Redis answers. `planetgen sectors --block <key>` does the
  same from the command line.
- `POST /api/systems` — generate and create a system, standalone or in
  an existing sector.
- `POST /api/uncharted-systems/<bright_star_id>/generate` — admin only
  (UX.87): generate that one scattered star's system by itself (standalone;
  the sector's other contents stay uncharted) and link the star to it.
  Answers `201` with `{"id": <system id>}` (or `202` with a job id when the
  queue is slow); `404` for an unknown star, `409` if it already has a system.
- `PATCH /api/systems/<id>` — rename a system (its stars, planets and
  moons follow; see "Renaming" below), and/or regenerate its contents
  in place (see "Regenerating a system" below).
- `PATCH /api/stars/<id>` — rename a star.
- `PATCH /api/planets/<id>` — rename a planet.
- `PATCH /api/moons/<id>` — rename a moon.
- `DELETE /api/systems/<id>` — remove a system.
- `POST /api/facilities` — add a starbase, colony or outpost (see
  "Facilities" below).
- `DELETE /api/facilities/<id>` — remove a facility.
- `DELETE /api/planets/<id>`, `/api/moons/<id>`, `/api/belts/<id>` and
  `POST .../<id>/regenerate` — delete one body, or roll it again in
  place (see "Deleting and regenerating" below).
- `DELETE /api/phenomena/<type>/<id>` and `POST
  /api/phenomena/<type>/<id>/regenerate` — the same for one phenomenon.
- `DELETE /api/sectors/<id>/contents` — remove a sector with everything
  in it; `POST /api/sectors/<id>/regenerate` — remove it and generate its
  galaxy slot again.
- `POST /api/systems/<id>/star` — change a single-star system's star (see
  "Changing a star" below). A body's class is derived from its conditions
  (GEN.646), so there is no class change; edits of the inputs belong to
  the 8.6 edit chain.

### Admin stats (admin auth required)

Both take `?db=` like the read endpoints, and need an admin past the
forced credential change. They back the admin stats page
(`/admin/stats`, `../src/planetgen/web/admin_pages.py`).

- `GET /api/admin/stats` — health and statistics for one database:
  `api` (version, Python version, process uptime, load average, memory),
  `mysql` (server version, uptime, connected threads; `null` fields when
  `SHOW GLOBAL STATUS` isn't allowed), `redis` (whether the queue's Redis `answering` at its `url`, without any
  password, and, when it does, `version`, `uptime_seconds`, `aof_enabled`, `maxmemory_policy`,
  `used_memory_bytes`, `maxmemory_bytes`, `last_save_ok`, `last_save_seconds_ago` and the live `workers`; `error`
  says why when it does not), and `database` (`reachable`,
  `schema_version`/`schema_expected`/`schema_current`, `generator_version` (OPS.37, the same counter as in `/api/health`), `size_bytes`,
  exact `counts` for `sectors`/`star_systems`, `sector_stats` (PERF.11:
  `measured` and `backfilled` sector counts, and `ratio`, the decaying
  average of the systems fills got against the systems expected, over
  `fills` fills; `null` before v53), `tables` with
  `information_schema`'s estimated rows and data/index bytes,
  `timestamps` with each v27 table's newest `created_at` and latest
  `modified_at`, and `name_collisions`: how many base names had to be
  made unique per level plus `distinct_base_names`). An unreachable
  database still returns 200, with `database.reachable: false` and a
  `detail`.
- `GET /api/admin/duplicate-names` — paginated (`limit`/`offset`, by base
  name, alphabetical) list of the names the uniqueness rules decorated:
  `{"total", "limit", "offset", "items": [{"base_name", "levels",
  "rows"}]}`. `levels` names the registries it collided in (`sector`,
  `system`, `body`); each row is `{"kind", "id", "name"}` with
  `sector_id` for a system and `star_system_id`/`system_name` for a
  planet or moon.
- `GET /api/admin/login-failures` — the newest refused sign-ins (wrong
  username or password, a locked login, a wrong current password on
  change-credentials), newest first: `{"items": [{"action", "username",
  "ip", "created_at"}]}`, `action` being `login.failed`, `login.locked`
  or `password.failed`. `?limit=` (1-200, default 20). Takes no `?db=`:
  the rows live in the control schema's `admin_audit_log`, kept 90 days.
  Each one is also an `AUTH` line in the activity log
  ([`config.md`](config.md#the-activity-log)).
- `GET /api/admin/generation-stats` — this server's measured generation
  speed and each galaxy's size: `{"buckets": [{"kind", "bucket",
  "density_low", "density_high", "samples", "seconds_per_task",
  "seconds_per_system", "systems_per_task", "stars_per_system",
  "max_density"}], "stages": [...], "sizes": {database: {"bytes_per_system", "systems",
  "total_bytes"}}, "object_types": {kind: mean seconds per object}
  (PERF.92, `{}` before any run saved them), "available"}` (`available` is
  false until `update.sh` has added control schema v6). `stages` is the last 100 stage runs
  (PERF.56, control schema v14): `stage_key`, `stage_n`, `stage_total`,
  `label`, `skipped`, `skip_reason`, `seconds`, `workers`, `settings` (mass
  limit, luminosity floor, ...), `metrics` (layers visited and modified,
  objects) and `finished_at`. Shown on the Stats page.
- `POST /api/admin/generation-stats/reset` — deletes every recorded generation
  rate (PERF.32); the next run records afresh and estimates from defaults until
  it has. Returns `{"deleted": n}`; audited as `generation_stats.reset`.
- `GET /api/admin/naming-key[?db=]` — the galaxy's naming key (GEN.70):
  `{"database", "key", "codec_version", "current_codec_version",
  "drawn_at", "changed_at", "changed_by"}`; `key` is `null` before the
  galaxy is planned. `POST /api/admin/naming-key[?db=]` `{"key": "0A1B2C3D"}`
  (8 hex digits) or `{"draw": true}` sets or draws a key; every object the
  phoneme codec names takes a new name at once and no row is rewritten.
  409 before the galaxy is planned. Audited as `naming-key.change`. On the
  Stats page.
- `GET /api/admin/work[?limit=&offset=]` — the work queue (ADM.10, PERF.99):
  `{"available", "status": {"workers_active", "runs_active", "paused",
  "paused_by", "paused_at", "holder", "holder_title", "heartbeat_age_s",
  "stale", "load": {"values", "kind", "text"}}, "items": [runs, newest
  first], "total", "limit", "offset"}`. A run is `{"id", "title", "kind"
  ("generation", "maintenance" or "check"), "who", "database", "job" (the
  Generate page job that started it), "argv", "state" ("running",
  "cancelling", "done", "failed" or "cancelled"), "error", "live", "stale",
  "created", "started", "updated", "finished", "seconds", "stage": {"index",
  "count", "label", "state"}, "units", "units_done", "percent", "eta":
  {"seconds", "low", "high", "source" ("live", "recorded", "earlier runs", "this run so far" or "not timed yet"),
  "from_earlier_runs", "partial", "stalled", "rough" (PERF.184: part of the time is an estimate from the host's workers and
  the stage's units, so `low` and `high` are wide)} or null}`, all read from the run record in Redis (`holder` is
  the run that owns the machine; `stale` says it has not been heard from
  for `max(120 s, 3 x its longest unit)`). `load.kind` is `load` (the 1, 5
  and 15 minute load average); `load.values` is null where the system has
  no load average. `available` is false when no Redis server answers.
- `GET /api/admin/work/warm?database=NAME` — what the cache warm-up has done on that database, phase by phase (PERF.141):
  `{"available", "run": {"state", "started_at", "finished_at", "reason", "error"} or null, "phases": [{"name", "text",
  "families", "state", "state_text", "line", "built", "cached", "left_out", "seconds", "started_at", "finished_at"}],
  "auto": {"enabled", "wanted", "holder"}, "index": {"size_bytes", "kinds": [{"kind", "rows", "mode", "built_at",
  "seconds"}]}}`. `index` is the object index as its last runs recorded it (PERF.140). A phase's `state` is `waiting`, `running`, `built`, `unchanged`, `stopped`,
  `nothing`, `failed` or `interrupted` (null before any run); `auto` says whether a change starts the warm-up by itself,
  which changes it waits to warm for, and the process that holds it.
- `GET /api/admin/work/workers` and `POST /api/admin/work/workers`
  `{"count": n}` or `{"ludicrous": true}` — the bulk queue's worker count
  (ADM.15, PERF.95). `{"wanted", "live", "ludicrous", "cpus"}`: the count
  asked for (the admin's, else 80% of the cores), the workers running now,
  the count "Ludicrous Speed" sets (95% of the cores) and the cores. A
  POST needs a fresh admin login, is audited as `work.workers.set`, and
  lowering the count asks the extra workers to finish their unit and exit.
  503 when Redis does not answer.
- `GET /api/admin/work/<id>` — one run (PERF.99): `{"run": the run above with
  "stages": [{"key", "label", "state" ("waiting", "running" or "done"),
  "units", "done", "cost", "cost_done", "started", "ended", "seconds",
  "longest_unit"}], "events": the last 300 lines of its log [{"id", "t",
  "level", "kind", "message"} and, where the event has them, "stage" (from 0), "stages", "unit" and "seconds"], "log": true}`; 404 for an unknown run.
  `GET /api/admin/work/<id>/log` answers the whole log as `{"text": ...}`,
  from Redis while the run's record lasts (30 days) and from the file in the
  jobs directory's `runs` folder after (the last 20 are kept).
- `POST /api/admin/work/<id>/control` `{"action": "cancel"|"clear"}` —
  cancel stops a live run now: its waiting units are taken off the queue,
  its started units are sent RQ's stop command, and the machine is free
  (measured under a second for a run of 4,000 units); clear ends a run that
  has stopped answering (its process is gone, or it is stale) as failed and
  frees the machine, and is refused for a run that still answers.
  Answers `{"ok": bool, "cancelled": n, "stopped": n}`.
  `POST /api/admin/work/<id>/delete` deletes a finished run's record;
  `POST /api/admin/work/clear-finished` deletes every finished run and
  answers `{"cleared": count}`. `POST /api/admin/work/queue`
  `{"action": "pause"|"resume"}` pauses the whole queue (a flag in Redis,
  `planetgen:queue:paused`: the workers of a run's units and of the maintenance work take no job
  until it is cleared, a unit already running finishes, and the two API lanes are never held; PERF.152)
  or resumes it; `503` when Redis does not answer. Each answers
  `{"ok": bool}` and writes the change to the audit and activity logs
  (`work.cancel`, `work.clear`, `work.delete`, `work.queue.pause`, ...).
- `GET /api/admin/settings` — the configuration page's data (ADM.84):
  `{"file_hash", "problem", "file_error", "groups", "options", "reauth",
  "pending_restart"}`. Each option is `{"path", "title", "description",
  "category", "kind" ("bool", "int", "float", "text", "list" or "choice"),
  "choices", "limits", "unit", "default", "value" (what config.json holds,
  the default where it holds nothing), "effective" (the value in use),
  "is_set" (a secret only), "source" ("default", "config.json" or
  "environment"), "env", "restart", "secret", "editable", "min_role",
  "allowed" (for this admin), "reason", "fixed_kind" ("" when allowed, else "server", "owner" or "reauth")}`. An
  option of the Owner tier (the mail server, `base_url`, the From address) is allowed only to the Owner signed in
  with a session whose password was typed within the last 10 minutes (ADM.86; an API key is never more than an
  admin); `reauth` is true when that is all the caller lacks (`POST /api/auth/confirm-password`). A secret's `value`, `default` and
  `effective` are `null`. `problem` says why the file can't be saved (missing
  or not writable by the web user). `pending_restart` lists the options
  saved since this process started that need a restart.
  `POST /api/admin/settings/plan` `{"changes": {"jobs.keep": 9,
  "site_name": {"reset": true}}}` answers what a save would do and writes
  nothing: `{"ok", "changes": [{"path", "old", "new", "secret", "restart",
  "reset"}], "errors": {option: message}, "shadowed": {option: variable},
  "restart", "file_hash"}` (a secret's `old` and `new` are `(changed)`).
  `POST /api/admin/settings` `{"changes": {...}, "expect_hash": "..."}`
  checks and saves: the answer has the same shape with `saved` and `backup`
  (the archive of the old file) and `audit_errors`; a value that doesn't fit
  is a 200 with `saved` false and the message at its option, a file changed
  since `expect_hash` or one that can't be written is a 409. One
  `config.update` audit row per changed option. The Owner's and the command
  line's options are refused with their reason until accounts (8.2).
- `GET /api/admin/backups` — the backup archives (ADM.85): `{"problem",
  "archives": [{"name", "taken_at", "taken_by", "reason", "version", "size_bytes",
  "complete", "readable", "parts": [{"name", "label", "size_bytes", "files",
  "known"}]}], "kept", "parts": [{"name", "label"}], "can_download",
  "file_hash"}`. `POST /api/admin/backups` (`{}`) makes one now, of every
  part (`backup.take` audit row; 409 when the folder can't be written).
  `GET /api/admin/backups/<name>` is one archive with `verified`,
  `problems` and `diff`: what restoring its settings part would change
  (the shape of the settings plan; secrets as `(changed)`; `errors` when the
  backup no longer fits the model). `POST /api/admin/backups/<name>/restore`
  `{"parts": ["settings"] or null for all, "expect_hash": "..."}` restores
  after verifying, taking a safety archive first (`safety` in the answer),
  with a `config.restore` row per changed option and a `backup.restore` row;
  a part this build can't restore is in `skipped`. 409 for a damaged archive,
  a backup that doesn't fit, or a settings file changed since `expect_hash`.
  `GET /api/admin/backups/<name>/download` is the archive file and is the
  Owner's (403 for an admin until accounts, 8.2).
  `GET /api/admin/settings/export` — the options set in config.json that an
  admin may change, as the file's own JSON, without secrets or
  command-line-only options. `POST /api/admin/settings/import`
  `{"document": {...}, "expect_hash", "check_only"}` checks an export like an
  edit (the settings plan's answer; `request` holds the changes) and, unless
  `check_only`, saves it as `config.import`.
- `GET /api/admin/images` — the site's icons and link-preview picture
  (ADM.44): `{"images": [{"kind", "title", "help", "formats", "setting",
  "name", "url", "mime", "width", "height"}], "built_in_preview",
  "max_bytes"}`; `name` and `url` are `null` for the built-in one (or no
  icon). `POST /api/admin/images/<kind>` `{"data": "<base64>"}` (kind is
  `favicon_svg`, `favicon_png`, `apple_touch_icon` or `default_image`)
  checks the picture (format, size, square where an icon must be, a
  sandboxed SVG; 422 with the reason when it can't be used), stores it as
  `<16 hex digits of its SHA-256>.<png|ico|jpg|svg>` (served at
  `/site-assets/<name>`) and saves the setting that names it through the
  settings writer, so it is backed up and audited as `config.update`.
  `DELETE /api/admin/images/<kind>` resets the setting and removes the
  file when no other setting names it. At most 1 MiB a picture.
- `GET /api/admin/lockouts` — every address and username locked right
  now: `{"items": [{"scope", "subject", "retry_after", "locked_until",
  "level"}], "proxy_warning"}`. `proxy_warning` is true when the site
  looks like it is behind a reverse proxy without `proxy_fix` (a private
  address is locked, or most recent failures share one private or
  loopback address).
- `POST /api/admin/lockouts/lift` `{"scope": "ip"|"user", "subject"}`, or
  `{"all": true}` — lifts lockouts and forgets their counts; returns
  `{"lifted": n}`; audited as `lockout.lift`. From a shell (for an admin
  locked out of the site itself): `python3 -m planetgen.cli.lockouts` lists
  them, `python3 -m planetgen.cli.lift_lockout` with `--ip <address>`,
  `--user <name>` or `--all` lifts them, and
  `python3 -m planetgen.cli.forget_devices <name>` revokes that admin's
  trusted-device cookies.

### Authentication

- `POST /api/auth/login` `{"username", "password"}` — sets the session
  cookie (`__Host-pg_admin_session`, or `pg_admin_session` when
  `admin_cookie_insecure` is on; `HttpOnly`/`Secure`/`SameSite=Strict`;
  the `__Host-` prefix is SEC.32) and
  a trusted-device cookie (`pg_admin_device`, the same flags, 90 days;
  SEC.22). A login or change-credentials from a browser holding a valid
  device cookie for that username skips the per-username lock below (it
  is neither held back nor counted by it), so someone failing on purpose
  can't lock the real admin out; the per-address lock and rate limit
  still apply. Each successful login replaces the device cookie, and a
  credentials change revokes every device of that admin.
  Returns `{"username", "must_change_credentials"}`. `401` for a wrong
  username or password (same message either way — this never reveals
  whether a username exists; an unknown username costs the same password
  hash check as a wrong password, so timing doesn't reveal it either).
  Rate-limited to 10/minute/IP. Failed logins are also counted per
  address and per username, in Redis (shared by every worker, kept
  across restarts; SEC.1, SEC.21, SEC.30). Three
  failures from one address lock it for 5 minutes, and each further
  lockout of it doubles that, up to 1 day; a successful login clears its
  count, but its doubling level only halves per day without a lockout.
  An IPv6 address counts by its /64; loopback and `login_allowlist`
  addresses are never locked. Per username, whatever address they come
  from: after 10 failures, each further one locks that username for
  twice as long as the last (1 s, 2 s, 4 s, ... up to 15 minutes);
  unknown usernames are counted the same way, and a successful login or
  an hour without failures clears the count. A login while either is
  locked is a `429` with `{"error", "retry_after", "scope"}` (`scope`
  `ip` or `user`) and a `Retry-After` header, before the password is
  checked. With `ratelimit.storage_uri` set to `memory://`, or while Redis
  can't be reached, the counts are kept in memory per worker process (with
  one warning in the error log).
- Two-factor sign-in (SEC.26, optional per admin): for an admin who
  has turned it on, a right password at `POST /api/auth/login` answers
  `{"totp_required": true, "pending"}` and sets no session yet. Then
  `POST /api/auth/login/totp` `{"pending", "code"}` (within 5 minutes;
  `code` is a current authenticator code or an unused recovery code)
  sets the session like a login. A wrong code is a `401` and counts as a
  failed login (`totp.failed`) for the address and username, exactly as
  a wrong password does; a stale or forged `pending`, or one made before
  a password change, is a `401` asking to sign in again. API keys never
  need a code.
- `POST /api/auth/confirm-password` `{"current_password"}` — session only; types the password again, guarded like a
  login, and keeps the Owner-tier settings open for `{"minutes": 10}`.
- `GET /api/auth/totp` — `{"enabled", "recovery_codes_left"}`.
- `POST /api/auth/totp/setup` `{"current_password"}` — `{"secret", "uri",
  "qr_svg"}` for an authenticator app; nothing changes at sign-in yet.
- `POST /api/auth/totp/confirm` `{"code"}` — turns it on once a code
  from the app matches; returns `{"recovery_codes": [...]}`, shown only
  this once (10 codes, each good for one sign-in).
- `POST /api/auth/totp/disable` `{"current_password", "code"}` — turns
  it off, and forgets every trusted device of that admin (the browser
  that turned it off gets a new device cookie). Lost the phone and the recovery codes? From a shell:
  `python3 -m planetgen.cli.reset_two_factor <user>`.
- `POST /api/auth/logout` — ends the current session, clears the cookie.
- `GET /api/auth/me` — the calling admin's identity: `{"username", "must_change_credentials", "avatar"}`; `avatar` is
  `/site-assets/<name>` or `null`.
- `POST /api/auth/avatar` `{"data": "<base64>"}` — sets the account's picture (SEC.36); a signed-in browser session
  and past the forced credential change. **The avatar rules.** A PNG, JPEG, ICO or GIF (animated too), at most 1 MiB,
  1,048,576 pixels a frame (1024 by 1024), 100 frames, and 16,000,000 pixels over all the frames (width times height
  times frames). Pillow is asked to try those four formats only, the size is read from the header before anything is
  decoded, a decompression-bomb warning is an error, and the frames are unpacked one at a time under the caps, so an
  upload over a cap is refused (422 with a sentence) having unpacked no more than fits. A GIF or PNG with several frames
  keeps its frames, durations and loop count; everything else the upload carried (comments, XMP, ICC profiles, EXIF,
  other application extensions) is dropped, because the stored file is Pillow's new encoding of the pixels, never the
  uploaded bytes (an ICO is stored as a PNG, a JPEG keeps its turn). It is named by the SHA-256 of what is stored and
  served at `/site-assets/<name>` with its own content type, `X-Content-Type-Options: nosniff` and a year's cache. The
  answer is `{"avatar", "width", "height", "frames"}`. The site's logo, icons and preview picture keep the SEC.35 rule
  and still refuse a GIF. Replacing the picture removes the old file unless another account or setting names it.
- `DELETE /api/auth/avatar` — the account has no picture; its file is removed unless another account holds it.
- `POST /api/auth/change-credentials`
  `{"current_password", "new_username", "new_password"}` — always
  requires the current password, regardless of whether
  `must_change_credentials` is set. `new_password` must be at least 12
  characters, not equal to the username, not on the bundled list of
  common and breached passwords, and not just the username, "planetgen",
  "password" or "admin" with fewer than 8 other characters (a `400`
  saying which). A wrong `current_password` counts as a failed login for
  the address and the username, is a `429` while either is locked, and
  the route is limited to 10/minute/IP like login. On success every one of that
  admin's sessions ends (other browsers are logged out) and the caller
  gets a fresh session cookie; API keys keep working.
- `GET /api/auth/api-keys` — the calling admin's own API keys (label,
  `key_prefix` (its first eight characters), `scopes`, timestamps and
  `expires_at`; never the key or its hash).
- `POST /api/auth/api-keys` `{"label", "scopes"?, "expires_days"?}` —
  creates a key, returning `{"id", "label", "scopes", "expires_days",
  "key"}`. `key` is the raw value, shown exactly once. `scopes` are from
  `read`, `generate`, `upload` and `admin` and default to `["admin"]`;
  `expires_days` (more than 0, at most 3650) defaults to never, and an
  expired key answers `401` like a revoked one.
- `DELETE /api/auth/api-keys/<id>` — revokes one of the calling admin's
  own keys.

Authenticate either with the session cookie (set by `/api/auth/login`,
used by the `../src/html/` admin pages) or an API key, sent as
`Authorization: Bearer <key>`, for programmatic callers.
A key holds scopes (API.9): `admin` does everything below except the
account routes, `generate` and `upload` each include `read`, and `read`
reads who the key is (`/api/auth/me`) and a queued job's state
(`/api/jobs/<id>`). A route the key's scopes don't reach answers `403
{"error": "this key lacks the 'admin' scope", "required_scope":
"admin"}`. Every route that changes anything needs `admin` until the
upload and generate routes (API.18, API.12) exist. Each key has a rate
limit bucket of its own, and the per-address defaults do not apply to a
valid key.

An API key can do anything a signed-in admin can except manage the
account: making API keys, changing credentials, setting up, confirming
or turning off two-factor sign-in and logging out answer `403` to a key
and need the session cookie, so a leaked key can't make itself a
replacement before it's revoked. A key can still revoke a key
(including itself) and read `/api/auth/me`, its keys and the two-factor
status. When a request carries both, the key is what counts.

### Pagination

`/api/sectors`, `/api/systems` and `/api/phenomena` return a paginated envelope rather than a
bare list — the galaxy is generated lazily at
galaxy scale (the original roadmap's phase 4, now done), so an unbounded listing endpoint would eventually
return an unbounded response:

```json
{
  "items": [ ... ],
  "total": 1234,
  "limit": 100,
  "offset": 0
}
```

`limit` defaults to 100 and is silently clamped to 500 (asking for "too
much" isn't an error, just more than one response will return); `offset`
defaults to 0. `limit` must be a positive integer (≥1) when given at all;
`offset` must be a non-negative integer (≥0) when given at all. Either
violation gets a `400`.

### Errors

A request body that fails validation (ADM.21: every body is a Pydantic
model in `planetgen/api/schemas.py`) is a `400` whose `{"error": "..."}`
message joins every problem and which also carries `"errors"`, the list of
those messages, so a client can show all of them at once.

Every error response is JSON, `{"error": "..."}`, regardless of what raised
it: a missing sector/system id is a `404`, a missing/invalid query parameter
or request body (including an out-of-range `limit`/`offset`, a
non-numeric/non-positive `radius`, or a write endpoint's body failing
validation — see below) is a `400`, an unmatched URL is a `404`, an
unsupported HTTP method is a `405`, a request body over 2 MB
(`MAX_CONTENT_LENGTH`) is a `413` (`{"error": "request body too large"}`;
a page answers its own HTML error page), exceeding a rate limit is a `429`
(`{"error": "rate limit exceeded", "detail": "..."}`, see "Rate limiting"),
a database that can't be opened (the content database, or the control
database an admin route checks credentials against) is a `503` whose
message is just `"database unavailable"`, a row a long job holds locked past
`WEB_LOCK_WAIT_SECONDS` is a `503` `{"error": "BUSY: ..."}` and a read that
ran past `mysql.statement_timeout_seconds` is a `504`
`{"error": "QUERY_TIMEOUT: ..."}` (both with a `Retry-After` header; a page
answers its own HTML page), and an unexpected server-side failure is a
`500` — the API never falls through to Flask's default HTML error page or
leaks a stack trace, a connection error or the database's user/host to
the client (the real detail still reaches Flask's own logger).

Every response also carries `X-Content-Type-Options: nosniff`,
`X-Frame-Options: DENY`, `Referrer-Policy: no-referrer` and
`Content-Security-Policy: default-src 'none'`, and, when the request came
in over HTTPS, `Strict-Transport-Security: max-age=31536000`.

## NAV

`GET /api/nav?from=<ref>&to=<ref>` returns a direct course
(distance, bearing and mark, warp and fold travel times) plus an optimal route via
adjacent systems (`planetgen.galaxy.nav_graph`, a k-nearest-neighbor adjacency
graph with Dijkstra shortest-path) between two endpoints. An endpoint is an
object reference (`GET /api/objects/<ref>`): a star system (`system:12`, or a
bare number), a body in one (`star:7`, `planet:40`, `moon:41`, `belt:5`,
`comet:8`), or a standalone phenomenon (`nebula:2`, `asteroid_field`,
`black_hole`, `neutron_star`, `supernova_remnant`, `rogue_planet`,
`interstellar_comet`, `quasar`). A sector is not an endpoint (`400`); a bad
reference or nonexistent id is a `404`. A phenomenon endpoint's own
`route.path`/`route.positions` id is a `"phenomenon:<type>:<id>"` string (a
plain int id for a system) -- only `route.path[0]`/`route.path[-1]` can ever
be a phenomenon; every intermediate hop is a system or a stop in unknown
space (below).

**Hops and unknown space (NAV.12).** There is no hop limit: a route always
finds the nearest star it can, across sector boundaries, and a pair in one
sector may route through the sectors next door when the sector is galaxy-placed
(the positions stay in the sector's frame). `route.hops` lists one
`{from, to, distance_ly, unknown_space}` per hop, and `route.longest_hop_ly`
is the longest. `unknown_space` is true when the hop's straight line crosses a
sector that has not been generated (`galaxyGeometry.sectors_along_segment`,
NAV.38). A route in a sector with no galaxy placement has the flag false.

**Stops in unknown space (NAV.47).** A route goes from system to system; only
where a hop crosses unknown space does it look for other stops, in the
uncharted cells beside the jump: scattered bright stars, black holes, neutron
stars, quasars and hypervelocity stars (not planetary nebulae or supernova
remnants), counting only rows nothing has been generated for yet. A chain of
them replaces the jump when it is no longer than the jump times 1.5
(`tuning.NAV_UNKNOWN_JUMP_PENALTY`); otherwise the single jump stays. Such a
stop's id in `route.path`, `route.positions` and `route.hops` is a
`"stone:<kind>:<row>"` string, its `route.stop_places` entry also has `kind`
(`bright-star`, `black-hole`, `neutron-star`, `quasar` or `hypervelocity-star`)
and `name`, every hop to or from one has `unknown_space` true, and a stay
(`stay`) is not spent at it. A hypervelocity star is where it is at the
database's orbit epoch. `GET /api/nav/chart` still reckons the direct jump.

**Courses along the route (NAV.42).** Each `route.hops` entry also has the
course from its stop to the next: `bearing_deg`, `mark_deg`, `elevation_deg`
and `frame`, worked out in the Sector Local Frame when both stops are systems
of one sector and in the Galactic Frame otherwise
(`docs/design/orbits/navigation-frames.md`). `route.stop_places` has one
`{node, sector_id, sector_position_ly}` per id in `path` (both `null` for a
phenomenon end). The NAV page shows each course as `045 mark 012, 3.2 ly`.

**Charting a course (NAV.48).** `GET /api/nav/chart?from=<ref>&to=<ref>[&border=1]`
works out the uncharted sectors that block a course, without generating
anything: `unknown_hops` (how many hops of the route cross unknown space),
`cells` (`[shell, band, wedge]` of the cells those hops cross that hold no
sector and lie inside the galaxy's outline, in the order the route enters
them; `border=1` adds the uncharted cells that share a face with them),
`count`, `outside_galaxy` (uncharted cells past the outline, never offered),
`route_distance_ly`, `bypass` (`{found, distance_ly, checked}`, `null` when no
hop is unknown: whether a route through charted space exists; `checked` false
means the search gave up after `tuning.NAV_BYPASS_MAX_EDGES` hops) and
`confirm_over` (`tuning.NAV_CHART_CONFIRM_SECTORS`, 5,000: past it charting
asks for a confirmation). The cells are those of the unknown hops only, not
every cell on the straight line. Errors are as for `/api/nav`. The NAV chart
page (`/nav/chart`, admins) starts the job through the Generate page, which
shows the size and time first and refuses a run the disk cannot hold.

**Travel times for the route (NAV.11).** Each `route.hops` entry also has
`warp_times` and `fold_times` for that hop alone (the same factors as the
top-level lists), and the route has `warp_times` and `fold_times` for the
whole trip: every hop timed from rest to rest at the factor's constant speed
(`warp_speed_c`, `fold_speed_c`; no acceleration model) plus a stay at each
stop between the two ends, a route of N hops having `route.stops` = N - 1 of
them. The optional `stay` query parameter is the minutes spent at each stop
(default 0, at most 10,000,000; anything else is a `400`), echoed as
`route.stay_minutes`. A one-hop route has no stop, so its total is the hop's.

**Bodies and legs (NAV.16).** `origin` and `destination` echo
`{ref, kind, name}`. `legs` lists the course's legs, each `{kind, from, to,
direct, warp_times, fold_times}`: `between` (the system-to-system or
phenomenon course, the same as the top-level `direct`), plus `out` when the
origin is a body (from it to its system's heliopause on the side facing the
destination) and `into` when the destination is one (from the heliopause on
the side facing the origin to the body), both in the System Local Frame
(`frame: "system"`). Two objects in one system (a system itself counts) are a
single `within` leg with `scope: "system"`, `route: null` and positions in
light-years from the system's origin. `note` says the in-system legs use the
same warp and fold tables as between systems (`null` when there is none).
A belt is placed at its mean radius.

**Availability.** NAV only applies to a pair of endpoints that satisfy one of:

- Both are systems in the SAME sector (`star_systems.sector_id IS NOT
  NULL`, matching) -- scoped to that sector's own systems (sector-local
  positions). Checked first/preferred over the galaxy-scope rule below.
- Otherwise, both have a galaxy-frame position -- a system whose sector has
  a galaxy placement (`sectors.center_x/y/z_pc IS NOT NULL`, i.e. placed by
  `planetgen galaxy`), or a phenomenon that's itself been placed in the galaxy
  (`--sector-id` at generation time -- see `planetgen phenomenon`'s own
  help). A phenomenon endpoint can only ever be resolved this way: it has
  no sector-local position of its own, regardless of which sector its own
  `sector_id` names as its nearest neighbor (that column is a browsing
  convenience, not real containment).
- Anything else (an unplaced system, a phenomenon never placed in the
  galaxy, or two systems in different sectors where either lacks a galaxy
  placement) is a `404` (unknown id) or a `400` with `"NAV requires both
  endpoints to share a sector, or both to have a galaxy placement"`.

**Course convention.** Courses use Boss's nested reference frames
(`docs/design/orbits/navigation-frames.md`, `planetgen/galaxy/navigation.py`):
`bearing_deg` is 0-360 with 0 pointing from `from` toward the frame's
center (flattened onto the galactic plane) and 90 to the East (North x Up);
`elevation_deg` is -90 to +90 above or below the galactic plane, and
`mark_deg` is that elevation mod 360 (0-90 up, 270-360 down). `frame` is
`"sector"` for `scope: "sector"` (center: the sector's center) or
`"galactic"` for `scope: "galaxy"` (center: the galactic core). The NAV
page shows this as "000 mark 000", rounded to whole degrees.

**Response:**

```json
{
  "scope": "sector",
  "direct": {
    "distance_ly": 4.0,
    "bearing_deg": 0.0,
    "mark_deg": 0.0,
    "elevation_deg": 0.0,
    "frame": "sector"
  },
  "warp_times": [
    {"warp_factor": 1, "velocity_multiple_of_c": 1.0, "years": 4.0, "formatted": "4 years"},
    {"warp_factor": 2, "velocity_multiple_of_c": 10.079, "years": 0.397, "formatted": "144 days 22 hours and 47 minutes"},
    "... warp 4, 8, 9, 9.5 and 9.9 ...",
    {"warp_factor": 9.995, "velocity_multiple_of_c": 12201.937, "years": 0.0003, "formatted": "2 hours and 52 minutes"}
  ],
  "fold_times": [
    {"fold_factor": 4, "velocity_multiple_of_c": 256.0, "years": 0.016, "formatted": "5 days 16 hours and 58 minutes"},
    "... fold 5, 6, 6.5, 7, 7.5 and 8 ...",
    {"fold_factor": 8.5, "velocity_multiple_of_c": 20880.25, "years": 0.0002, "formatted": "1 hour and 41 minutes"}
  ],
  "origin_position": [0.0, 0.0, 0.0],
  "destination_position": [4.0, 0.0, 0.0],
  "route": {
    "path": [1, 3],
    "distance_ly": 4.0,
    "positions": {"1": [0.0, 0.0, 0.0], "3": [4.0, 0.0, 0.0]}
  }
}
```

- `scope`: `"sector"` (same-sector, sector-local positions) or `"galaxy"`
  (cross-sector, absolute galaxy-frame positions).
- `direct`: straight-line course from `from` to `to`.
- `warp_times`: travel time for `direct.distance_ly` at warp 1, 2, 4, 8, 9,
  9.5, 9.9 and 9.995 (`program_constants.WARP_FACTORS_FOR_NAV`), on the warp
  curve `w^(10/3) + 1 / (1 + e^(-9.3575 (w - 9.5))) * (198.9 / (10 - w)^0.75
  + 1721.7 - w^(10/3))` in multiples of c (`navigation.warp_speed_c`),
  formatted via the same duration formatter used elsewhere in this project.
- `fold_times`: the same for dimensional fold 4, 5, 6, 6.5, 7, 7.5, 8 and
  8.5 (`FOLD_FACTORS_FOR_NAV`), at `6 F^4 / (10 - F)` times c
  (`navigation.fold_speed_c`).
- `origin_position`/`destination_position`: the `[x, y, z]` light-year
  positions `direct` was computed from, in `scope`'s frame (sector-local for
  `"sector"`, absolute galaxy-frame for `"galaxy"`) — what the NAV page's (`/nav`)
  NAV Map plot (`planetgen/web/maps/navmap.py`) draws.
- `route`: the shortest path via adjacent systems (nodes: every system in
  scope, plus a phenomenon endpoint's own one-off node when `from`/`to` is
  one; edges: each node's `k`-nearest neighbors, symmetrized, then each
  island that leaves linked to its nearest few islands by their closest
  pair of systems, until one is left -- NAV.34), as a list of
  node ids from `from` to `to` inclusive (a plain int for a system, a
  `"phenomenon:<type>:<id>"` string for a phenomenon -- only `path[0]`/
  `path[-1]` can ever be the latter), its total distance, and `positions`
  (one `[x, y, z]` entry per id in `path`, keyed as a string either way
  since JSON object keys always are, same frame as `origin_position`/
  `destination_position`). `null` only when `from` and `to` are the same
  node: with the islands joined, any two placed endpoints have a route,
  however far apart the generated areas around them are.

## Write endpoints

Every write route requires an authenticated admin (session cookie or
`Authorization: Bearer <api-key>` — see "Authentication" above) whose
`must_change_credentials` flag is clear — the seeded `admin` bootstrap
admin (random first password, printed once by `planetgen.cli.migrate`) can log in and call `/api/auth/change-credentials`, but
nothing else, until it changes its own credentials (`403` otherwise). A
missing/invalid credential is `401`. Every write route runs against the
same `PLANETGEN_MYSQL_*` database account every read route above uses —
see `config.py` and [`deployment/README.md`](../deployment/README.md#mysql-accounts) — and
records one row in the control schema's `admin_audit_log` (who, what,
when) after the write actually succeeds.

### Sectors — request body

`POST /api/sectors` (all fields required) takes:

```json
{
  "name": "Voranthis Kelmoor",
  "edge_ly": 11.5
}
```

`PATCH /api/sectors/<id>` takes any non-empty subset of the same two fields:

```json
{
  "name": "Voranthis Kelmoor",
  "edge_ly": 11.5
}
```

- `name`: non-empty string.
- `edge_ly`: number, greater than 0 — the sector's cube edge, in
  light-years (matches `sectors.edge_mpc` after unit conversion; see
  `database-schema.md`).
An unrecognized field, a missing required field (`POST` only), a wrong
type, or a value failing the constraints above is a `400`.

### Systems — request body

`POST /api/systems` takes a generation **"recipe"**, the same shape
`planetgen system --system-file` already takes (`SystemConfig.to_dict()`/
`from_dict()`) — every field is optional (`null`/omitted means "let the
generator decide"), and the server generates a brand-new, real system from
it the same way `planetgen system` does, via `StarSystem(system_config=...)`:

```json
{
  "star_type": "G2V",
  "name": "Voranthis Vesta",
  "age": "old",
  "habitable_world": true,
  "asteroid_belt": null,
  "large_star": false,
  "moons": null,
  "max_planets": null,
  "planets": null,
  "intelligent_life": null,
  "binary_system": false,
  "wide_binary": null,
  "num_orbits": 6
}
```

Accepted fields: `star_type`, `name`, `age` (`"young"`,
`"old"`, or `null`), `num_orbits` (a positive integer up to 500, or `null`), and the
tri-state booleans `habitable_world`/`asteroid_belt`/`large_star`/`moons`/
`max_planets`/`planets`/`intelligent_life`/`binary_system`/`wide_binary`
(`true`, `false`, or `null`) — `wide_binary` selects an S-type (wide) vs.
P-type (close) binary and is only meaningful together with
`binary_system: true`. An unrecognized field (including `slots`, or a
fully-specified object graph shaped like `StarSystem.to_dict()`) is a
`400` — not silently ignored.

**Placing it in a sector.** Two more optional fields: `sector_id` (a
sector's id) and `position` (`[x, y, z]` light-years from that sector's
center; needs `sector_id`). Without `sector_id` the system is standalone
(`sector_id = NULL`, like `planetgen system`). With it, the system is
placed clear of every system already in the sector's Hill sphere
(`SpaceSector.add_system`, via `_db.add_system_to_sector`), or exactly at
`position`, which must lie inside the sector. A sector in the galaxy also
gives the system its stellar population (and, once bright stars have been
scattered, keeps it below the bright-star threshold), and its
containment and nearest-systems rows are filled in like any generated
sector's. The response then also carries `sector_id` and `position`. An
unknown sector is a `404`, a position outside it a `400`, and a sector
too full to place one more a `409`.

A generation failure (an internally-inconsistent recipe, e.g. an
impossible `num_orbits`/class combination) is reported as a `400`, not a
`500` — the request body was the problem, not the server.

### Regenerating a system

`PATCH /api/systems/<id>` with `{"regenerate": recipe}` (the recipe
fields above except `name`; `{}` rolls a fresh system with defaults)
replaces the system's stars, planets, moons, asteroid belts and comets
with a newly generated set. The system keeps its id, name, sector,
position and location, and the new bodies are named from the
system's name. A system in a galaxy sector keeps that sector's stellar
population. It is a `409` when the system was built around a pre-placed
bright star (its star is fixed by the galaxy scatter), or when facilities
are hosted on it, unless `"drop_facilities": true` is sent too (they are
deleted with the bodies). `name` may be sent in the same request (it is
applied first). Success returns `{"status": "ok", "id", "name",
"regenerated"}`.

### Deleting and regenerating

`src/planetgen/api/edits.py` (ADM.8). Every one of these needs an admin whose
credentials are current and writes an audit-log row.

- **Planets, moons and asteroid belts.** `DELETE /api/planets/<id>`
  (with its moons), `/api/moons/<id>` and `/api/belts/<id>` remove one
  body; `POST .../<id>/regenerate` rolls it again at the same orbit: a
  planet gets a new class (any that fits its zone), size, atmosphere,
  life and moons; a moon a new moon class that fits its planet; a belt a
  new density and composition over the same span. A regenerated body
  keeps its name and its row (so facilities on it stay). The rest of the
  system is then re-validated from the moons outward
  (`planetgen/generation/validation.py`): bodies are moved outward until the
  orbits are stable, never removed. The answer is `{"status": "ok",
  "summary", "moved", "reclassified", "removed", "warnings",
  "star_system_id"}`: the names of the bodies that moved or changed class
  to fit, and a plain sentence for anything still unstable. It is a `409`
  when facilities would be deleted (on a deleted body or its moons, or on
  a regenerated planet's old moons) unless the optional body
  `{"drop_facilities": true}` is sent.
- **Phenomena.** `<type>` is the phenomenon page's type (`nebula`,
  `black_hole`, `asteroid_field`, ...). Regenerating rolls the same type
  again, keeping its id, name, sector, galaxy position and what it sits
  inside. A black hole or neutron star that is a star system's star is a
  `409`: delete or regenerate the system instead.
- **Sectors.** `DELETE /api/sectors/<id>/contents` deletes the sector
  with its star systems, the phenomena filed under it and facilities
  parked in it; its slot in the galaxy is left empty, so it can be
  generated again. (`DELETE /api/sectors/<id>` still deletes only the
  sector row and keeps its systems as standalone ones.) Answers
  `{"systems", "phenomena"}` deleted. `POST /api/sectors/<id>/regenerate`
  does the same, then generates the slot again from the galaxy's density
  plan. It is queued (see "Queued jobs"): `202` with a job id, and the
  job's `result` has `deleted` and the new sector's id and name
  (`sector_id`, `sector_name`; `null` when the slot is outside the
  galaxy's outline). It is a `404` for an unknown sector and a `409` for a
  sector off the galaxy grid or before `planetgen plan`, both before
  anything is queued or deleted.

### Changing a star

Also `src/planetgen/api/edits.py`, with the same admin, audit and answer
shape as above.

- **Star (ADM.7).** `POST /api/systems/<id>/star` `{"star_type": "K2V"}`
  replaces a single star with a new one of that spectral type (class
  letter, subclass digit, Yerkes class; case doesn't matter), keeping
  its name, row and place in the galaxy. Planets, moons and belts keep
  their conditions; their orbits scale with the square root of the change
  in luminosity (so each keeps its place against the habitable zone),
  the innermost is moved clear of a larger star's surface, and then
  everything is re-spaced for the new mass. Bodies left past the new
  star's farthest stable orbit, and moons a planet that moved inward can
  no longer hold, are removed and listed in `removed`; comets' orbits
  scale the same way. The answer adds `star_type`. It is a `400` for a
  malformed type and a `409` for a binary system, a black hole or
  neutron star, a system built around a pre-placed bright star, or when
  removed bodies host facilities and `"drop_facilities": true` isn't
  sent.

### Renaming

`PATCH /api/stars/<id>`, `/api/planets/<id>` and `/api/moons/<id>` each
accept only `{"name": str}` (`PATCH /api/systems/<id>` also takes
`regenerate`, above). Runs of whitespace
collapse to one space; a blank name, one over 255 characters, or any other
field is a `400`, and an unknown id is a `404`. A name any other sector,
system or star already has is a `409`
(`{"error": "a star is already named 'Sirius'"}`). Planet and moon names
aren't checked: they come from their star's name, so only uniquely named
objects are searched. Success returns
`{"status": "ok", "id", "name"}` (plus `star_system_id` for a star,
planet or moon).

Generated names derive from the system (`src/planetgen/names/bodies.py`):
a single star shares the system's name (`Voranthis`), a close pair's stars
are its A and B (`Voranthis A`, `Voranthis B`), a wide pair's stars are the
system name's first word plus their own word (`Voranthis Kelmoor`,
`Voranthis Pikkita`), planets are numbered in orbit order (`Voranthis II`;
around a wide pair's secondary, after its own word: `Pikkita II`; planet
names may repeat across systems), and moons add a
letter (`Voranthis IIa`). Renames keep that in step:

- **System:** every star, planet and moon whose name starts with the old
  system name is renamed with it (for a wide pair, the old first word
  becomes the new first word). Names set by hand are left alone.
- **Star:** a single star shares its system's name, so this renames the
  system (as above). A binary's star is renamed on its own, along with the
  planets and moons named after it (a wide pair's primary's planets follow
  its first word, its secondary's follow its last).
- **Planet or moon:** just that body. A renamed planet's moons keep their
  names.

### Facilities

`POST /api/facilities` takes `name`, `kind`, `placement`, `host_type` and
`host_id` (all required), plus `distance_km` and `phase_deg` for an orbital
facility, `offset_ly` (`[x, y, z]` light-years from the sector's center,
along its own axes; the center if left out) for a stand-alone one, and an
optional `description`. It returns `{"id"}` with a `201`.

| `placement` | `host_type` | Kinds |
|---|---|---|
| `terrestrial` | `planet`, `moon` (terrestrial only) | `colony`, `outpost` |
| `orbital` | `star`, `planet`, `moon` | `outpost`, `station`, `starbase` |
| `asteroid` | `asteroid_belt` | `outpost`, `mining-colony` |
| `asteroid` | `asteroid_field` | `outpost` |
| `standalone` | `space` (`host_id` is a sector) | `outpost`, `station`, `starbase` |

Anything else is a `400` naming the rule, and an unknown host is a `404`.
An orbital facility's period and speed come from its host's mass the way a
moon's do (a close binary's pair is orbited as one); a gas giant takes
orbital facilities only. A facility reads back with every column of its
`facilities` row plus `host_id` and `host_name`.

**Type, affiliation and ratings (POP.10).** `POST /api/facilities` and
`PATCH /api/facilities/<id>` also take:

- `type_id` — a facility type of the facility's own kind that is not
  retired. `fields` — `{key: value}` for that type's fields, each checked
  against its definition (a whole number for `int`, one of the options for
  a `choice`, and so on). A field left out takes its default; a required
  field with neither is a `400`. Changing the type drops the old type's
  values.
- `affiliation_polity_id` (a polity; the facility keeps the polity's name
  in `affiliation_name`, so a rescan that renumbers polities loses nothing)
  or `affiliation_name` (free text, no polity); not both.
- `ratings` — `{crime, housing, resources, maintenance, health}`, each
  `green`, `yellow`, `red` or `null`. Saving ratings sets `ratings_source`
  to `admin` and stamps `ratings_set_at`. `ratings_note` is free text.

`PATCH` changes only what it names, and `null` clears a value (a required
field cannot be cleared). A facility's place never changes. It returns the
facility.

A facility reads back with every column of its `facilities` row plus
`host_id` and `host_name`, and `type_id`, `type_name`, `type_icon`,
`affiliation_polity_id`, `affiliation_name`, `fields` (`{key: value}`),
`ratings` (a colour or `null` for each of the five), `rating_words` (what
the colour means for that dimension: Crime Low, Elevated, High; Housing
Ample, Tight, Overcrowded; Resources Sufficient, Strained, Critical;
Maintenance Sound, Deferred, Failing; Health Healthy, At risk, Crisis),
`overall` (the worst assessed rating), `ratings_source`, `ratings_note` and
`ratings_set_at`.

`GET /api/facilities` lists facilities, a page at a time (`limit`, `offset`;
`total` counts what passes the filters). Filters: `kind`, `type_id`,
`affiliation_polity_id`, `affiliation_name`, `rating_crime`,
`rating_housing`, `rating_resources`, `rating_maintenance`, `rating_health`
(each `green`, `yellow` or `red`) and `overall`. `sort` is `name` (the
default), `kind`, `type` or `worst`, which lists the worst overall rating
first and the unassessed last; `order=desc` reverses the ratings.

### Facility types

A facility type refines a kind ("Agricultural colony" is a `colony`) and
carries up to 24 typed fields an admin defines. Anyone may read them; an
admin writes (scope `admin`, current credentials).

- `GET /api/facility-types[?kind=][&active=1]` returns `{"items", "icons"}`:
  every type in sort order with its `fields` and how many `facilities` use
  it, and the icon names a type may take. `active=1` leaves out retired
  types. `GET /api/facility-types/<id>` returns one type.
- `POST /api/facility-types` takes `name` (unique) and `kind` (required),
  `icon`, `description`, `sort_order`, `active` and `fields` (a list of
  field definitions). It returns `{"id"}` with a `201`.
- `PATCH /api/facility-types/<id>` changes any of `name`, `kind`, `icon`,
  `description`, `sort_order`, `active` and returns the type. `DELETE`
  removes it. A type that facilities use keeps its kind and is not deleted
  (`409`); set `active` false to retire it.
- `POST /api/facility-types/<id>/fields` takes `key` (a lowercase letter,
  then up to 31 lowercase letters, digits or underscores; never changes),
  `label` and `datatype` (`int`, `decimal`, `text`, `bool` or `choice`),
  plus `unit`, `required`, `min_value` and `max_value` (numbers only),
  `choices` (a `choice`'s 2 to 32 distinct options), `default_text` (the
  value as text, checked like any value) and `sort_order`. It returns
  `{"id"}`. `PATCH` and `DELETE /api/facility-types/<id>/fields/<field_id>`
  change a definition and remove a field.

Values already stored stay valid: a field some facility holds a value for
cannot be deleted, change datatype, drop an option in use or narrow its
range past a stored value (`409`). A required field added to a type that
facilities already use needs a `default_text`, which every one of them takes.
A name or key that is already taken is a `409`; anything else wrong is a `400`
naming the problem.

**Suggested ratings.** `POST /api/facilities/<id>/suggest-ratings` (no body)
fills the ratings nobody set, from what is stored: the host's habitability
for a surface facility (a sealed orbital one counts as 0.6), the six tech
indices of the affiliated species, the kind of facility, what the host
offers, and how far the system is from its polity's capital. Each score
(0 to 1, plus a stable draw of up to 0.1 either way seeded from the facility)
is Green at 0.66 or more, Yellow at 0.33 or more, Red below. It never changes
a rating that is set; `ratings_source` becomes `generated` when no rating was
set by hand, and stays `admin` otherwise. It returns `{"filled", "facility"}`,
where `filled` lists the dimensions it set. `POST /api/facilities/suggest-ratings`
does the same for every facility with an unset rating (at most 5,000, lowest
ids first) and returns `{"facilities", "ratings"}`.

`POST /api/facility-types/examples` (no body) adds the example type, "Research
outpost", unless it is there already, and returns `{"added": [names]}`.

Editing individual generated bodies (stars/planets/moons/belts) isn't
supported via this API beyond renaming them; regenerate the whole system
in place instead (see "Regenerating a system" above).

## Queued jobs

A route that would hold a web worker for minutes queues its work on Redis
(PERF.24) and answers `202 Accepted` with `{"status": "accepted",
"job_id": ..., "status_url": "/api/jobs/<id>"}` and a `Location` header.
Each job has its own queue and a burst worker started for it, detached
from the web server, so a proxy timeout or a reload doesn't stop it.

Quick generating edits queue the same way but wait up to eight seconds
for the job first, so they normally still answer in the same response
(`201` or `200` as before); only a slow one answers `202`. These are
`POST /api/systems`, `PATCH /api/systems/<id>` with `regenerate`,
`POST /api/{planets,moons,belts}/<id>/regenerate`, `POST
/api/phenomena/<type>/<id>/regenerate`,
and `POST /api/systems/<id>/star` (`201` when it finishes in time). Deletes and plain renames stay in the
request. A refusal comes back with its own status (404, 409, 400); `503`
means no Redis server answered.

- `GET /api/jobs/<id>` (admin) — `{"id", "state", "result", "error",
  "error_status", "made_url", "waiting"}`; `error_status` is the HTTP status of a
  refusal; `waiting` is, for a running job that waits in line for the generation run in front of it
  (PERF.153, PERF.154: the one waiting list holds up to `queue.max_waiting` entries, 8 by default, and refuses one more with the list named; `queue.wait_minutes`, 60 by default, then it fails naming that run), `{"run", "position", "text"}`
  (`text` reads "Waiting at position 1 for run ..."), else `null`; `made_url` is, for a finished job that generated sectors
  (a neighborhood or a regenerate), the Galaxy Map fitted to them and
  highlighting them, else `null`.
  `state` is `queued`, `running`, `succeeded` or `failed`; `result` is
  what the work returned once it succeeded; `error` is the last line of
  its error text once it failed. `404` for an unknown or expired id
  (results are kept for a day), `503` when no Redis server answers.

### Slow page answers (PERF.105)

A page whose answer is slow to build does not build it in the request: it queues one job that stores the answer in the
page cache and polls these routes (design in [slow-page-answers.md](../design/operations/slow-page-answers.md)).

- `GET /api/table-estimates?db=` — `{"tables": {name: rows}}`, the database's own estimate (from `information_schema`;
  instant, within some tens of percent) of every table's row count. The pages read it about once a minute to tell which
  requests are slow.
- `GET /api/page-jobs/<id>` — `{"id", "state", "done", "total", "unit", "estimated", "eta_seconds", "seconds",
  "cancelling", "error", "cache_key", "target"}`. `state` is `queued`, `running`, `succeeded`, `failed` or `cancelled`.
  `done` of `total` are rows and an estimate (`estimated` is always `true`): the share of the expected time that has
  passed, never past 95 % before `succeeded`; `eta_seconds` is the time left, `null` when it is past the expected time
  or unknown. `cache_key` is where the finished answer is stored and `target` the request it builds. Public, like the
  pages it serves, and only for page jobs. `404` for an unknown or expired id (kept an hour), `503` without Redis.
- `POST /api/page-jobs/<id>/cancel` (admin) — stops the job if its query has not started, and answers the `GET` body. A
  query already running finishes and its answer is kept. Every visitor waiting on the request shares the job.

## Rate limiting

Every route is rate-limited per client IP via
[Flask-Limiter](https://flask-limiter.readthedocs.io/), on top of which
every write endpoint applies its own stricter limit
(`routes.WRITE_RATE_LIMIT`, currently 10/minute):

- Default: **200 requests/day, 50/hour**, per client IP — Flask-Limiter's
  own quickstart example limit, a reasonable starting point for a
  low-traffic public API with no other usage data to tune against yet.
  Override with `PLANETGEN_RATELIMIT_DEFAULT` (semicolon-separated, e.g.
  `"1000 per day;200 per hour"`).
- Storage backend: the Redis server in `redis.url` by default
  (`ratelimit.storage_uri` empty; `PLANETGEN_RATELIMIT_STORAGE_URI`
  overrides it), shared by every worker process, and by the login
  lockouts. `memory://` is right only for a single-process deployment
  (Flask's dev server, or `mod_wsgi`/`gunicorn` with exactly one
  process): each worker would otherwise track its own separate counters
  and the real, aggregate request rate could exceed the configured limit
  by roughly the worker count. A Redis outage falls back to counting in
  memory until it returns.
- `/api/health` has its own limit instead of the default:
  `ratelimit.pages.health` in `config.json`, default **60/minute** (so a
  monitor can poll it freely but nobody can use it to hammer the database).
  The HTML pages have per-IP limits of their own too (`ratelimit.pages`,
  see [`html-interface.md`](html-interface.md#rate-limits)); the API calls
  a page makes in-process count against neither.
- Exceeding a limit returns `429` with a `Retry-After` header and
  `X-RateLimit-*` headers (`RATELIMIT_HEADERS_ENABLED`).
- "Client IP" is `request.remote_addr`. Behind a separate reverse proxy
  (nginx, Caddy, Apache `mod_proxy`) that is the proxy's address
  unless `config.json`'s `proxy_fix` is set (see
  [`config.md`](config.md) and
  [`deployment/README.md`](../deployment/README.md#behind-a-reverse-proxy-proxy_fix));
  without it every visitor shares one budget. Under Apache + `mod_wsgi`
  it is already the client's address.

## Running locally

```bash
pip install -e .[api]
python src/html/wsgi.py
```

Connects to the same MySQL database every other tool in this project
defaults to (`PLANETGEN_MYSQL_*` env vars, `config.json`'s `mysql` section,
or their built-in defaults — see `planetgen.db.store.MySQLConfig` and
[`config.md`](config.md)). Point it at a different database with either a
`config.json` at the repo root or:

```bash
PLANETGEN_MYSQL_HOST=db.example.com PLANETGEN_MYSQL_DATABASE=planetgen_alpha python src/html/wsgi.py
```

Every other `PLANETGEN_*` variable mentioned throughout this document
(rate limits, the control-schema name, the admin cookie's `Secure` flag)
has a matching `config.json` field too — see [`config.md`](config.md) for
the full list; env vars still take precedence over `config.json` when
both are set.

Every read *and* write (sector/system create/update/delete) goes through
the same `PLANETGEN_MYSQL_USER`/`PLANETGEN_MYSQL_PASSWORD` account —
`WRITE_MYSQL_CONFIG` simply reuses `MYSQL_CONFIG` (see
`planetgen/api/config.py`), there's no separate write-capable override. Point
it at an account with `INSERT`/`UPDATE`/`DELETE`/`SELECT` grants in
production — see [`deployment/README.md`](../deployment/README.md#mysql-accounts).

Admin logins/sessions/API keys/audit log live in a separate **control
schema** (`PLANETGEN_CONTROL_DATABASE`, default `planetgen_control`),
global to the deployment rather than per-galaxy-database — see
`planetgen/db/control_schema.sql`'s header comment and
[`database-schema.md`](database-schema.md). `planetgen.cli.migrate` creates and
seeds it alongside its usual content-schema migration.

### The first admin login

There is no published default password. The first time `planetgen.cli.migrate`
runs against an empty control schema (so on the first `install.sh`), it
creates the user `admin` with a random password and prints both once, in
a boxed block in the installer's output. Only the password's hash is
stored, so nothing can show it again. Log in at `/login` with it; you are
sent straight to `/account` to choose your own username and password, and
every other admin page and write endpoint refuses to work until you do.

#### Resetting the admin login

If that password is lost (or every admin is locked out), delete the admin
rows and let `planetgen.cli.migrate` seed a new one:

    mysql -e 'DELETE FROM planetgen_control.admin_users'
    python3 -m planetgen.cli.migrate        # or ./update.sh

(use your `PLANETGEN_CONTROL_DATABASE` name if you changed it). Deleting
the rows also deletes their sessions and API keys; the audit log keeps
its entries. The new random password is printed as above.

## Deploying

`src/html/wsgi.py` exposes the standard WSGI `application` object: the
same app serves the API under `/api` and every HTML page. Every hosting
setup loads that one file:

- **Apache + `mod_wsgi`** (the reference setup, `examples/apache/`): a
  `WSGIScriptAlias` for `/` pointing at `src/html/wsgi.py`, in its own
  `WSGIDaemonProcess`. See
  [`deployment/apache.md`](../deployment/apache.md).
- **gunicorn** behind nginx or Caddy (Linux, macOS):
  `gunicorn --pythonpath <checkout>/src/html wsgi:application`.

[`deployment/README.md`](../deployment/README.md) compares them and links
each guide. Behind nginx, Caddy or Apache's `mod_proxy`, set
`config.json`'s `proxy_fix` so the rate limits see the client's address
and the app knows a request came over HTTPS; under `mod_wsgi` leave it
off. No separate vhost or `ServerName` is needed for the API: the pages
call it in-process (`planetgen/web/transport.py`), not over HTTP.

Point the app at the deployed MySQL database with `config.json`'s `mysql`
section, or `PLANETGEN_MYSQL_*` in the process environment. That account
needs write grants, not just `SELECT`: the admin pages and every write
endpoint write through it (see "Running locally" above and
[`deployment/README.md`](../deployment/README.md#mysql-accounts)). Under
`mod_wsgi`, environment variables must be in the Apache service's own
process environment (e.g. `/etc/apache2/envvars`, or an `Environment=`
line on the apache2 systemd unit), **not** the vhost's `SetEnv`
directives: those never reach `os.environ` under `mod_wsgi` --
`planetgen/api/config.py` reads its config from `os.environ` once, at process
startup, and `SetEnv` values only ever show up in a request's `environ`
dict, which doesn't exist yet at that point. Under gunicorn,
the service's own environment (systemd `EnvironmentFile`, launchd
`EnvironmentVariables`) works the normal way. If running
more than one process, leave `ratelimit.storage_uri` empty so the counts
are on Redis (see "Rate limiting").

**The admin session cookie requires HTTPS.** It's set `Secure` by default
(`config.SESSION_COOKIE_SECURE`) -- the browser never sends it over plain
HTTP, so `/login` and `/admin` won't work on an HTTP-only site. Terminate
TLS in front of the app (e.g. `certbot --apache` or `--nginx`, or Caddy's
automatic HTTPS) before using the admin pages; only set
`PLANETGEN_ADMIN_COOKIE_INSECURE=1` for local development without TLS in
front, never in production.

## Not done yet

Planned API changes, as of 2026-10-15. None is built; the
items and their full text are in `docs/TODO.md`, and the phases in
`docs/plan/`.

- **Every API call logged (API.15, phase 1).** By default every call is
  logged with the time, the route, the account that made it (the key's
  owner for an API key, the signed-in admin for the web site, and `god`
  for the console), how it came in (API key, web session or console),
  and the HTTP response code it got back. The log uses HTTP response
  codes for parity with the web server's own logs. It needs no user
  accounts: it logs today's admin accounts, and any account once user
  accounts exist (USR.1). The activity log
  ([`config.md`](config.md#the-activity-log)) stays as it is.
- **Remote generation (API.3 and its parts, phases 2 and 3).** The
  download of the seed, skeleton and naming key (API.12), run
  reservations (API.10), staging tables (API.11), compressed batch
  uploads (API.14) and checks on them (API.8); a remote run with the same
  seed and release produces what the server would, checked by
  fingerprint (API.17).
- **Generate by recipe (API.18, phase 2; API.19, phase 3).** A JSON
  recipe, where any field can be fixed, ranged or left random, generates
  a sector, a system, a planet, a moon or any phenomenon, and later a
  whole galaxy built up region by region. A request that makes no sense
  gets 400; one that fails validation gets 422; both return a JSON error
  object with details and the command's log output.
- **Keys owned by accounts (API.6, phase 3+)**, after user accounts
  (USR.2).

Editing a generated star beyond changing its class or renaming it is also
still open (planets, moons and belts can be deleted and regenerated one at a
time, see "Deleting and regenerating"). The pages in `../src/planetgen/web/` (see
[`html-interface.md`](html-interface.md)) are this API's own
server-rendered frontend.
