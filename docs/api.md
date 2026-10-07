# planetGen API

A JSON API over the planetGen database (`src/stellarObjects/schema.sql`),
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
publish to a wiki, generate a sector's neighborhood) requires an authenticated admin whose credentials aren't still the
seeded first login — see "Authentication" and "Write endpoints" below.

## Why Flask

Comparison against FastAPI/Django REST Framework: the persistence layer
(`stellarObjects/_db.py`) is deliberately plain SQL over a small `pymysql`
wrapper with no ORM, this API is read-heavy with no concurrency pressure yet,
and it needs to deploy onto a plain Apache2/VPS setup (or nginx, Caddy,
IIS or macOS; see [`deployment/`](deployment/README.md)).
Flask has no opinion
about the data layer (route handlers call straight into `queryDb.py`'s and
`stellarObjects._db`'s existing functions), deploys as a plain WSGI app
(`mod_wsgi`, gunicorn or waitress), and
lives at `../src/html/api/`, mounted at `/api/` by the same Flask app that
serves the HTML pages (`../src/html/web/`, see "Deploying"
below). Those pages are this API's own frontend, calling it in-process. FastAPI's
headline advantages (async, auto-generated OpenAPI docs) still don't pay
for themselves: this API is read-heavy and low-concurrency regardless of
framework, and the pages are plain server-rendered HTML with no use for
generated API-contract docs the way a JS single-page app
would. Worth revisiting if a richer JS frontend is ever built against
this API instead.

## Endpoints

All under `/api/`, all JSON in, JSON out. Every endpoint below except
`/databases` accepts an optional `db=<name>` query
parameter selecting *which* MySQL schema on the configured server to
read from (validated against the same prefix-filtered list
`/api/databases` itself returns — an unrecognized name is a `404`, same
as an unrecognized sector/system id); omitted, it falls back to
`MYSQL_CONFIG`'s own configured default database (`config.py`). This is
what lets one API process serve several databases — see
`stellarObjects._db.list_databases`/`resolve_database`. (The HTML pages
always show the one configured database.) `/api/health` also honors `db=` — passing it checks
connectivity to that specific schema rather than the default one.

### Read

- `GET /api/health` — liveness/readiness check: confirms the process is up
  and the configured (or `db=`-selected) database can actually be opened.
  Returns `{"status": "ok", "schema_version": <n>, "schema_current":
  <bool>}` (plus a `detail` naming the fix when the schema is behind the
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
  listed or selectable with `db=`) (`stellarObjects._db.list_databases`),
  each with `name`, `size_bytes`, `modified_at`, and a quick-glance
  `sector_count`/`system_count` (`null` for a matching schema missing this
  project's own tables, e.g. mid-migration, rather than failing the whole
  listing).
- `GET /api/sectors?limit=<n>&offset=<n>` — every sector, paginated (see
  "Pagination" below), each with `id`, `name`, `edge_mpc`, `edge_ly`,
  `system_count`, and its galaxy placement (`center_x_pc`/`center_y_pc`/
  `center_z_pc`/`galactic_radius_pc`/`galactic_radius_ly`/`ring_index`/
  `layer_index`/`ring_slot_index` (the cylindrical grid address), all
  `null` together for an unplaced sector, plus a `placed` bool), nearest
  the galactic core first, then unplaced sectors by name
  (`queryDb.list_sectors`/`count_sectors`).
- `GET /api/sectors/<id>` — one sector's full display detail: the same
  fields as the listing above, plus `systems` (every system placed in
  it, nearest the sector's center first: `id`, `name`, `quadrant`, `location`, `is_binary`, `binary_type`,
  `position_x_mpc`/`position_y_mpc`/`position_z_mpc`, `center_distance_ly`,
  `runaway_class`/`runaway_speed_kms`, `inside` (the nebula or supernova
  remnant around it, or `null`), `nearest` (its stored nearest systems,
  `{id, name, distance_ly}`), and `stars`, each
  with `role`/`name`/`star_type`/`temperature_k`/`radius_km`/`luminosity_w`),
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
  relative to this sector's own, and `home` (`true` when generated as
  part of this sector, `false` for a neighbor's cloud) — `queryDb.phenomena_near_sector`, empty
  for an unplaced sector; see `schema.sql`'s "v18"/"v21"/"v28" header
  notes),
  `wiki_url` (`null` until this sector has a wiki page — see "Wiki
  publishing" below), and `stats` (PERF.11, GEN.44: the sector's
  `sector_stats` row — `bright_level_sol`, `relative_density`,
  `expected_systems`, `actual_systems`, `actual_stars`,
  `mean_temperature_k`, `mean_luminosity_sol`, `fill_share` and its Galaxy
  Map color `color_r`/`color_g`/`color_b` (MAP.86) — or `null` for a sector off
  the grid) (`queryDb.sector_detail`). Distinct from
  `stellarObjects._db.load_sector(...).to_dict()`'s *generation* object
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
  from omitting `sector_id` entirely (no sector filter at all).
- `GET /api/systems/<id>` — one system's full display detail: `id`,
  `name`, `sector_id`, `quadrant`, `location`, `is_binary`, `binary_type`,
  `binary_configuration` (`"close"`/`"wide"`/`null`),
  `binary_mutual_position_x/y/z_km` (the secondary's position relative to
  the primary — `null` for a single star; the System Map's own real
  binary-star placement is derived from this plus each star's `mass_kg`),
  `runaway_class` (`"runaway"`/`"hypervelocity"`/`null`) and
  `runaway_speed_kms`, `wikijs_url`/`mediawiki_url`
  (each `null` until this system has been uploaded to that wiki — see
  "Wiki publishing" below), `stars`, `planets` (each with
  its own nested `moons`; every planet and moon also carries `habitable`,
  `life_stage` — the furthest evolutionary milestone its timeline reached,
  or `null` — and `inhabited`), `belts` (each with its `composition`,
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
  `stellarObjects._db.load_star_system(...).to_dict()` as `/api/sectors/<id>`
  above.
- `GET /api/systems/<id>/text?format=wikitext|markdown` — the system's
  full wiki page, `{"id", "format", "content"}`, rendered from its
  database rows on each request (`stellarObjects/systemRender.py`; no page
  text is stored since schema v29). `format` defaults to `wikitext`; any
  other value is a `400`.
- `GET /api/systems/<id>/sections` — the same page as Markdown split for
  the system page's expandable list: `overview` (a binary pair's own data
  table, the system summary, flavor text) plus `stars`, `planets`,
  `moons`, `belts` and `comets`, each an object mapping a row `id` (as a
  string) to that body's own Markdown, without its heading (and, for a
  planet, without its moons' sections).
- `GET /api/systems/<id>/near?radius=<ly>` — other systems in the same
  sector within `radius` light-years (`queryDb.systems_within_radius`).
- `GET /api/nav?from=<id>&to=<id>` — course, distance, and an optimal route
  between two systems (`queryDb.nav_between`) — see "NAV" below.
- `GET /api/galaxy/sectors` — every galaxy-placed sector (non-`null`
  galaxy placement), each with `id`, `name`, `x`/`y`/`z`
  (`center_x/y/z_pc`), `galactic_radius_pc`, `ring_index`, and
  `system_count` (`queryDb.galaxy_placed_sectors`) — the data
  the Galaxy Map page (`/galaxy`, `../src/html/web/galaxy_views.py`) plots. Not paginated: bounded by
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
  density-skeleton shape (`generate.py plan`'s output): every
  `planetgen.galaxy.density.GalaxyShape` field plus `edge_pc`,
  `outer_ring_index`, and `expected_system_count_at_density_1`
  (`queryDb.galaxy_density_shape`). `shape` is `null` when the skeleton
  has never been built. The Galaxy Map shades its "expected density"
  cloud from this real model (falling back to a generic illustrative
  gradient when `null`) instead of a placeholder.
- `GET /api/galaxy/cell?ring=<i>&layer=<j>&slot=<k>` (or `?x=&y=&z=`,
  parsecs, for the cell holding that point) — one cell of the cylindrical
  sector grid, whether or not anything was generated there
  (`galaxyGeometry.describe_sector_cell`): `ring_index`/`layer_index`/
  `ring_slot_index`, `designation`, `cartesian_pc`, `cylindrical`
  (`r_pc`, `theta_rad`, `z_pc`), `spherical` (`r_pc`, `theta_rad`,
  `polar_rad` from galactic north), `bounds`, `mean_arc_length_pc`,
  `volume_pc3`, `vertices_pc` (8 corners), `edge_pc`, and `sector_id`
  (the generated sector there, or `null`). `400` for a bad or incomplete
  query.
- `GET /api/galaxy/tiles?tiles=<key>,<key>,...` — the 3D
  Galaxy Map's data, one fixed cube of space ("tile") at a time
  (`queryDb.galaxy_tiles`). Space is an octree: level 0 is one cube
  65,536 pc on a side centered on the galactic origin, each level halves
  the edge down to 16 pc at level 12, and a key is `level/ix/iy/iz`
  (`ix` counts cubes along x from the root cube's −x face). Returns
  `{"tiles": {"<key>": {"placed": [...], "planned": [...], "filled": {...},
  "clouds": [...], "stars": [...], "generated": [...], "points": [...]}},
  "edge_pc": ...,
  "has_shape": ...}`:
  `placed` is every placed sector whose center is in the tile's half-open
  box (the `/api/galaxy/sectors` shape plus `layer_index`, `ring_slot_index`,
  `designation` and `edge_ly`, this sector's real edge length, `null` if it
  predates per-sector edge tracking), lowest id first,
  at most 250; `planned` lists the tile's qualifying not-yet-generated
  slots, only for 16 pc tiles (empty otherwise); `filled` counts every
  placed sector in the tile, nothing sampled (`queryDb.galaxy_filled_in_box`):
  `{"g": 1, "cells": [[ring, layer, slot, id, system_count, name], ...]}`
  one per sector, or, for large tiles or more than 5,000 sectors, `{"g": g,
  "cells": [[ring, layer, wedge, count], ...]}` per cell `g` sectors a side
  (`g` a power of 3; cell ring `I` has `max(3, round(2π(I + ½)))` wedges
  from +x, and a sector counts in the wedge holding its center angle).
  The Galaxy Map sums these into its blocks. `clouds` lists the placed
  nebulae and supernova remnants whose sphere reaches into the tile
  (`queryDb.galaxy_clouds_in_box`), largest first, at most 200:
  `{type, id, name, descriptor, class, radius_pc, x, y, z}`. `stars`
  lists the tile's pre-placed bright stars (`bright_stars`, filled or
  not; `queryDb.galaxy_bright_stars_in_box`), most luminous first, at
  most 400: `{id, x, y, z, luminosity_sol, temperature_k, radius_sol,
  star_type, yerkes_class, ring_index, layer_index, ring_slot_index,
  system_id}` (`system_id` is `null` until the star's sector is filled).
  `generated` lists the stars of the tile's generated systems (by their
  sector's center; `queryDb.galaxy_generated_stars_in_box`) at or above
  the tile level's luminosity floor (`queryDb.generated_star_floor_sol`:
  every star at level 12, four times brighter per coarser level, none
  past 1,000 L☉), most luminous first, at most 1,000 (4,000 at level 12,
  so the Galaxy Map zoomed to a sector shows all of its stars, MAP.80),
  read from at most
  1,500 of its sectors (a hashed sample past that); a star already in
  `stars` isn't repeated: `{id, name, x, y, z, luminosity_sol,
  temperature_k, radius_sol, star_type, ring_index, layer_index,
  ring_slot_index, system_id}` (`id` is the `stars` row). `points` lists
  the placed black holes, neutron stars and quasars centered in the tile
  (`queryDb.galaxy_point_phenomena_in_box`), only at level 10 (64 pc) and
  finer, most luminous first, at most 200: `{type, id, name, descriptor,
  luminosity_sol, x, y, z}` (`type` is `black_hole`, `neutron_star` or
  `quasar`). Predicted density isn't
  served: the page evaluates the galaxy's shape itself
  (`static/galaxyprisms.js`). At most 128 keys per request; a
  malformed key is a 400. Every part depends only on its key and the
  database's contents, so callers cache it by key and `/api/galaxy/stamp`.
- `GET /api/galaxy/stage?at=<m.ring.wedge.slab>` — one stage of the
  Galaxy Map's drill-down (`queryDb.galaxy_stage`, design in
  `design/galaxy-drilldown-navigation.md`): how many generated sectors
  each child block of block `at` holds, as `{"at", "child_m", "children":
  [{"ring", "wedge", "slab", "generated", "look"}], "sectors"}`. `look`
  is what the generated sectors hold, from `sector_stats` (MAP.86):
  `{share, color, colored}`, their mean `fill_share` (0 for one without
  stats), the mean sRGB `[r, g, b]` of those with stars (`null` if none)
  and how many had one. Blocks are 243,
  27 and 3 sectors a side (`planetgen.galaxy.drill`); with no `at`,
  the children are the galaxy's level-243 blocks. Children with nothing
  generated are left out (the page computes totals itself). At a level-3
  block the children are sectors (`wedge` is the slot, `slab` the layer)
  and `sectors` lists each one as `{ring, layer, slot, id, name,
  system_count, look}`; otherwise `sectors` is `null`. A malformed or impossible
  key is a 400.
- `GET /api/galaxy/locate?q=<part of a name>` — the Galaxy Map address
  bar's name lookup (`queryDb.galaxy_locate`): `{"matches": [{"kind"
  (`"sector"` or `"system"`), "id", "name", "sector_id", "sector_name",
  "ring", "layer", "slot"}]}`, at most 8, exact names first, then names
  that start with the term. Sectors with no address and systems outside a
  sector are left out, since the map can't fly to them.
- `GET /api/galaxy/stamp` — `{"stamp": "<16 hex characters>", "state":
  "<token>"}` (`queryDb.galaxy_content_stamp`). `stamp` changes whenever
  tile contents could: sectors placed, edited or removed (their
  `modified_at`, schema v27), new star systems, a re-planned galaxy shape,
  or a new planetGen release. `state` is what `/api/galaxy/changes` takes.
- `GET /api/galaxy/changes?since=<state>` — `{"stamp", "state", "full",
  "tiles", "stages"}` (`queryDb.galaxy_changes`): the keys of the tiles
  that changed since that `state`, found from `sectors.modified_at`, new
  sector ids and new star-system ids, and the drill-down stages holding
  those sectors (`"galaxy"` and each one's level-243, 27 and 3 block
  keys). `full` is `true` (and both lists empty) when that can't
  be pinned to tiles: a deleted sector, a new shape or release, more than
  1,000 changed sectors, or a missing or unreadable `since`.
  `../src/html/lib/tilecache.py` (the web layer's disk cache) calls it
  about once a minute and deletes only the listed tiles, and passes the
  list on to the map's browser cache.
- `GET /api/phenomena?limit=<n>&offset=<n>` — every exotic phenomenon,
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
  The data the `/phenomena` page shows.
- `GET /api/phenomena/<type>/<id>` — one phenomenon's full detail (every
  column its own table has, e.g. a nebula's `composition`/
  `formation_cause`, a black hole's `mass_solar`/`spin`/
  `has_accretion_disk`, a supernova remnant's `morphology`/`progenitor_type`/
  `age_years`), plus `type`, `sector_name` and `nearest` (its three
  nearest star systems, `{id, name, distance_ly}`). An asteroid field also
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
  now including a to-scale AU diagram (`../src/html/lib/phenomenonmap.py`).
- `GET /api/search?sector_q=&system_q=&star_q=&planet_q=&moon_q=&<facet>=<value>...` —
  the faceted search behind the `/search` page: click-to-filter tags
  (object type; star spectral/luminosity class; planet/moon class, body
  type, supported life chemistry; asteroid belt density; phenomenon type
  and phenomenon class, the class tag written `<type>:<class>` such as
  `phenomenon_class=nebula:D` — repeat a facet
  name for multiple active values, e.g. `class=M&class=K`), a per-entity
  name search, and a min/max size range per entity —
  `star_min_radius_km`/`star_max_radius_km` (likewise `planet_`/`moon_`),
  in km, either bound optional (a size range is a continuous quantity,
  not a discrete facet value, so it's its own pair of query parameters
  rather than a tag). Returns `facets` (one `{value, label, count,
  tooltip}` list per facet, built from the distinct values actually
  present), `autocomplete` (`sectors`/`systems`/`stars`/`planets`/`moons`
  name lists), `facet_labels` (`"facet:value"` -> label, for an
  active-filter chip), and `results` (`sectors`/`systems`/`stars`/
  `planets`/`moons`/`belts`/`phenomena` -> `{"rows": [...], "total", "limit",
  "offset", "truncated"}`, or
  `null` for an object type with no active reason to query it — see
  `queryDb.search`'s docstring for the exact inclusion rule; a size range
  alone is reason enough, same as a tag or name term). `stars`/`planets`/
  `moons` result rows each include their own `radius_km`. Each result
  panel is paged on its own: `limit` sets the rows per panel (default
  300, clamped to 500 like "Pagination" below), and `sectors_offset`/
  `systems_offset`/`stars_offset`/`planets_offset`/`moons_offset`/
  `belts_offset` pick each panel's page; `total` is that panel's full
  match count and `truncated` is true when `rows` isn't all of them. An
  offset past the last match returns the last page (with its real
  `offset`).
- `GET /api/wiki-config` — `{"wikijs": bool, "mediawiki": bool}`, whether
  each wiki backend has a `base_url` plus credentials configured
  deployment-wide (`config.json`'s `wiki` section, or the matching
  `PLANETGEN_WIKIJS_*`/`PLANETGEN_MEDIAWIKI_*` env vars — see
  `docs/config.md`) and so is offered as an upload target at all. Never
  exposes any of those credentials, only the two booleans.
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
  there, which `updateOrbits.py` advances like an orbital facility's.
- `GET /api/galaxy/bright-stars?ring=<i>&layer=<j>&slot=<k>[&all=1]` —
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
- `GET /api/species?spacefaring=true|false&limit=<n>&offset=<n>` — the
  dominant species of every life world, by name, paginated (see
  docs/design/population-and-politics.md). Each item has its homeworld
  (`homeworld_planet_id`, `homeworld_name`, `star_system_id`,
  `system_name`), `life_chemical`, `life_stage` (`multicellularity` or
  `technological_civilization`), `build`/`climate`/`size`,
  `civilization_age_years` and `era` (`null` without a civilization),
  `spacefaring`, and its `polity_id`/`polity_name` (`null` unless
  spacefaring). `400` for any other `spacefaring` value.
- `GET /api/species/<id>` — one species; `404` if unknown.
- `GET /api/planets/<id>/species` — the species whose homeworld that
  planet is; `404` when it has none.
- `GET /api/polities?limit=<n>&offset=<n>` — every polity (one per
  spacefaring species), by name: `name`, `government`, `color`
  (`#rrggbb`), `reach_ly`, its species, `era`, capital and
  `system_count`.
- `GET /api/polities/<id>?limit=<n>&offset=<n>` — one polity plus a page
  of the systems it owns (`id`, `name`, `distance_ly` from the capital),
  nearest first; `404` if unknown.
- `GET /api/systems/<id>/owner` — `{"owner": {polity_id, polity_name,
  color, distance_ly}}`, or `{"owner": null}` when no polity holds it;
  `404` for an unknown system.
- `GET /api/territories` — `points`: up to 20,000 owned systems with
  galaxy-frame positions in parsecs (`x`, `y`, `z`) and their polity's
  `color`, nearest their capitals first; `polities`: each polity's `id`,
  `capital_pc` and `reach_ly`. For a 3D territory overlay.

### Write (admin auth required — see "Authentication" and "Write endpoints")

- `POST /api/sectors` — create a sector.
- `PATCH /api/sectors/<id>` — modify a sector (including manually
  setting/clearing its `wiki_url`).
- `DELETE /api/sectors/<id>` — remove a sector.
- `POST /api/sectors/<id>/wiki` — publish a sector-summary page to a
  wiki (see "Wiki publishing" below).
- `POST /api/sectors/<id>/generate-neighborhood` — generate every
  not-yet-generated sector within `radius_ly` (optional JSON body field,
  default 12 pc, at most the generator's own radius cap) of this
  galaxy-placed sector, synchronously. Returns counts: `generated`,
  `already_existed`, `skipped`, `candidates` and `outside_galaxy`, plus
  the size and time `estimate` (see [`cli.md`](cli.md#size-and-time-estimates)).
  `"estimate_only": true` returns the counts and `estimate` without
  generating anything; a run the database disk can't hold is refused
  with `507` and nothing written. `404`
  for an unknown or unplaced sector, `409`
  when the galaxy has never been planned (`generate.py plan`). Every new
  sector also gets the bright stars (100 L_sun and up) within 100 ly of it
  (GEN.23). A large radius (100 ly) covers thousands of candidate slots,
  so this can run for a long time; a reverse proxy's timeout may need raising for it. The sector
  page's admin form calls it.
- `POST /api/systems` — generate and create a system, standalone or in
  an existing sector.
- `PATCH /api/systems/<id>` — rename a system (its stars, planets and
  moons follow; see "Renaming" below), and/or regenerate its contents
  in place (see "Regenerating a system" below).
- `PATCH /api/stars/<id>` — rename a star.
- `PATCH /api/planets/<id>` — rename a planet.
- `PATCH /api/moons/<id>` — rename a moon.
- `DELETE /api/systems/<id>` — remove a system.
- `POST /api/systems/<id>/wiki` — publish a system's already-generated
  page to a wiki (see "Wiki publishing" below).
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
- `GET /api/systems/<id>/class-options`, `POST /api/planets/<id>/class`,
  `/api/moons/<id>/class` — change a planet's or moon's class; `POST
  /api/systems/<id>/star` — change a single-star system's star (see
  "Changing a class or a star" below).

### Admin stats (admin auth required)

Both take `?db=` like the read endpoints, and need an admin past the
forced credential change. They back the admin stats page
(`/admin/stats`, `../src/html/web/admin_pages.py`).

- `GET /api/admin/stats` — health and statistics for one database:
  `api` (version, Python version, process uptime, load average, memory),
  `mysql` (server version, uptime, connected threads; `null` fields when
  `SHOW GLOBAL STATUS` isn't allowed), and `database` (`reachable`,
  `schema_version`/`schema_expected`/`schema_current`, `size_bytes`,
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
  "max_density"}], "sizes": {database: {"bytes_per_system", "systems",
  "total_bytes"}}, "available"}` (`available` is false until `update.sh`
  has added control schema v6). Shown on the Stats page.
- `GET /api/admin/work[?limit=&offset=]` — the work queue (ADM.10):
  `{"available", "status": {"workers_active", "runs_active", "paused",
  "paused_by", "paused_at", "holder", "job_id", "heartbeat_age_s",
  "stale", "load": {"values", "kind", "text"}}, "items": [job tree
  roots, newest first, each with "totals"], "total", "limit",
  "offset"}`. `load.kind` is `load` (the 1, 5 and 15 minute load
  average) or, on Windows, `cpu` (CPU percent over the same windows).
  `available` is false until `update.sh` has added control schema v7.
- `GET /api/admin/work/<id>` — the whole job tree holding node `<id>`
  (ADM.12): every node's `kind`, `title`, `state`, `status` (`interrupted`
  for a run that stopped answering), Unix-time `created_at`,
  `started_at`, `finished_at`, `seconds`, `control`, `children`, a queue
  node's `tasks` (up to 200, failed first) and `totals` (`tasks`,
  `queued`, `running`, `done`, `failed`, `cancelled`, `work_seconds`,
  `eta_seconds`) added up from the children.
- `POST /api/admin/work/<id>/control` `{"action": "pause"|"resume"|"cancel"}`
  — asks a live node and everything under it to pause (running tasks
  finish, then it stands by without the lease, so other jobs run),
  resume, or cancel. `POST /api/admin/work/<id>/delete` deletes a
  finished tree by its root. `POST /api/admin/work/queue`
  `{"action": "pause"|"resume"}` pauses the whole queue (no run takes the
  lease or a task) or resumes it. `POST /api/admin/work/lease/clear`
  frees a lease whose holder stopped refreshing it. Each answers
  `{"ok": bool}` (`{"cleared": holder}` for the lease) and writes the
  change to the audit and activity logs (`work.pause`, `work.queue.pause`,
  ...).
- `GET /api/admin/lockouts` — every address and username locked right
  now: `{"items": [{"scope", "subject", "retry_after", "locked_until",
  "level"}], "proxy_warning"}`. `proxy_warning` is true when the site
  looks like it is behind a reverse proxy without `proxy_fix` (a private
  address is locked, or most recent failures share one private or
  loopback address).
- `POST /api/admin/lockouts/lift` `{"scope": "ip"|"user", "subject"}`, or
  `{"all": true}` — lifts lockouts and forgets their counts; returns
  `{"lifted": n}`; audited as `lockout.lift`. From a shell (for an admin
  locked out of the site itself): `python3 src/loginLockouts.py` lists
  them, `--ip <address>`, `--user <name>` or `--all` lifts them, and
  `--forget-devices <name>` revokes that admin's trusted-device cookies.

### Authentication

- `POST /api/auth/login` `{"username", "password"}` — sets the session
  cookie (`pg_admin_session`; `HttpOnly`/`Secure`/`SameSite=Strict`) and
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
  address and per username, in the control database's `login_throttle`
  (shared by every worker, kept across restarts; SEC.1, SEC.21). Three
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
  checked. Until `update.sh` has created `login_throttle`, the counts are
  kept in memory per worker process (with one warning in the error log).
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
- `GET /api/auth/totp` — `{"enabled", "recovery_codes_left"}`.
- `POST /api/auth/totp/setup` `{"current_password"}` — `{"secret", "uri",
  "qr_svg"}` for an authenticator app; nothing changes at sign-in yet.
- `POST /api/auth/totp/confirm` `{"code"}` — turns it on once a code
  from the app matches; returns `{"recovery_codes": [...]}`, shown only
  this once (10 codes, each good for one sign-in).
- `POST /api/auth/totp/disable` `{"current_password", "code"}` — turns
  it off, and forgets every trusted device of that admin (the browser
  that turned it off gets a new device cookie). Lost the phone and the recovery codes? From a shell:
  `python3 src/loginLockouts.py --reset-two-factor <user>`.
- `POST /api/auth/logout` — ends the current session, clears the cookie.
- `GET /api/auth/me` — the calling admin's identity.
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
- `GET /api/auth/api-keys` — the calling admin's own API keys (label/
  timestamps only, never the key or its hash).
- `POST /api/auth/api-keys` `{"label"}` — creates a key, returning
  `{"id", "label", "key"}`. `key` is the raw value, shown exactly once.
- `DELETE /api/auth/api-keys/<id>` — revokes one of the calling admin's
  own keys.

Authenticate either with the session cookie (set by `/api/auth/login`,
used by the `../src/html/` admin pages) or an API key, sent as
`Authorization: Bearer <key>`, for programmatic callers.
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
message is just `"database unavailable"`, and an unexpected server-side failure is a
`500` — the API never falls through to Flask's default HTML error page or
leaks a stack trace, a connection error or the database's user/host to
the client (the real detail still reaches Flask's own logger).

Every response also carries `X-Content-Type-Options: nosniff`,
`X-Frame-Options: DENY`, `Referrer-Policy: no-referrer` and
`Content-Security-Policy: default-src 'none'`, and, when the request came
in over HTTPS, `Strict-Transport-Security: max-age=31536000`.

## NAV

`GET /api/nav?from=<system_id>&to=<system_id>` returns a direct course
(distance, bearing and mark, warp and fold travel times) plus an optimal route via
adjacent systems (`planetgen.galaxy.nav_graph`, a k-nearest-neighbor adjacency
graph with Dijkstra shortest-path) between two endpoints -- each either a
star system (the default) or a standalone phenomenon (nebula/asteroid
field/black hole/neutron star/supernova remnant/rogue planet/interstellar
comet/quasar).

**Phenomenon endpoints.** Pass `from_kind=phenomenon&from_type=<type>`
(and/or the `to_*` equivalents) to route to/from a phenomenon instead of a
system -- `from`/`to` then names that phenomenon's own row id, and
`from_type`/`to_type` is one of `nebula`, `asteroid_field`, `black_hole`,
`neutron_star`, `supernova_remnant`, `rogue_planet`, `interstellar_comet`,
`quasar` (an unrecognized type or a nonexistent id is a `404`). A phenomenon
endpoint's own `route.path`/`route.positions` id is a
`"phenomenon:<type>:<id>"` string (a plain int id, same as always, for a
system) -- only `route.path[0]`/`route.path[-1]` can ever be a phenomenon;
every intermediate hop is always a system.

**Availability.** NAV only applies to a pair of endpoints that satisfy one of:

- Both are systems in the SAME sector (`star_systems.sector_id IS NOT
  NULL`, matching) -- scoped to that sector's own systems (sector-local
  positions). Checked first/preferred over the galaxy-scope rule below.
- Otherwise, both have a galaxy-frame position -- a system whose sector has
  a galaxy placement (`sectors.center_x/y/z_pc IS NOT NULL`, i.e. placed by
  `generate.py galaxy`), or a phenomenon that's itself been placed in the galaxy
  (`--sector-id` at generation time -- see `generate.py phenomenon`'s own
  help). A phenomenon endpoint can only ever be resolved this way: it has
  no sector-local position of its own, regardless of which sector its own
  `sector_id` names as its nearest neighbor (that column is a browsing
  convenience, not real containment).
- Anything else (an unplaced system, a phenomenon never placed in the
  galaxy, or two systems in different sectors where either lacks a galaxy
  placement) is a `404` (unknown id) or a `400` with `"NAV requires both
  endpoints to share a sector, or both to have a galaxy placement"`.

**Course convention.** Courses use Boss's nested reference frames
(`docs/design/navigation-frames.md`, `planetgen/galaxy/navigation.py`):
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
  NAV Map plot (`html/lib/navmap.py`) draws.
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
admin (random first password, printed once by `migrateDb.py`) can log in and call `/api/auth/change-credentials`, but
nothing else, until it changes its own credentials (`403` otherwise). A
missing/invalid credential is `401`. Every write route runs against the
same `PLANETGEN_MYSQL_*` database account every read route above uses —
see `config.py` and [`deployment/README.md`](deployment/README.md#mysql-accounts) — and
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

`PATCH /api/sectors/<id>` (any non-empty subset) additionally accepts
`wiki_url`:

```json
{
  "name": "Voranthis Kelmoor",
  "edge_ly": 11.5,
  "wiki_url": "https://wiki.example.com/Voranthis_Kelmoor"
}
```

- `name`: non-empty string.
- `edge_ly`: number, greater than 0 — the sector's cube edge, in
  light-years (matches `sectors.edge_mpc` after unit conversion; see
  `database-schema.md`).
- `wiki_url` (`PATCH` only): an absolute `http://` or `https://` URL with
  a host (at most 2048 characters; anything else, e.g. `javascript:` or
  `data:`, is a `400`), or `null` to clear it back
  to "no page yet" — the manual "set the wiki link directly" admin
  affordance (the `/admin` page); the same column `POST
  /api/sectors/<id>/wiki` (below) writes automatically on a successful
  upload. Rejected as an unrecognized field on `POST` — a brand-new
  sector has never been uploaded anywhere.

An unrecognized field, a missing required field (`POST` only), a wrong
type, or a value failing the constraints above is a `400`.

### Systems — request body

`POST /api/systems` takes a generation **"recipe"**, the same shape
`generate.py system --system-file` already takes (`SystemConfig.to_dict()`/
`from_dict()`) — every field is optional (`null`/omitted means "let the
generator decide"), and the server generates a brand-new, real system from
it the same way `generate.py system` does, via `StarSystem(system_config=...)`:

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
  "num_orbits": 6,
  "markdown": false
}
```

Accepted fields: `markdown`, `star_type`, `name`, `age` (`"young"`,
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
(`sector_id = NULL`, like `generate.py system`). With it, the system is
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
position, location and wiki links, and the new bodies are named from the
system's name. A system in a galaxy sector keeps that sector's stellar
population. It is a `409` when the system was built around a pre-placed
bright star (its star is fixed by the galaxy scatter), or when facilities
are hosted on it, unless `"drop_facilities": true` is sent too (they are
deleted with the bodies). `name` may be sent in the same request (it is
applied first). Success returns `{"status": "ok", "id", "name",
"regenerated"}`.

### Deleting and regenerating

`src/html/api/edits.py` (ADM.8). Every one of these needs an admin whose
credentials are current and writes an audit-log row.

- **Planets, moons and asteroid belts.** `DELETE /api/planets/<id>`
  (with its moons), `/api/moons/<id>` and `/api/belts/<id>` remove one
  body; `POST .../<id>/regenerate` rolls it again at the same orbit: a
  planet gets a new class (any that fits its zone), size, atmosphere,
  life and moons; a moon a new moon class that fits its planet; a belt a
  new density and composition over the same span. A regenerated body
  keeps its name and its row (so facilities on it stay). The rest of the
  system is then re-validated from the moons outward
  (`stellarObjects/validation.py`): bodies are moved outward until the
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
  plan; the new sector has a new id and name (`sector_id`,
  `sector_name`; `null` when the slot is outside the galaxy's outline). It is a `409`
  for a sector off the galaxy grid or before `generate.py plan`.

### Changing a class or a star

Also `src/html/api/edits.py`, with the same admin, audit and answer
shape as above.

- **Class (ADM.6).** `GET /api/systems/<id>/class-options` answers
  `{"recommended": {"planet:<id>" or "moon:<id>": [classes]}, "all":
  [classes]}`: for each planet and moon, the classes it could take where
  it is without moving any other planet (valid in its zone, at a typical mass
  for the class still clear of its neighbors and able to hold its moons;
  for a moon, a class its planet can hold), most common first, its own
  class left out (a planet's own moons may still be re-spaced for its new
  size). `POST /api/planets/<id>/class` (or `/api/moons/<id>/class`)
  `{"class": "M", "force": false}` re-rolls the body as that class at
  the same orbit, keeping its name, row and moons. Without `force` only a
  recommended class is accepted (`409` otherwise); with `"force": true`
  any class is, even one that can't form where the body is. The body
  keeps the class it was given; the rest of the system is re-spaced
  around it, nothing is removed, and what still doesn't validate comes
  back in `warnings`.
- **Star (ADM.7).** `POST /api/systems/<id>/star` `{"star_type": "K2V"}`
  replaces a single star with a new one of that spectral type (class
  letter, subclass digit, Yerkes class; case doesn't matter), keeping
  its name, row and place in the galaxy. Planets, moons and belts keep
  their classes; their orbits scale with the square root of the change
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

### Wiki publishing — request body

`POST /api/systems/<id>/wiki` and `POST /api/sectors/<id>/wiki` share the
same request shape:

```json
{
  "backend": "wikijs",
  "path": "systems/voranthis-vesta"
}
```

- `backend`: `"wikijs"` or `"mediawiki"` — required. `501` if that
  backend has no `base_url`/credentials configured deployment-wide (see
  `GET /api/wiki-config` and `docs/config.md`'s `wiki.*` fields).
- `path`: the target page's path/slug. **Required for `"wikijs"`**,
  which addresses a page separately from its title (`src/wikiClient/wikijs.py`);
  accepted but ignored for `"mediawiki"`, whose title — the system's or
  sector's own name — is its address instead (`src/wikiClient/mediawiki.py`).

A system publishes its page rendered from the database at upload time
(Markdown to `wikijs`, wikitext to `mediawiki` — the same text `GET
/api/systems/<id>/text` returns). A sector's page is likewise built fresh at upload time from its own current detail (name,
edge, and a table of its systems — `routes.py`'s `_sector_wiki_content`).

On success (`201`), both return the new page's
`{"id", "path", "title", "url"}` and record `url` on the matching column
(`star_systems.wikijs_url`/`mediawiki_url`, or `sectors.wiki_url`) — see
`database-schema.md`'s "Rendered wiki text and URLs". `409` if a page
already exists at that path/title (every `wikiClient` backend is
create-only); `502` if the wiki instance rejected the credentials or
couldn't be reached.
Editing individual generated bodies (stars/planets/moons/belts) isn't
supported via this API beyond renaming them; regenerate the whole system
in place instead (see "Regenerating a system" above).

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
- Storage backend: in-memory by default (`PLANETGEN_RATELIMIT_STORAGE_URI`,
  default `memory://`) — correct for a single-process deployment (Flask's
  dev server, or `mod_wsgi`/`gunicorn`/waitress with exactly one process). **A
  multi-worker deployment needs a shared backend** (e.g. Redis:
  `PLANETGEN_RATELIMIT_STORAGE_URI=redis://host:6379/0`), since each
  worker otherwise tracks its own separate counters and the real,
  aggregate request rate can exceed the configured limit by roughly the
  worker count.
- `/api/health` has its own limit instead of the default:
  `ratelimit.pages.health` in `config.json`, default **60/minute** (so a
  monitor can poll it freely but nobody can use it to hammer the database).
  The HTML pages have per-IP limits of their own too (`ratelimit.pages`,
  see [`html-interface.md`](html-interface.md#rate-limits)); the API calls
  a page makes in-process count against neither.
- Exceeding a limit returns `429` with a `Retry-After` header and
  `X-RateLimit-*` headers (`RATELIMIT_HEADERS_ENABLED`).
- "Client IP" is `request.remote_addr`. Behind a separate reverse proxy
  (nginx, Caddy, IIS, Apache `mod_proxy`) that is the proxy's address
  unless `config.json`'s `proxy_fix` is set (see
  [`config.md`](config.md) and
  [`deployment/README.md`](deployment/README.md#behind-a-reverse-proxy-proxy_fix));
  without it every visitor shares one budget. Under Apache + `mod_wsgi`
  it is already the client's address.

## Running locally

```bash
pip install -e .[api]
python src/html/wsgi.py
```

Connects to the same MySQL database every other tool in this project
defaults to (`PLANETGEN_MYSQL_*` env vars, `config.json`'s `mysql` section,
or their built-in defaults — see `stellarObjects._db.MySQLConfig` and
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
`html/api/config.py`), there's no separate write-capable override. Point
it at an account with `INSERT`/`UPDATE`/`DELETE`/`SELECT` grants in
production — see [`deployment/README.md`](deployment/README.md#mysql-accounts).

Admin logins/sessions/API keys/audit log live in a separate **control
schema** (`PLANETGEN_CONTROL_DATABASE`, default `planetgen_control`),
global to the deployment rather than per-galaxy-database — see
`stellarObjects/control_schema.sql`'s header comment and
[`database-schema.md`](database-schema.md). `migrateDb.py` creates and
seeds it alongside its usual content-schema migration.

### The first admin login

There is no published default password. The first time `migrateDb.py`
runs against an empty control schema (so on the first `install.sh`), it
creates the user `admin` with a random password and prints both once, in
a boxed block in the installer's output. Only the password's hash is
stored, so nothing can show it again. Log in at `/login` with it; you are
sent straight to `/account` to choose your own username and password, and
every other admin page and write endpoint refuses to work until you do.

#### Resetting the admin login

If that password is lost (or every admin is locked out), delete the admin
rows and let `migrateDb.py` seed a new one:

    mysql -e 'DELETE FROM planetgen_control.admin_users'
    python3 src/migrateDb.py        # or ./update.sh

(use your `PLANETGEN_CONTROL_DATABASE` name if you changed it). Deleting
the rows also deletes their sessions and API keys; the audit log keeps
its entries. The new random password is printed as above.

## Deploying

`src/html/wsgi.py` exposes the standard WSGI `application` object: the
same app serves the API under `/api` and every HTML page. Every hosting
setup loads that one file:

- **Apache + `mod_wsgi`** (the reference setup, `examples/apache/`): a
  `WSGIScriptAlias` for `/` pointing at `src/html/wsgi.py`, in its own
  `WSGIDaemonProcess`, plus the `<Directory>` blocks that deny direct
  requests into `html/api/` and `html/lib/`. See
  [`deployment/apache.md`](deployment/apache.md).
- **gunicorn** behind nginx or Caddy (Linux, macOS):
  `gunicorn --pythonpath <checkout>/src/html wsgi:application`.
- **waitress** behind IIS, Caddy or Apache (Windows):
  `waitress-serve wsgi:application` run from `src/html`.

[`deployment/README.md`](deployment/README.md) compares them and links
each guide. Behind nginx, Caddy, IIS or Apache's `mod_proxy`, set
`config.json`'s `proxy_fix` so the rate limits see the client's address
and the app knows a request came over HTTPS; under `mod_wsgi` leave it
off. No separate vhost or `ServerName` is needed for the API: the pages
call it in-process (`html/web/transport.py`), not over HTTP.

Point the app at the deployed MySQL database with `config.json`'s `mysql`
section, or `PLANETGEN_MYSQL_*` in the process environment. That account
needs write grants, not just `SELECT`: the admin pages and every write
endpoint write through it (see "Running locally" above and
[`deployment/README.md`](deployment/README.md#mysql-accounts)). Under
`mod_wsgi`, environment variables must be in the Apache service's own
process environment (e.g. `/etc/apache2/envvars`, or an `Environment=`
line on the apache2 systemd unit), **not** the vhost's `SetEnv`
directives: those never reach `os.environ` under `mod_wsgi` --
`html/api/config.py` reads its config from `os.environ` once, at process
startup, and `SetEnv` values only ever show up in a request's `environ`
dict, which doesn't exist yet at that point. Under gunicorn or waitress,
the service's own environment (systemd `EnvironmentFile`, launchd
`EnvironmentVariables`, WinSW `<env>`) works the normal way. If running
more than one process, also set `PLANETGEN_RATELIMIT_STORAGE_URI` to a
shared backend (see "Rate limiting").

**The admin session cookie requires HTTPS.** It's set `Secure` by default
(`config.SESSION_COOKIE_SECURE`) -- the browser never sends it over plain
HTTP, so `/login` and `/admin` won't work on an HTTP-only site. Terminate
TLS in front of the app (e.g. `certbot --apache` or `--nginx`, or Caddy's
automatic HTTPS) before using the admin pages; only set
`PLANETGEN_ADMIN_COOKIE_INSECURE=1` for local development without TLS in
front, never in production.

## Not done yet

Planned API changes, as of 2026-10-07. None is built; the
items and their full text are in `docs/TODO.md`, and the phases in
`docs/plan/`.

- **Request validation and rate limits on libraries (ADM.21, SEC.30,
  phase 0).** Request bodies are Pydantic models, and limits move to
  Flask-Limiter with Redis storage; jobs started through the API run on
  RQ (PERF.24). See [`design/library-migration.md`](design/library-migration.md).
- **Every API call logged (API.15, phase 1).** By default every call is
  logged with the time, the route, the account that made it (the key's
  owner for an API key, the signed-in admin for the web site, and `god`
  for the console), how it came in (API key, web session or console),
  and the HTTP response code it got back. The log uses HTTP response
  codes for parity with the web server's own logs. It needs no user
  accounts: it logs today's admin accounts, and any account once user
  accounts exist (USR.1). The activity log
  ([`config.md`](config.md#the-activity-log)) stays as it is.
- **Routes with no hop limit (NAV.12, phase 1).** `/api/nav` always
  finds a route between two placed endpoints, a same-sector route may
  leave the sector, the answer gives the longest hop, and each hop
  carries a flag saying whether its line crosses unfilled (unknown)
  sectors. Travel times per hop and for the route follow (NAV.11). See
  [`design/course-routing.md`](design/course-routing.md).
- **Key scopes (API.9, phase 1).** API keys get a scope: read, admin or
  upload.
- **The galaxy's seed, version and run history (API.16, phase 2).** A
  route returns the galaxy's 128-bit seed, its 22-hex-digit version key
  and the run history; API.12's download uses the same fields. See
  [`design/reproducible-galaxies.md`](design/reproducible-galaxies.md).
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

Editing individual generated bodies (stars, planets, moons, belts)
beyond renaming them or regenerating the whole system is also still
open. The pages in `../src/html/web/` (see
[`html-interface.md`](html-interface.md)) are this API's own
server-rendered frontend.
