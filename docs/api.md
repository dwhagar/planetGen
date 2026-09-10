# planetGen API

A JSON API over the planetGen database (`src/stellarObjects/schema.sql`),
built with [Flask](https://flask.palletsprojects.com/). This is `TODO.md`'s
Phase 5 backend — backed by the same MySQL database every other tool in this
project uses (see [`database-schema.md`](database-schema.md) for the MySQL
port). It now has a frontend: the interim `../src/html/` browser
([`html-interface.md`](html-interface.md)) is a thin client over this API
rather than a direct database consumer — see that doc's own note on the
switch, and "Deploying behind Apache" below for how both are mounted on
one vhost.

Every read endpoint is fully implemented. The write endpoints (create/modify/
delete a sector or system) are **stubs**: routed, rate-limited, and
validating their request body against the schema below, but not yet wired
up to the database — see "Write endpoints" below before building anything
against them.

## Why Flask

Comparison against FastAPI/Django REST Framework: the persistence layer
(`stellarObjects/_db.py`) is deliberately plain SQL over a small `pymysql`
wrapper with no ORM, this API is read-heavy with no concurrency pressure yet,
and it needs to deploy onto the same Apache2/VPS setup that already serves the
interim `../src/html/` CGI browser (see [`apache-deployment.md`](apache-deployment.md)).
Flask has no opinion
about the data layer (route handlers call straight into `queryDb.py`'s and
`stellarObjects._db`'s existing functions), deploys via `mod_wsgi` in the
same Apache process model the CGI scripts already use, and lives at
`../src/html/api/` -- served from the same tree/DocumentRoot as the CGI
browser rather than a separately-deployed package -- mounted at `/api/`
alongside `../src/html/` on one vhost (see "Deploying behind Apache"
below). `../src/html/` is now this API's own frontend (see that
directory's docs), not a separate database consumer anymore. FastAPI's
headline advantages (async, auto-generated OpenAPI docs) still don't pay
for themselves: this API is read-heavy and low-concurrency regardless of
framework, and `../src/html/` is a plain server-rendered CGI client with
no use for generated API-contract docs the way a JS single-page app
would. Worth revisiting if a richer JS frontend is ever built against
this API instead.

## Endpoints

All under `/api/`, all JSON in, JSON out. Every endpoint below except
`/health` and `/databases` accepts an optional `db=<name>` query
parameter selecting *which* MySQL schema on the configured server to
read from (validated against the same prefix-filtered list
`/api/databases` itself returns — an unrecognized name is a `404`, same
as an unrecognized sector/system id); omitted, it falls back to
`MYSQL_CONFIG`'s own configured default database (`config.py`). This is
what lets one API process back `../src/html/`'s multi-database picker
(`index.py`'s `?db=`) — see `stellarObjects._db.list_databases`/
`resolve_database`.

### Read

- `GET /api/health` — liveness/readiness check: confirms the process is up
  and the configured database can actually be opened. Returns
  `{"status": "ok"}`, or `{"status": "error", "detail": "..."}` with a `503`
  if the database can't be reached. Exempt from rate limiting.
- `GET /api/databases` — every MySQL schema on the configured server
  matching this deployment's prefix (`stellarObjects._db.list_databases`),
  each with `name`, `size_bytes`, `modified_at`, and a quick-glance
  `sector_count`/`system_count` (`null` for a matching schema missing this
  project's own tables, e.g. mid-migration, rather than failing the whole
  listing). `../src/html/index.py`'s database picker.
- `GET /api/sectors?limit=<n>&offset=<n>` — every sector, paginated (see
  "Pagination" below), each with `id`, `name`, `edge_mpc`, `edge_ly`,
  `system_count`, and its galaxy placement (`center_x_pc`/`center_y_pc`/
  `center_z_pc`/`shell_index`/`shell_slot_index`, all `null` together for
  an unplaced sector, plus a `placed` bool) (`queryDb.list_sectors`/
  `count_sectors`).
- `GET /api/sectors/<id>` — one sector's full display detail: the same
  fields as the listing above, plus `systems` (every system placed in
  it — `id`, `name`, `quadrant`, `location`, `is_binary`, `binary_type`,
  `position_x_mpc`/`position_y_mpc`/`position_z_mpc`, and `stars`, each
  with `role`/`star_type`/`temperature_k`/`radius_km`/`luminosity_w`)
  (`queryDb.sector_detail`). Distinct from
  `stellarObjects._db.load_sector(...).to_dict()`'s *generation* object
  graph (config/provenance, no database ids) — this is the flat,
  ids-and-display-fields shape `../src/html/sector.py`'s systems table
  and Sector Map actually need.
- `GET /api/systems?star_type=<prefix>&sector_id=<id|none>&limit=<n>&offset=<n>` —
  filtered, paginated system listing (`queryDb.list_systems`/
  `count_systems`), each with `id`, `name`, `sector_id`, `is_binary`, and
  `star_summary` (the single star's `star_type`, or a binary's
  `binary_type`). `sector_id=none` matches only standalone systems
  (`sector_id IS NULL`, `../src/html/browse.py`'s own table) — distinct
  from omitting `sector_id` entirely (no sector filter at all).
- `GET /api/systems/<id>` — one system's full display detail: `id`,
  `name`, `sector_id`, `quadrant`, `location`, `is_binary`, `binary_type`,
  `markdown_content`, `wikitext_content`, `stars`, `planets` (each with
  its own nested `moons`), `belts`, and `sector_siblings` (`{id, name}`
  for every other system in the same sector, for linkifying `location`'s
  "nearest: ..." names) (`queryDb.system_detail`) — same "flat display
  shape, not the generation object graph" relationship to
  `stellarObjects._db.load_star_system(...).to_dict()` as `/api/sectors/<id>`
  above.
- `GET /api/systems/<id>/near?radius=<ly>` — other systems in the same
  sector within `radius` light-years (`queryDb.systems_within_radius`).
- `GET /api/nav?from=<id>&to=<id>` — course, distance, and an optimal route
  between two systems (`queryDb.nav_between`) — see "NAV" below.
- `GET /api/galaxy/sectors` — every galaxy-placed sector (non-`null`
  galaxy placement), each with `id`, `name`, `x`/`y`/`z`
  (`center_x/y/z_pc`), `galactic_radius_pc`, `shell_index`, and
  `system_count` (`queryDb.galaxy_placed_sectors`) — the data
  `../src/html/galaxy.py`'s Galaxy Map plots. Not paginated: bounded by
  how much of the galaxy has actually been generated (see `TODO.md`'s
  Phase 4 lazy-generation design), not by the addressable galaxy's own
  scale.
- `GET /api/search?sector_q=&system_q=&star_q=&planet_q=&moon_q=&<facet>=<value>...` —
  the faceted search behind `../src/html/search.py`: click-to-filter tags
  (object type; star spectral/luminosity class; planet/moon class, body
  type, supported life chemistry; asteroid belt density — repeat a facet
  name for multiple active values, e.g. `class=M&class=K`) plus a
  per-entity name search. Returns `facets` (one `{value, label, count,
  tooltip}` list per facet, built from the distinct values actually
  present), `autocomplete` (`sectors`/`systems`/`stars`/`planets`/`moons`
  name lists), `facet_labels` (`"facet:value"` -> label, for an
  active-filter chip), and `results` (`sectors`/`systems`/`stars`/
  `planets`/`moons`/`belts` -> `{"rows": [...], "truncated": bool}`, or
  `null` for an object type with no active reason to query it — see
  `queryDb.search`'s docstring for the exact inclusion rule).

### Write (stubs — see below)

- `POST /api/sectors` — create a sector.
- `PATCH /api/sectors/<id>` — modify a sector.
- `DELETE /api/sectors/<id>` — remove a sector.
- `POST /api/systems` — create a system.
- `PATCH /api/systems/<id>` — modify a system.
- `DELETE /api/systems/<id>` — remove a system.

### Pagination

`/api/sectors` and `/api/systems` return a paginated envelope rather than a
bare list — this project's own roadmap (`docs/TODO.md`, Phase 4) plans
galaxy-scale generation, so an unbounded listing endpoint would eventually
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
defaults to 0. Both must be non-negative integers when given at all, or the
request gets a `400`.

### Errors

Every error response is JSON, `{"error": "..."}`, regardless of what raised
it: a missing sector/system id is a `404`, a missing/invalid query parameter
or request body (including an out-of-range `limit`/`offset`, a
non-numeric/non-positive `radius`, or a write endpoint's body failing
validation — see below) is a `400`, an unmatched URL is a `404`, an
unsupported HTTP method is a `405`, exceeding a rate limit is a `429`
(`{"error": "rate limit exceeded", "detail": "..."}`, see "Rate limiting"),
and an unexpected server-side failure is a `500` — the API never falls
through to Flask's default HTML error page or leaks a stack trace to the
client (the real detail still reaches Flask's own logger).

## NAV

`GET /api/nav?from=<system_id>&to=<system_id>` returns a direct course
(distance, azimuth, altitude, warp travel times) plus an optimal route via
adjacent systems (`stellarObjects.navGraph`, a k-nearest-neighbor adjacency
graph with Dijkstra shortest-path) between two systems.

**Availability.** NAV only applies to a pair of systems that satisfy all of:

- Both systems are assigned to a sector (`star_systems.sector_id IS NOT
  NULL`) — a `404` (unknown system id) or a `400` with `"NAV requires both
  systems to be assigned to a sector"` otherwise.
- If the two systems are in different sectors, both of those sectors must
  have a galaxy placement (`sectors.center_x/y/z_pc IS NOT NULL`, i.e. both
  were placed by `galaxyGen.py`) — otherwise a `400` with `"NAV between
  different sectors requires both sectors to have a galaxy placement"`.
  Same-sector NAV never needs this — it works even in a database with no
  galaxy generated at all.

**Course convention.** Azimuth and altitude are both galactic-plane-relative
(see `stellarObjects/navigation.py`'s module docstring): azimuth is the
angle in the galactic X-Y plane measured counterclockwise from +X (0-360°),
altitude is elevation above (+) or below (-) that plane (-90° to +90°) —
not a bearing relative to any particular ship heading.

**Response:**

```json
{
  "scope": "sector",
  "direct": {
    "distance_ly": 4.0,
    "azimuth_deg": 0.0,
    "altitude_deg": 0.0
  },
  "warp_times": [
    {"warp_factor": 1, "velocity_multiple_of_c": 1.0, "years": 4.0, "formatted": "4 years and 1 day"},
    {"warp_factor": 3, "velocity_multiple_of_c": 38.94, "years": 0.1, "formatted": "37 days 12 hours and 27 minutes"},
    {"warp_factor": 6, "velocity_multiple_of_c": 392.5, "years": 0.01, "formatted": "3 days 17 hours and 20 minutes"},
    {"warp_factor": 9, "velocity_multiple_of_c": 1516.38, "years": 0.003, "formatted": "23 hours and 7 minutes"}
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
- `warp_times`: travel time for `direct.distance_ly` at warp 1, 3, 6, and 9
  (`velocity_multiple_of_c = warp_factor ** (10/3)`), formatted via the same
  duration formatter used elsewhere in this project.
- `origin_position`/`destination_position`: the `[x, y, z]` light-year
  positions `direct` was computed from, in `scope`'s frame (sector-local for
  `"sector"`, absolute galaxy-frame for `"galaxy"`) — what `html/nav.py`'s
  NAV Map plot (`html/lib/navmap.py`) draws.
- `route`: the shortest path via adjacent systems (nodes: every system in
  scope; edges: each system's `k`-nearest neighbors, symmetrized), as a list
  of system ids from `from` to `to` inclusive, its total distance, and
  `positions` (one `[x, y, z]` entry per id in `path`, same frame as
  `origin_position`/`destination_position`). `null` if no path exists
  through the adjacency graph (only possible for the `"galaxy"` scope — the
  `"sector"` scope's graph is always fully reachable since every system in a
  sector gets an edge once `k` is at least the sector's own system count
  minus one).

## Write endpoints

**These are stubs.** Every one validates its request body against the
schema below and applies the rate limit, but always responds
`501 {"error": "... is not implemented yet"}` — no row is ever inserted,
updated, or deleted. They exist now so the request/response contract is
settled and something can already build/test against the validation and
error shape before the actual database logic lands.

Before filling them in for real, two things this stub stage deliberately
doesn't need yet still have to land first:

- **A write-capable database account.** `MYSQL_CONFIG` (see `config.py`)
  is the same connection every read endpoint uses, and this project's own
  docs recommend a `SELECT`-only account for it. A real write needs a
  second, write-capable account/config — reusing the read-only one will
  simply fail with a permissions error.
- **Authentication/authorization.** Nothing in this API currently checks
  *who* is calling — fine when every route only reads, not once a request
  can create/modify/delete data. Decide on an auth scheme (API key, JWT,
  mTLS, ...) before wiring real logic behind these routes.

### Sectors — request body

`POST /api/sectors` (all fields required) and `PATCH /api/sectors/<id>`
(any non-empty subset) both take:

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

**Not yet decided.** `POST`/`PATCH /api/systems` today only check that the
body is a JSON object — a system is a much richer nested object (stars,
planets, moons, belts) than a sector, and this project hasn't settled
whether creating one via the API should take a generation "recipe" (shaped
like `SystemConfig` — `star_type`, `planets`, `moons`, ... — closer to what
`systemGen.py --star-type ...` takes and letting the server generate the
system), a fully-specified object graph (shaped like
`StarSystem.to_dict()`, with every star/planet/moon/belt spelled out), or
both. Settle this alongside filling in the stub itself, and document the
chosen shape here.

## Rate limiting

Every route (`/api/health` excepted) is rate-limited via
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
  dev server, or `mod_wsgi`/`gunicorn` with exactly one worker). **A
  multi-worker deployment needs a shared backend** (e.g. Redis:
  `PLANETGEN_RATELIMIT_STORAGE_URI=redis://host:6379/0`), since each
  worker otherwise tracks its own separate counters and the real,
  aggregate request rate can exceed the configured limit by roughly the
  worker count.
- Exceeding a limit returns `429` with a `Retry-After` header and
  `X-RateLimit-*` headers (`RATELIMIT_HEADERS_ENABLED`).

## Running locally

```bash
pip install -e .[api]
python src/html/wsgi.py
```

Connects to the same MySQL database every other tool in this project
defaults to (`PLANETGEN_MYSQL_*` env vars, or their built-in defaults — see
`stellarObjects._db.MySQLConfig`). Point it at a different database with:

```bash
PLANETGEN_MYSQL_HOST=db.example.com PLANETGEN_MYSQL_DATABASE=planetgen_alpha python src/html/wsgi.py
```

Every read goes through `PLANETGEN_MYSQL_USER`/`PLANETGEN_MYSQL_PASSWORD` —
point those at a database account with `SELECT`-only grants in production,
same recommendation as `queryDb.py`'s (see that script's module docstring).
This is safe today even with the write endpoints present, since they're
stubs that never touch the database — revisit once they're implemented for
real (see "Write endpoints" above).

## Deploying behind Apache (mod_wsgi)

`src/html/wsgi.py` exposes the standard `application` object `mod_wsgi`
expects, and now lives inside the same `html/` tree the vhost's
`DocumentRoot` already points at (see "Why Flask" above -- moved there
from `src/wsgi.py`/`src/api/` so the API is served from the same
checkout/deployment tree as the CGI browser instead of a second,
separately-tracked location). Add a `WSGIScriptAlias` for `/api` pointing
at `src/html/wsgi.py` to the existing vhost config in `examples/apache/`
(see [`apache-deployment.md`](apache-deployment.md) for the vhost this
project already deploys, including `set-permissions.sh` and the
`<Directory>` block that denies direct requests into `html/api/` the same
way it already does for `html/lib/`), or run it behind `gunicorn` +
`mod_proxy`/`mod_proxy_http` if `mod_wsgi` isn't available. No separate
vhost/`ServerName` is needed either way: `../src/html/`'s own pages read
`PLANETGEN_API_BASE_URL` (default `http://127.0.0.1/api`, i.e. this same
vhost) to find the API -- point it at wherever `gunicorn` ends up
listening if using that instead. Either way, set `PLANETGEN_MYSQL_*` in
the process environment to point at the deployed MySQL database, ideally
via a read-only account (see above) -- under `mod_wsgi` this means the
Apache service's own process environment (e.g. `/etc/apache2/envvars`, or
an `Environment=` line on the apache2 systemd unit), **not** the vhost's
`SetEnv` directives: those work for the CGI browser (`mod_cgi`/`mod_cgid`
copies them into each script's real process environment) but never reach
`os.environ` under `mod_wsgi` -- `html/api/config.py` reads its config
from `os.environ` once, at process startup, and `SetEnv` values only ever
show up in a request's `environ` dict, which doesn't exist yet at that
point. Under `gunicorn`, its own service's environment file works the
normal way. If running more than one `mod_wsgi`/`gunicorn` worker, also set
`PLANETGEN_RATELIMIT_STORAGE_URI` to a shared backend (see "Rate
limiting").

## Not done yet

See `TODO.md`'s Phase 5 section — everything under "Write endpoints"
above (the real insert/update/delete logic, a write-capable database
account, and authentication/authorization). The frontend gap that
section used to describe is closed: `../src/html/` (see
[`html-interface.md`](html-interface.md)) is this API's own server-
rendered frontend now.
