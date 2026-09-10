# planetGen API

A JSON API over the planetGen database (`src/stellarObjects/schema.sql`),
built with [Flask](https://flask.palletsprojects.com/). This is `TODO.md`'s
Phase 5 backend — backed by the same MySQL database every other tool in this
project uses (see [`database-schema.md`](database-schema.md) for the MySQL
port), no frontend yet.

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
same Apache process model the CGI scripts already use, and can be mounted
at `/api/` alongside `../src/html/` for an incremental rollout rather than a hard
cutover. FastAPI's headline advantages (async, auto-generated OpenAPI docs)
don't pay for themselves yet: this API is read-heavy and low-concurrency
regardless of framework, and there's no separate frontend consuming this API yet to
benefit from generated docs. Worth revisiting if/when a dedicated frontend
(Phase 5's other open item) makes API-contract docs valuable.

## Endpoints

All under `/api/`, all JSON in, JSON out.

### Read

- `GET /api/health` — liveness/readiness check: confirms the process is up
  and the configured database can actually be opened. Returns
  `{"status": "ok"}`, or `{"status": "error", "detail": "..."}` with a `503`
  if the database can't be reached. Exempt from rate limiting.
- `GET /api/sectors?limit=<n>&offset=<n>` — every sector, with system count
  (`queryDb.list_sectors`/`count_sectors`), paginated (see "Pagination"
  below).
- `GET /api/sectors/<id>` — one sector's full nested detail
  (`stellarObjects._db.load_sector(...).to_dict()`).
- `GET /api/systems?star_type=<prefix>&sector_id=<id>&limit=<n>&offset=<n>` —
  filtered, paginated system listing (`queryDb.list_systems`/`count_systems`).
- `GET /api/systems/<id>` — one system's full nested detail (stars, planets,
  moons, belts) (`stellarObjects._db.load_star_system(...).to_dict()`).
- `GET /api/systems/<id>/near?radius=<ly>` — other systems in the same
  sector within `radius` light-years (`queryDb.systems_within_radius`).
- `GET /api/nav?from=<id>&to=<id>` — course, distance, and an optimal route
  between two systems (`queryDb.nav_between`) — see "NAV" below.

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
python src/wsgi.py
```

Connects to the same MySQL database every other tool in this project
defaults to (`PLANETGEN_MYSQL_*` env vars, or their built-in defaults — see
`stellarObjects._db.MySQLConfig`). Point it at a different database with:

```bash
PLANETGEN_MYSQL_HOST=db.example.com PLANETGEN_MYSQL_DATABASE=planetgen_alpha python src/wsgi.py
```

Every read goes through `PLANETGEN_MYSQL_USER`/`PLANETGEN_MYSQL_PASSWORD` —
point those at a database account with `SELECT`-only grants in production,
same recommendation as `queryDb.py`'s (see that script's module docstring).
This is safe today even with the write endpoints present, since they're
stubs that never touch the database — revisit once they're implemented for
real (see "Write endpoints" above).

## Deploying behind Apache (mod_wsgi)

`src/wsgi.py` exposes the standard `application` object `mod_wsgi`
expects. Add a `WSGIScriptAlias` for `/api` pointing at `src/wsgi.py` to
the existing vhost config in `examples/apache/` (see
[`apache-deployment.md`](apache-deployment.md) for the vhost this project
already deploys, including `set-permissions.sh`), or run it
behind `gunicorn` + `mod_proxy`/`mod_proxy_http` if `mod_wsgi` isn't
available. Either way, set `PLANETGEN_MYSQL_*` in the process environment
(e.g. the vhost's `SetEnv` directives, or the `gunicorn` service's
environment file) to point at the deployed MySQL database, ideally via a
read-only account (see above). If running more than one `mod_wsgi`/
`gunicorn` worker, also set `PLANETGEN_RATELIMIT_STORAGE_URI` to a shared
backend (see "Rate limiting").

## Not done yet

See `TODO.md`'s Phase 5 section — a frontend, plus everything under
"Write endpoints" above (the real insert/update/delete logic, a
write-capable database account, and authentication/authorization).
