# planetGen API

A read-only JSON API over the planetGen database (`src/stellarObjects/schema.sql`),
built with [Flask](https://flask.palletsprojects.com/). This is `TODO.md`'s
Phase 5 backend — no write path (generation still happens through
`sectorGen.py`/`systemGen.py`, which persist directly), no frontend yet, and
now backed by the same MySQL database every other tool in this project uses
(see [`database-schema.md`](database-schema.md) for the MySQL port).

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

All under `/api/`, all read-only, all JSON:

- `GET /api/health` — liveness/readiness check: confirms the process is up
  and the configured database can actually be opened. Returns
  `{"status": "ok"}`, or `{"status": "error", "detail": "..."}` with a `503`
  if the database can't be reached.
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
(including an out-of-range `limit`/`offset`, or a non-numeric/non-positive
`radius`) is a `400`, an unmatched URL is a `404`, an unsupported HTTP method
is a `405`, and an unexpected server-side failure is a `500` — the API never
falls through to Flask's default HTML error page or leaks a stack trace to
the client (the real detail still reaches Flask's own logger).

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

This API never writes — point `PLANETGEN_MYSQL_USER`/`PLANETGEN_MYSQL_PASSWORD`
at a database account with `SELECT`-only grants in production, same
recommendation as `queryDb.py`'s (see that script's module docstring).

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
read-only account (see above).

## Not done yet

See `TODO.md`'s Phase 5 section — a frontend is the remaining open item.
