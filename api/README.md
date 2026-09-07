# planetGen API

A read-only JSON API over the planetGen database (`stellarObjects/schema.sql`),
built with [Flask](https://flask.palletsprojects.com/). This is the start of
`TODO.md`'s Phase 5 backend — groundwork only: no write path (generation
still happens through `sectorGen.py`/`systemGen.py`, which persist directly),
no frontend yet, and still backed by the same SQLite database every other
tool in this project uses today.

## Why Flask

Comparison against FastAPI/Django REST Framework: the existing persistence
layer (`stellarObjects/_db.py`) is deliberately raw `sqlite3`/plain SQL with
no ORM, this API is read-heavy with no concurrency pressure yet, and it
needs to deploy onto the same Apache2/VPS setup that already serves the
interim `html/` CGI browser (see `apache/README.md`). Flask has no opinion
about the data layer (route handlers call straight into `queryDb.py`'s and
`stellarObjects._db`'s existing functions), deploys via `mod_wsgi` in the
same Apache process model the CGI scripts already use, and can be mounted
at `/api/` alongside `html/` for an incremental rollout rather than a hard
cutover. FastAPI's headline advantages (async, auto-generated OpenAPI docs)
don't pay for themselves yet: `sqlite3`'s driver is synchronous regardless
of framework, and there's no separate frontend consuming this API yet to
benefit from generated docs. Worth revisiting if/when a dedicated frontend
(Phase 5's other open item) makes API-contract docs valuable.

## Endpoints

All under `/api/`, all read-only, all JSON:

- `GET /api/sectors` — every sector, with system count (`queryDb.list_sectors`).
- `GET /api/sectors/<id>` — one sector's full nested detail
  (`stellarObjects._db.load_sector(...).to_dict()`).
- `GET /api/systems?star_type=<prefix>&sector_id=<id>` — filtered system
  listing (`queryDb.list_systems`).
- `GET /api/systems/<id>` — one system's full nested detail (stars, planets,
  moons, belts) (`stellarObjects._db.load_star_system(...).to_dict()`).
- `GET /api/systems/<id>/near?radius=<ly>` — other systems in the same
  sector within `radius` light-years (`queryDb.systems_within_radius`).

A missing sector/system id returns `404` with `{"error": "..."}`; a missing
required query parameter returns `400`.

## Running locally

```bash
pip install -e .[api]
python -m api.app
```

Defaults to `db/planetgen.db` (same default as every other tool). Point it
at a different database with:

```bash
PLANETGEN_DB_PATH=/path/to/other.db python -m api.app
```

## Deploying behind Apache (mod_wsgi)

`wsgi.py` (repo root) exposes the standard `application` object `mod_wsgi`
expects. Add a `WSGIScriptAlias` for `/api` pointing at `wsgi.py` to the
existing vhost config in `apache/` (see `apache/README.md` for the vhost
this project already deploys, including `set-permissions.sh`), or run it
behind `gunicorn` + `mod_proxy`/`mod_proxy_http` if `mod_wsgi` isn't
available. Either way, set `PLANETGEN_DB_PATH` in the process environment
(e.g. the vhost's `SetEnv` directive, or the `gunicorn` service's
environment file) to point at the deployed `db/planetgen.db`.

## Not done yet

See `TODO.md`'s Phase 5 section — a frontend, and the eventual move off
SQLite to MySQL (the database server already present on the deployment
host) for real concurrent multi-user access. That migration is a
significantly larger effort than this scaffold and is deliberately not
started here; see `TODO.md` for the concrete list of what it will require.
