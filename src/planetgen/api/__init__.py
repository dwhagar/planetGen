# planetgen/api/__init__.py

"""
The JSON API over the planetGen database (served under `/api`).

The reads are `routes.py` (sectors, systems, planets and moons, display
detail with `planetgen.db.query`'s `sector_detail` and `system_detail`, a
faceted `search`, the Galaxy Map's tiles, stages and windows, navigation
courses), `population.py` (species, polities and territories) and
`facility_types.py` (facility types and their fields). The writes are
`edits.py` (an admin's edits to stored systems), `workqueue.py` (the work
queue page's controls and queued generation jobs), `auth.py` (sign-in, accounts and API keys, with
`authz.py` resolving the caller from a session cookie or a Bearer key) and
`admin.py` (server health and database statistics for the admin pages).
Every route is rate limited by Flask-Limiter (`limiter.py`, with
`loginguard.py` around each password check), validated with Pydantic
(`schemas.py`), versioned (`version.py`), and answers object IDs as printed
IDs (`ids.py`, `idpaths.py`, `naming.py`). `common.py` holds what the route
modules share and `config.py` the API's configuration. `docs/api.md` lists every
route, its scopes and its errors.

`planetgen.db.store`'s `load_sector` and `load_star_system` are different:
they rebuild the *generation* object graph, not the display shape these
routes return.

The app itself is `planetgen/web/app.py`'s `create_app`, which serves these
routes together with the HTML pages (`planetgen/web/`); `src/html/wsgi.py` is
the mod_wsgi and gunicorn entry point. The pages fetch their data from the
routes in process (`planetgen/web/lib/apiclient.py`,
`planetgen/web/transport.py`), not over HTTP and not from MySQL directly. See
`docs/api.md` for how to run it and `docs/html-interface.md` for the pages.
"""
