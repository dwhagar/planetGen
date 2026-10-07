# planetgen/api/__init__.py

"""
JSON API over the planetGen database (Phase 5, TODO.md).

A Flask app exposing the same read queries `planetgen.db.query` already provides
as a CLI, plus full sector/system display detail (`queryDb.sector_detail`/
`system_detail` -- distinct from `planetgen.db.store`'s `load_sector`/
`load_star_system`, which reconstruct the *generation* object graph
instead) and a faceted search (`queryDb.search`), rate-limited via
Flask-Limiter. Write endpoints (create/modify/delete a sector or system)
exist as stubs -- routed and validated, but not yet wired up to the
database; see `docs/api.md`'s "Write endpoints" section for what's still
needed before they do anything real.

Lives under `html/` (moved here from `src/api/`) next to the HTML pages
(`planetgen/web/`), which the same Flask app serves -- see `html/wsgi.py` for
the mod_wsgi entry point that imports this package, and
`examples/apache/planetgen.conf.example` for the vhost. The pages fetch
their data from these routes in-process (`planetgen/web/lib/apiclient.py`,
`planetgen/web/transport.py`) instead of querying MySQL directly -- see
`docs/api.md` for how to run this and `docs/html-interface.md` for that
side of it. `docs/TODO.md`'s "Open items" section covers what's still
open (the write endpoints above, mainly).
"""
