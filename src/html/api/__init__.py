# html/api/__init__.py

"""
JSON API over the planetGen database (Phase 5, TODO.md).

A Flask app exposing the same read queries `queryDb.py` already provides
as a CLI, plus full system/sector detail via `stellarObjects._db`'s
existing `load_star_system`/`load_sector` read path, rate-limited via
Flask-Limiter. Write endpoints (create/modify/delete a sector or system)
exist as stubs -- routed and validated, but not yet wired up to the
database; see `docs/api.md`'s "Write endpoints" section for what's still
needed before they do anything real.

Lives under `html/` (moved here from `src/api/`) so it's served from the
same tree/DocumentRoot as the interim CGI browser -- see `html/wsgi.py`
for the mod_wsgi entry point that imports this package, and
`examples/apache/planetgen.conf.example` for the vhost `<Directory>`
block that denies direct requests into this package the same way it
already does for `html/lib/`.

See `docs/api.md` for how to run this and `docs/TODO.md`'s Phase 5
section for what's still open (a frontend is the remaining major item;
the MySQL migration this package's config once called "eventual" is
done).
"""
