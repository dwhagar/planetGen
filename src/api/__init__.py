# api/__init__.py

"""
JSON API over the planetGen database (Phase 5, TODO.md).

A Flask app exposing the same read queries `queryDb.py` already provides
as a CLI, plus full system/sector detail via `stellarObjects._db`'s
existing `load_star_system`/`load_sector` read path, rate-limited via
Flask-Limiter. Write endpoints (create/modify/delete a sector or system)
exist as stubs -- routed and validated, but not yet wired up to the
database; see `docs/api.md`'s "Write endpoints" section for what's still
needed before they do anything real.

See `docs/api.md` for how to run this and `docs/TODO.md`'s Phase 5
section for what's still open (a frontend is the remaining major item;
the MySQL migration this package's config once called "eventual" is
done).
"""

# TODO: whether this package should eventually move into `../src/html/` (so
# the API is served from the same tree/DocumentRoot as the interim browser)
# is still an open question, not a decision -- see docs/TODO.md,
# "Investigate Further" for the full reasoning already written up there.
