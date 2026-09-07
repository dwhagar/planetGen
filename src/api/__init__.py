# api/__init__.py

"""
Read-only JSON API over the planetGen database (Phase 5, TODO.md).

Groundwork only: a Flask app exposing the same read-only queries
`queryDb.py` already provides as a CLI, plus full system/sector detail via
`stellarObjects._db`'s existing `load_star_system`/`load_sector` read path.
No write path -- generation still happens through `sectorGen.py`/
`systemGen.py`, which already persist to the database directly.

See `api/README.md` for how to run this and `TODO.md`'s Phase 5 section for
what's still open (framework choice is settled -- Flask, see the README for
why -- but the frontend and the eventual MySQL migration are not).
"""
