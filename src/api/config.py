# api/config.py

"""
Configuration for the planetGen API.

`MYSQL_CONFIG` is a `stellarObjects._db.MySQLConfig`, itself built from
the same `PLANETGEN_MYSQL_*` environment variables every other entry
point in this project (`sectorGen.py`, `systemGen.py`, `queryDb.py`)
reads, so a WSGI deployment (see `wsgi.py`) points this API at a specific
database without editing code, set via the vhost's `SetEnv` directives or
the `gunicorn` service's environment file.

Every read in `routes.py` goes through this same connection today, and
this project's docs recommend a `SELECT`-only database account for it
(see `queryDb.py`'s module docstring) -- true for every route currently
implemented, including the write stub endpoints (`POST`/`PATCH`/`DELETE`
on `/sectors`/`/systems`), which validate their request body but don't
touch the database yet (see `routes.py`). **Once those stubs are filled
in with real inserts/updates/deletes, `MYSQL_CONFIG` will need a
write-capable account instead** -- also revisit authentication/
authorization at that point: nothing in this module or `routes.py` checks
who's calling, which is fine for a read-only API but not once a request
can mutate data.

`RATELIMIT_DEFAULT`/`RATELIMIT_STORAGE_URI` configure Flask-Limiter (see
`limiter.py`) -- applied to every route by default; `routes.py`'s write
endpoints additionally layer a stricter per-route limit on top. The
in-memory storage default suits a single-process deployment (Flask's own
dev server, or `mod_wsgi`/`gunicorn` running exactly one worker); a
multi-worker deployment needs a shared backend (e.g. Redis, via
`PLANETGEN_RATELIMIT_STORAGE_URI=redis://...`) since each worker would
otherwise track its own separate request counts, letting the true
request rate exceed the configured limit by roughly a factor of the
worker count.
"""

import os

from stellarObjects._db import MySQLConfig

DEFAULT_RATE_LIMITS = "200 per day;50 per hour"
"""str: Flask-Limiter's own quickstart uses this exact pair as its
"standard" example limit -- a reasonable default for a low-traffic public
API with no other guidance, and easy to override per-deployment via
`PLANETGEN_RATELIMIT_DEFAULT` without editing code. A single
semicolon-separated *string*, not a list -- Flask-Limiter's own
`RATELIMIT_DEFAULT` config value must be a string (or a callable), and
silently mis-parses (raising `ValueError` on every request, confirmed by
testing) if handed a list of strings instead, however natural that looks
in Python."""


class Config:
    MYSQL_CONFIG = MySQLConfig()

    RATELIMIT_DEFAULT = os.environ.get("PLANETGEN_RATELIMIT_DEFAULT", DEFAULT_RATE_LIMITS)
    RATELIMIT_STORAGE_URI = os.environ.get("PLANETGEN_RATELIMIT_STORAGE_URI", "memory://")
    RATELIMIT_HEADERS_ENABLED = True
