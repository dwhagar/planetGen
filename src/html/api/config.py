# html/api/config.py

"""
Configuration for the planetGen API.

`MYSQL_CONFIG` is a `stellarObjects._db.MySQLConfig`, itself built from
the same `PLANETGEN_MYSQL_*` environment variables (or `config.json`'s
`mysql` section -- see `stellarObjects.appconfig` and `docs/config.md`)
every other entry point in this project (`sectorGen.py`, `systemGen.py`,
`queryDb.py`) reads, so a WSGI deployment (see `wsgi.py`) points this API
at a specific database without editing code, set via the vhost's `SetEnv`
directives, the `gunicorn` service's environment file, or a single shared
`config.json`.

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

from stellarObjects._db import MySQLConfig, control_mysql_config
from stellarObjects.appconfig import load_config

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


def _write_mysql_config(read_only, write_defaults):
    """
    Builds the write-capable `MySQLConfig` from `PLANETGEN_MYSQL_WRITE_*`,
    falling back field-by-field to `config.json`'s `mysql_write` section
    (`write_defaults`) and then to `read_only`'s own values (the existing
    `SELECT`-only `PLANETGEN_MYSQL_*` config) when a given field isn't set
    anywhere -- so a single-account local/dev setup (one set of
    credentials for everything) keeps working with no extra
    configuration, while a production deployment is documented
    (`docs/apache-deployment.md`) to point `PLANETGEN_MYSQL_WRITE_*` (or
    `config.json`'s `mysql_write`) at a distinct account with `INSERT`/
    `UPDATE`/`DELETE` (but not `CREATE`/`DROP`) grants on both the content
    schemas and the control schema (`stellarObjects._db.control_mysql_config`),
    instead of reusing the read-only one -- reusing it would simply fail
    every write with a permissions error.

    Args:
        read_only (MySQLConfig): The existing `PLANETGEN_MYSQL_*` config,
            used as the final fallback for any field left unset.
        write_defaults (dict): `config.json`'s `mysql_write` section --
            empty-string fields mean "inherit from `read_only`".

    Returns:
        MySQLConfig: Ready to pass to `stellarObjects._db.open_write`.
    """
    return MySQLConfig(
        host=os.environ.get("PLANETGEN_MYSQL_WRITE_HOST") or write_defaults["host"] or read_only.host,
        port=int(os.environ.get("PLANETGEN_MYSQL_WRITE_PORT") or write_defaults["port"] or read_only.port),
        user=os.environ.get("PLANETGEN_MYSQL_WRITE_USER") or write_defaults["user"] or read_only.user,
        password=os.environ.get("PLANETGEN_MYSQL_WRITE_PASSWORD") or write_defaults["password"] or read_only.password,
        database=os.environ.get("PLANETGEN_MYSQL_WRITE_DATABASE") or write_defaults["database"] or read_only.database,
    )


def _session_cookie_secure(admin_cookie_insecure):
    """
    `PLANETGEN_ADMIN_COOKIE_INSECURE` (if set) wins outright; otherwise
    `config.json`'s `admin_cookie_insecure` decides. See `Config.SESSION_COOKIE_SECURE`'s
    own comment for why this defaults to secure.
    """
    env_value = os.environ.get("PLANETGEN_ADMIN_COOKIE_INSECURE")
    if env_value is not None:
        return env_value != "1"
    return not admin_cookie_insecure


_config_file = load_config()


class Config:
    MYSQL_CONFIG = MySQLConfig()

    # Write-capable config for the (now real, no longer stub) sector/
    # system write endpoints -- see `_write_mysql_config` above and
    # `routes.py`'s write handlers. `CONTROL_MYSQL_CONFIG` reuses this
    # same account's host/user/password against the separate control
    # schema (`control_schema.sql`'s header comment) that holds admin
    # logins/sessions/API keys/audit log -- see `auth.py`.
    WRITE_MYSQL_CONFIG = _write_mysql_config(MYSQL_CONFIG, _config_file["mysql_write"])
    CONTROL_MYSQL_CONFIG = control_mysql_config(WRITE_MYSQL_CONFIG)

    RATELIMIT_DEFAULT = os.environ.get("PLANETGEN_RATELIMIT_DEFAULT", _config_file["ratelimit"]["default"])
    RATELIMIT_STORAGE_URI = os.environ.get("PLANETGEN_RATELIMIT_STORAGE_URI", _config_file["ratelimit"]["storage_uri"])
    RATELIMIT_HEADERS_ENABLED = True

    # The admin session cookie (auth.py) is Secure by default -- never
    # sent over plain HTTP -- matching this project's documented
    # deployment (Apache + Let's Encrypt, docs/apache-deployment.md).
    # Only ever disable this for local development over plain HTTP (e.g.
    # `python src/html/wsgi.py` without TLS in front of it); a production
    # deployment must never set PLANETGEN_ADMIN_COOKIE_INSECURE or
    # config.json's admin_cookie_insecure.
    SESSION_COOKIE_SECURE = _session_cookie_secure(_config_file["admin_cookie_insecure"])
