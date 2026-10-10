# planetgen/api/config.py

"""
Configuration for the planetGen API.

`MYSQL_CONFIG` is a `planetgen.db.store.MySQLConfig`, itself built from
the same `PLANETGEN_MYSQL_*` environment variables (or `config.json`'s
`mysql` section -- see `planetgen.util.settings` and `docs/config.md`)
every other entry point in this project (`sectorGen.py`, `systemGen.py`,
`planetgen.db.query`) reads, so a WSGI deployment (see `wsgi.py`) points this API
at a specific database without editing code, set via the vhost's `SetEnv`
directives, the `gunicorn` service's environment file, or a single shared
`config.json`.

Every read and write in `routes.py` goes through this same connection --
`WRITE_MYSQL_CONFIG`/`CONTROL_MYSQL_CONFIG` both simply reuse
`MYSQL_CONFIG`, there's no separate write-capable account to configure.
Give this one account whatever grants the most demanding caller needs
(`INSERT`/`UPDATE`/`DELETE` on the content schemas and the control schema,
at minimum, once any write endpoint is reachable) -- see
`docs/deployment/README.md`'s "MySQL accounts" section.

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

import secrets

from planetgen.db.store import MySQLConfig, control_mysql_config
from planetgen.util.settings import get_settings

def _wiki_config(wiki):
    """
    The resolved wiki-publishing config from the settings model's `wiki`
    section (`planetgen.util.settings.Wiki`; the `PLANETGEN_WIKIJS_*`/
    `PLANETGEN_MEDIAWIKI_*` variables are applied there).

    Neither backend has a separate on/off flag: a backend counts as
    configured -- and so is offered as an "Upload to Wiki" target by
    `routes.py`'s `POST /api/systems/<id>/wiki`/`POST /api/sectors/<id>/wiki`
    -- purely by having a non-empty `base_url` plus every credential field
    its own `planetgen.wiki` backend requires (see `_configured` below); an
    empty `base_url` alone already means "don't offer this one", so a
    separate flag would only ever duplicate that check.

    Returns:
        dict: `{"wikijs": {"base_url", "api_token", "configured"},
            "mediawiki": {"base_url", "username", "password", "configured"}}`.
    """
    wikijs = wiki.wikijs.model_dump()
    wikijs["configured"] = bool(wikijs["base_url"] and wikijs["api_token"])
    mediawiki = wiki.mediawiki.model_dump()
    mediawiki["configured"] = bool(mediawiki["base_url"] and mediawiki["username"] and mediawiki["password"])
    return {"wikijs": wikijs, "mediawiki": mediawiki}


_settings = get_settings()


def _secret_key():
    """
    `secret_key` (`PLANETGEN_SECRET_KEY`), else a random per-process key
    (`SECRET_KEY_IS_EPHEMERAL` is then True and `create_app` logs a
    warning). Used to sign the Flask-served pages' CSRF tokens
    (`planetgen/web/csrf.py`); nothing else in the app signs anything
    with it today.
    """
    if _settings.secret_key:
        return _settings.secret_key, False
    return secrets.token_hex(32), True


_SECRET_KEY, _SECRET_KEY_IS_EPHEMERAL = _secret_key()


RATELIMIT_KEY_PREFIX = "planetgen:limiter"
"""str: What Flask-Limiter's Redis keys and the login lockouts' (SEC.30)
start with; a test sets its own to keep its counts apart."""


def ratelimit_storage_uri(settings=None):
    """
    Where Flask-Limiter and the login lockouts (`loginguard.py`) count:
    `ratelimit.storage_uri`, and when that is empty the Redis server
    (`redis.url`), SEC.30.
    """
    settings = get_settings() if settings is None else settings
    return settings.ratelimit.storage_uri or settings.redis.url


class Config:
    MYSQL_CONFIG = MySQLConfig()

    # No separate write-capable account -- the sector/system write
    # endpoints (`routes.py`'s write handlers) use this same
    # `MYSQL_CONFIG`. `CONTROL_MYSQL_CONFIG` reuses its host/user/password
    # against the separate control schema (`control_schema.sql`'s header
    # comment) that holds admin logins/sessions/API keys/audit log -- see
    # `auth.py`.
    WRITE_MYSQL_CONFIG = MYSQL_CONFIG
    CONTROL_MYSQL_CONFIG = control_mysql_config(WRITE_MYSQL_CONFIG)

    RATELIMIT_DEFAULT = _settings.ratelimit.default
    RATELIMIT_STORAGE_URI = ratelimit_storage_uri(_settings)
    # A Redis outage falls back to counting in memory instead of failing requests.
    RATELIMIT_IN_MEMORY_FALLBACK_ENABLED = True
    RATELIMIT_KEY_PREFIX = RATELIMIT_KEY_PREFIX
    RATELIMIT_HEADERS_ENABLED = True
    # Per-IP limits on the HTML pages and /api/health (`web/ratelimits.py`):
    # `search`, `galaxy`, `galaxy_tiles`, `health`, and `other` for every
    # other page. An empty value turns that limit off.
    RATELIMIT_PAGES = _settings.ratelimit.pages.model_dump()
    # Per-address lockout and per-username backoff (`loginguard.py`).
    # Always on in a real deployment; only a test harness that fails many
    # logins on purpose turns it off.
    LOGIN_BACKOFF_ENABLED = True
    # Addresses and networks never locked out (loopback never is either):
    # `config.json`'s `login_allowlist`, or PLANETGEN_LOGIN_ALLOWLIST
    # (comma- or space-separated).
    LOGIN_ALLOWLIST = list(_settings.login_allowlist)

    # See `_wiki_config` above -- read by `routes.py`'s
    # `POST /api/systems/<id>/wiki`/`POST /api/sectors/<id>/wiki` to build
    # a `planetgen.wiki.WikiClient(backend=..., base_url=..., ...)` per
    # request, and by `GET /api/wiki-config` so the CGI browser knows
    # which backend(s) to offer without duplicating this resolution.
    WIKI_CONFIG = _wiki_config(_settings.wiki)

    # The admin session cookie (auth.py) is Secure by default -- never
    # sent over plain HTTP -- matching this project's documented
    # deployment (HTTPS in front of the app, docs/deployment/).
    # Only ever disable this for local development over plain HTTP (e.g.
    # `python src/html/wsgi.py` without TLS in front of it); a production
    # deployment must never set PLANETGEN_ADMIN_COOKIE_INSECURE or
    # config.json's admin_cookie_insecure.
    SESSION_COOKIE_SECURE = not _settings.admin_cookie_insecure

    # See `_proxy_fix` above. All 0 (the default) leaves the app as it is.
    PROXY_FIX = _settings.proxy_fix.model_dump()

    # The largest request body the app reads (SEC, TEST.47): a bigger one
    # is a 413 before any route parses it. The biggest real body is a
    # system's text sent back for download (well under 1 MB); without a
    # limit a multi-megabyte JSON or form body would be read into memory
    # whole.
    MAX_CONTENT_LENGTH = 2 * 1024 * 1024

    # See `_secret_key` above.
    SECRET_KEY = _SECRET_KEY
    SECRET_KEY_IS_EPHEMERAL = _SECRET_KEY_IS_EPHEMERAL

    # The one database the Flask-served pages (`planetgen/web/`) show. Empty
    # (the default) means `MYSQL_CONFIG.database`, i.e. `config.json`'s
    # `mysql.database` or `PLANETGEN_MYSQL_DATABASE` -- see
    # `web/helpers.db_name`. It never travels in a page URL.
    WEB_DATABASE = ""
