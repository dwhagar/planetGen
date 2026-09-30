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

Every read and write in `routes.py` goes through this same connection --
`WRITE_MYSQL_CONFIG`/`CONTROL_MYSQL_CONFIG` both simply reuse
`MYSQL_CONFIG`, there's no separate write-capable account to configure.
Give this one account whatever grants the most demanding caller needs
(`INSERT`/`UPDATE`/`DELETE` on the content schemas and the control schema,
at minimum, once any write endpoint is reachable) -- see
`docs/apache-deployment.md`'s "MySQL accounts" section.

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
import secrets

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


def _wiki_config(wiki_defaults):
    """
    Builds the resolved wiki-publishing config from `PLANETGEN_WIKIJS_*`/
    `PLANETGEN_MEDIAWIKI_*` environment variables, falling back to
    `config.json`'s `wiki` section (`wiki_defaults`) -- same precedence
    (explicit env var, then `config.json`, then `""` = "not set") every
    other layered setting in this file already follows.

    Neither backend has a separate on/off flag: a backend counts as
    configured -- and so is offered as an "Upload to Wiki" target by
    `routes.py`'s `POST /api/systems/<id>/wiki`/`POST /api/sectors/<id>/wiki`
    -- purely by having a non-empty `base_url` plus every credential field
    its own `wikiClient` backend requires (see `_configured` below); an
    empty `base_url` alone already means "don't offer this one", so a
    separate flag would only ever duplicate that check.

    Returns:
        dict: `{"wikijs": {"base_url", "api_token", "configured"},
            "mediawiki": {"base_url", "username", "password", "configured"}}`.
    """
    wikijs = {
        "base_url": os.environ.get("PLANETGEN_WIKIJS_BASE_URL") or wiki_defaults["wikijs"]["base_url"],
        "api_token": os.environ.get("PLANETGEN_WIKIJS_API_TOKEN") or wiki_defaults["wikijs"]["api_token"],
    }
    wikijs["configured"] = bool(wikijs["base_url"] and wikijs["api_token"])

    mediawiki = {
        "base_url": os.environ.get("PLANETGEN_MEDIAWIKI_BASE_URL") or wiki_defaults["mediawiki"]["base_url"],
        "username": os.environ.get("PLANETGEN_MEDIAWIKI_USERNAME") or wiki_defaults["mediawiki"]["username"],
        "password": os.environ.get("PLANETGEN_MEDIAWIKI_PASSWORD") or wiki_defaults["mediawiki"]["password"],
    }
    mediawiki["configured"] = bool(mediawiki["base_url"] and mediawiki["username"] and mediawiki["password"])

    return {"wikijs": wikijs, "mediawiki": mediawiki}


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


def _secret_key():
    """
    `PLANETGEN_SECRET_KEY`, else `config.json`'s `secret_key`, else a
    random per-process key (`SECRET_KEY_IS_EPHEMERAL` is then True and
    `create_app` logs a warning). Used to sign the Flask-served pages'
    CSRF tokens (`html/web/csrf.py`); nothing else in the app signs
    anything with it today.
    """
    configured = os.environ.get("PLANETGEN_SECRET_KEY") or _config_file.get("secret_key") or ""
    if configured:
        return configured, False
    return secrets.token_hex(32), True


_SECRET_KEY, _SECRET_KEY_IS_EPHEMERAL = _secret_key()


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

    RATELIMIT_DEFAULT = os.environ.get("PLANETGEN_RATELIMIT_DEFAULT", _config_file["ratelimit"]["default"])
    RATELIMIT_STORAGE_URI = os.environ.get("PLANETGEN_RATELIMIT_STORAGE_URI", _config_file["ratelimit"]["storage_uri"])
    RATELIMIT_HEADERS_ENABLED = True

    # See `_wiki_config` above -- read by `routes.py`'s
    # `POST /api/systems/<id>/wiki`/`POST /api/sectors/<id>/wiki` to build
    # a `wikiClient.WikiClient(backend=..., base_url=..., ...)` per
    # request, and by `GET /api/wiki-config` so the CGI browser knows
    # which backend(s) to offer without duplicating this resolution.
    WIKI_CONFIG = _wiki_config(_config_file["wiki"])

    # The admin session cookie (auth.py) is Secure by default -- never
    # sent over plain HTTP -- matching this project's documented
    # deployment (Apache + Let's Encrypt, docs/apache-deployment.md).
    # Only ever disable this for local development over plain HTTP (e.g.
    # `python src/html/wsgi.py` without TLS in front of it); a production
    # deployment must never set PLANETGEN_ADMIN_COOKIE_INSECURE or
    # config.json's admin_cookie_insecure.
    SESSION_COOKIE_SECURE = _session_cookie_secure(_config_file["admin_cookie_insecure"])

    # See `_secret_key` above.
    SECRET_KEY = _SECRET_KEY
    SECRET_KEY_IS_EPHEMERAL = _SECRET_KEY_IS_EPHEMERAL

    # The one database the Flask-served pages (`html/web/`) show. Empty
    # (the default) means `MYSQL_CONFIG.database`, i.e. `config.json`'s
    # `mysql.database` or `PLANETGEN_MYSQL_DATABASE` -- see
    # `web/helpers.db_name`. It never travels in a page URL.
    WEB_DATABASE = ""
