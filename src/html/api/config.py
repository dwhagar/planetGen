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


PROXY_FIX_HEADERS = ("x_for", "x_proto", "x_host")
"""tuple: The `proxy_fix` fields, each werkzeug `ProxyFix`'s argument of
the same name: `x_for` (`X-Forwarded-For`, the client address),
`x_proto` (`X-Forwarded-Proto`, http or https) and `x_host`
(`X-Forwarded-Host`)."""


def _proxy_fix(section):
    """
    How many proxies to trust for each `X-Forwarded-*` header:
    `PLANETGEN_PROXY_FIX_X_FOR`/`_X_PROTO`/`_X_HOST` when set, else
    `config.json`'s `proxy_fix` section, else 0 (don't trust the header).

    Behind nginx, Caddy, IIS or Apache's `mod_proxy`, the WSGI server only
    sees the proxy's connection: without this every visitor has the
    proxy's address (so they all share one rate-limit bucket) and every
    request looks like plain HTTP (so no `Strict-Transport-Security`).
    `create_app` wraps the app in `ProxyFix` when any count is above 0.
    Under Apache + mod_wsgi the connection already is the client's, so
    leave them all at 0 there: a trusted header a client can set itself
    lets it pick its own rate-limit address.

    Raises:
        ValueError: A count that isn't a whole number 0 or above -- a
            typo here should stop the app, not quietly turn trust off
            (or on).
    """
    if not isinstance(section, dict):
        section = {}
    counts = {}
    for name in PROXY_FIX_HEADERS:
        env_name = f"PLANETGEN_PROXY_FIX_{name.upper()}"
        raw = os.environ.get(env_name)
        source = env_name
        if raw is None or raw.strip() == "":
            raw = section.get(name, 0)
            source = f"config.json proxy_fix.{name}"
        if isinstance(raw, bool):
            raise ValueError(f"{source} must be a whole number of proxies (0 = off), not {raw!r}")
        try:
            value = int(str(raw).strip())
        except ValueError:
            raise ValueError(f"{source} must be a whole number of proxies (0 = off), not {raw!r}") from None
        if value < 0:
            raise ValueError(f"{source} must be 0 or more, not {value}")
        counts[name] = value
    return counts


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
    # Per-IP limits on the HTML pages and /api/health (`web/ratelimits.py`):
    # `search`, `galaxy`, `galaxy_tiles`, `health`, and `other` for every
    # other page. An empty value turns that limit off.
    RATELIMIT_PAGES = dict(_config_file["ratelimit"].get("pages") or {})
    # Per-address lockout and per-username backoff (`loginguard.py`).
    # Always on in a real deployment; only a test harness that fails many
    # logins on purpose turns it off.
    LOGIN_BACKOFF_ENABLED = True
    # Addresses and networks never locked out (loopback never is either):
    # `config.json`'s `login_allowlist`, or PLANETGEN_LOGIN_ALLOWLIST
    # (comma- or space-separated).
    LOGIN_ALLOWLIST = (os.environ.get("PLANETGEN_LOGIN_ALLOWLIST")
                       or _config_file.get("login_allowlist") or [])

    # See `_wiki_config` above -- read by `routes.py`'s
    # `POST /api/systems/<id>/wiki`/`POST /api/sectors/<id>/wiki` to build
    # a `wikiClient.WikiClient(backend=..., base_url=..., ...)` per
    # request, and by `GET /api/wiki-config` so the CGI browser knows
    # which backend(s) to offer without duplicating this resolution.
    WIKI_CONFIG = _wiki_config(_config_file["wiki"])

    # The admin session cookie (auth.py) is Secure by default -- never
    # sent over plain HTTP -- matching this project's documented
    # deployment (HTTPS in front of the app, docs/deployment/).
    # Only ever disable this for local development over plain HTTP (e.g.
    # `python src/html/wsgi.py` without TLS in front of it); a production
    # deployment must never set PLANETGEN_ADMIN_COOKIE_INSECURE or
    # config.json's admin_cookie_insecure.
    SESSION_COOKIE_SECURE = _session_cookie_secure(_config_file["admin_cookie_insecure"])

    # See `_proxy_fix` above. All 0 (the default) leaves the app as it is.
    PROXY_FIX = _proxy_fix(_config_file.get("proxy_fix"))

    # See `_secret_key` above.
    SECRET_KEY = _SECRET_KEY
    SECRET_KEY_IS_EPHEMERAL = _SECRET_KEY_IS_EPHEMERAL

    # The one database the Flask-served pages (`html/web/`) show. Empty
    # (the default) means `MYSQL_CONFIG.database`, i.e. `config.json`'s
    # `mysql.database` or `PLANETGEN_MYSQL_DATABASE` -- see
    # `web/helpers.db_name`. It never travels in a page URL.
    WEB_DATABASE = ""
