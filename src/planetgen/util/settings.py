# planetgen/util/settings.py

"""
The settings model (ADM.42)
===========================

One Pydantic model, `Settings`, describes every `config.json` option: its
type, default, help text, unit, whether it is a secret, whether changing it
needs a restart, which environment variable overrides it and whether the
web may edit it. The loader (`get_settings`), the generated tables in
`docs/config.md`, `config.json.example` and `config.schema.json`
(`planetgen.cli.config`) all read the model, and tests fail when they
drift or when code reads a `PLANETGEN_*` variable the model doesn't name.

Two files, one model (docs/design/settings-seo-and-accounts.md A2):
`config.json` (repo root; root-owned, read-only to the web process) and
`settings.json` (web-owned overlay that ADM.43's page will write; only
fields marked editable may appear in it). Load order, lowest to highest:
the model's defaults, `config.json`, `settings.json`, the `PLANETGEN_*`
environment variables, then an explicit argument at the call site.

Field metadata sits in `json_schema_extra` under `x-` keys, written with
`opt` below:

- `x-env`: the environment variable that overrides it
- `x-env-empty`: an environment variable that is set but empty still
  overrides (a blank MySQL password), instead of counting as unset
- `x-restart`: a change only applies after a restart
- `x-secret`: never shown in a page or log
- `x-editable`: `false` means file-only; `config.json` only
- `x-min_role`: `admin` (the default) or `owner`: who may save it
- `x-unit`, `x-category`: the form's suffix and section

The log options (`debug`, `log_file`, `log_dir`, `log_rotation`) are read
by `planetgen.util.logpaths`, which needs only the standard library (the
installers run it as root before any virtual environment exists); their
defaults come from there, and a test checks the two agree.

Python 3.9: annotations are `Optional[...]`/`List[...]`, never `X | None`.
"""

import hashlib
import json
import os
import sys
import threading
from typing import Any, List, Literal, get_args, get_origin

from pydantic import BaseModel, ConfigDict, Field, ValidationError, field_validator

from planetgen.util import logpaths
from planetgen.util.logpaths import PROJECT_ROOT

CONFIG_VERSION = 1
"""int: The `config_version` this code writes and reads. A file without one
is version 1; a newer one is refused."""

ENV_PREFIX = "PLANETGEN_"

_FALSE_STRINGS = ("", "0", "false", "no", "off")


def opt(default, description, env=None, restart=False, secret=False, editable=True, role="admin",
        unit=None, category="general", env_empty=False, **constraints):
    """A model field with its `x-` metadata (see the module docstring)."""
    extra = {"x-restart": restart, "x-secret": secret, "x-editable": editable, "x-min_role": role,
             "x-category": category}
    if env:
        extra["x-env"] = env
    if env_empty:
        extra["x-env-empty"] = True
    if unit:
        extra["x-unit"] = unit
    return Field(default, description=description, json_schema_extra=extra, **constraints)


def _loose_bool(value):
    """`false`, `0`, `no`, `off` and empty text mean off, any other text on
    (what `PLANETGEN_DEBUG` always did)."""
    if isinstance(value, str):
        return value.strip().lower() not in _FALSE_STRINGS
    return value


class _Section(BaseModel):
    model_config = ConfigDict(extra="ignore", validate_default=True)


class MySQL(_Section):
    host: str = opt("127.0.0.1", "The MySQL or MariaDB server's host name or address.", env="PLANETGEN_MYSQL_HOST",
                    restart=True, editable=False, category="database", env_empty=True)
    port: int = opt(3306, "The server's port.", env="PLANETGEN_MYSQL_PORT", restart=True, editable=False,
                    category="database", ge=1, le=65535, env_empty=True)
    user: str = opt("planetgen", "The account every entry point connects with; give it whatever grants the most "
                    "demanding caller needs.", env="PLANETGEN_MYSQL_USER", restart=True, editable=False,
                    category="database", env_empty=True)
    password: str = opt("", "That account's password.", env="PLANETGEN_MYSQL_PASSWORD", restart=True, secret=True,
                        editable=False, category="database", env_empty=True)
    database: str = opt("planetgen", "The galaxy database the pages show and the generators fill.",
                        env="PLANETGEN_MYSQL_DATABASE", restart=True, editable=False, category="database",
                        env_empty=True)
    database_prefix: str = opt("planetgen", "The schema-name prefix the database picker lists, for a deployment "
                               "with several galaxies on one server.", env="PLANETGEN_MYSQL_DATABASE_PREFIX",
                               restart=True, editable=False, category="database")
    statement_timeout_seconds: float = opt(
        10, "The longest any statement on the web interface's and API's read-only connections may run (0 turns "
        "the limit off). A query it stops gives a 504 page, or a QUERY_TIMEOUT error from the API.",
        unit="seconds", category="database", ge=0)


class Redis(_Section):
    url: str = opt("redis://127.0.0.1:6379/0", "The Redis server the work queue, the rate limits and the login "
                   "lockouts use.", env="PLANETGEN_REDIS_URL", restart=True, editable=False, category="queue")


class RatePages(_Section):
    search: str = opt("30 per minute", "Limit per client address on /search and /galaxy/locate. Empty turns it off.",
                      restart=True, category="rate limits")
    galaxy: str = opt("60 per minute", "Limit on /galaxy, the Galaxy Map page. Empty turns it off.", restart=True,
                      category="rate limits")
    galaxy_tiles: str = opt("600 per minute", "Limit on /galaxy/tiles, /galaxy/stage and /galaxy/territories, "
                            "fetched as the map's camera moves. Empty turns it off.", restart=True,
                            category="rate limits")
    health: str = opt("60 per minute", "Limit on /api/health. Empty turns it off.", restart=True,
                      category="rate limits")
    other: str = opt("300 per minute", "Limit on every other page, counted together. Empty turns it off.",
                     restart=True, category="rate limits")


class RateLimit(_Section):
    default: str = opt("200 per day;50 per hour", "Flask-Limiter's default limit for the API.",
                       env="PLANETGEN_RATELIMIT_DEFAULT", restart=True, category="rate limits")
    storage_uri: str = opt("", "Where limits count; empty means the Redis server in redis.url, memory:// counts "
                           "per process (right with exactly one worker).", env="PLANETGEN_RATELIMIT_STORAGE_URI",
                           restart=True, category="rate limits")
    pages: RatePages = RatePages()


class ProxyFix(_Section):
    @field_validator("x_for", "x_proto", "x_host", mode="before")
    @classmethod
    def _not_a_bool(cls, value):
        if isinstance(value, bool):
            raise ValueError("must be a whole number of proxies (0 = off), not true or false")
        return value

    x_for: int = opt(0, "Reverse proxies to trust for X-Forwarded-For, the client address (0 = off).",
                     env="PLANETGEN_PROXY_FIX_X_FOR", restart=True, editable=False, category="proxy", ge=0)
    x_proto: int = opt(0, "Reverse proxies to trust for X-Forwarded-Proto, the scheme (0 = off).",
                       env="PLANETGEN_PROXY_FIX_X_PROTO", restart=True, editable=False, category="proxy", ge=0)
    x_host: int = opt(0, "Reverse proxies to trust for X-Forwarded-Host (0 = off).",
                      env="PLANETGEN_PROXY_FIX_X_HOST", restart=True, editable=False, category="proxy", ge=0)


class TileCache(_Section):
    dir: str = opt("", "Where the Galaxy Map's tiles are cached on disk; empty means /var/cache/planetgen/tiles.",
                   env="PLANETGEN_TILE_CACHE_DIR", editable=False, category="cache")
    max_mb: float = opt(200, "Roughly how big the tile cache may grow before its oldest files are pruned (0 turns "
                        "the disk cache off).", env="PLANETGEN_TILE_CACHE_MAX_MB", unit="MB", category="cache", ge=0)


class PageCache(_Section):
    enabled: bool = opt(True, "Keep the API's public answers in memory, per WSGI process.",
                        env="PLANETGEN_PAGE_CACHE", restart=True, category="cache")
    max_entries: int = opt(2000, "Most answers kept.", restart=True, category="cache", ge=1)
    max_mb: float = opt(64, "Most memory the answers may take.", restart=True, unit="MB", category="cache", gt=0)
    stamp_seconds: float = opt(15, "How often the galaxy's content stamp is checked.", restart=True, unit="seconds",
                               category="cache", ge=0)
    max_age_seconds: float = opt(300, "Nothing is kept longer than this.", restart=True, unit="seconds",
                                 category="cache", gt=0)

    _loose = field_validator("enabled", mode="before")(_loose_bool)


class Jobs(_Section):
    dir: str = opt("", "Where each background job's command lines, status and output are kept; empty means "
                   "/var/lib/planetGen/jobs.", env="PLANETGEN_JOBS_DIR", editable=False, category="jobs")
    keep: int = opt(20, "How many finished jobs are kept.", category="jobs", ge=1)
    python: str = opt("", "The interpreter that runs planetgen and planetgen.cli.reset for the jobs; empty means "
                      "the web app's own Python. Editing it from the web would be code execution, so it is file-only.",
                      env="PLANETGEN_PYTHON", editable=False, category="jobs")


class WikiJs(_Section):
    base_url: str = opt("", "The Wiki.js instance's root URL. Empty means it isn't offered as an upload target.",
                        env="PLANETGEN_WIKIJS_BASE_URL", role="owner", category="wiki")
    api_token: str = opt("", "A Wiki.js Personal API Token (Admin, API Access).", env="PLANETGEN_WIKIJS_API_TOKEN",
                         secret=True, role="owner", category="wiki")


class MediaWiki(_Section):
    base_url: str = opt("", "The MediaWiki API entry point's directory (everything up to, not including, api.php). "
                        "Empty means it isn't offered as an upload target.", env="PLANETGEN_MEDIAWIKI_BASE_URL",
                        role="owner", category="wiki")
    username: str = opt("", "A Bot Password user name, in User@BotName form.", env="PLANETGEN_MEDIAWIKI_USERNAME",
                        role="owner", category="wiki")
    password: str = opt("", "That Bot Password.", env="PLANETGEN_MEDIAWIKI_PASSWORD", secret=True, role="owner",
                        category="wiki")


class Wiki(_Section):
    wikijs: WikiJs = WikiJs()
    mediawiki: MediaWiki = MediaWiki()


class Settings(_Section):
    """Every `config.json` option."""

    config_version: int = opt(CONFIG_VERSION, "The settings file format version; a file without one is version 1.",
                              editable=False, category="file", ge=1)
    site_name: str = opt("planetGen", "Display name for this deployment, shown in the header, the home page "
                         "heading, every page title and the footer.", category="site", min_length=1)
    base_url: str = opt("http://localhost/", "The base URL this deployment is served from. Not yet read by any "
                        "page; reserved for absolute links that a request alone can't give.", role="owner",
                        category="site")
    api_base_url: str = opt("http://127.0.0.1/api", "Base URL of the Flask API's /api mount point, used by "
                            "web/lib/apiclient.py only outside the Flask app.", env="PLANETGEN_API_BASE_URL",
                            restart=True, editable=False, category="site")
    debug: bool = opt(False, "Write a verbose debug log to log_file: every decision the generator makes, every "
                      "random roll, every SQL statement and request. The log grows fast.", env="PLANETGEN_DEBUG",
                      restart=True, editable=False, category="logging")
    log_file: str = opt(logpaths.DEFAULT_LOG_FILE, "Where the debug log goes (set it explicitly on Windows).",
                        env="PLANETGEN_LOG_FILE", restart=True, editable=False, category="logging")
    log_dir: str = opt("", "The folder of the always-on activity log; empty means /var/log/planetgen on Linux, "
                       "/Library/Logs/planetgen on macOS, logs under the checkout on Windows.",
                       env="PLANETGEN_LOG_DIR", restart=True, editable=False, category="logging")
    log_rotation: Literal["auto", "system", "app"] = opt(
        logpaths.DEFAULT_LOG_ROTATION, "How the activity log is rotated: system (logrotate or newsyslog), app (the "
        "program, 100 MB per file, 30 old copies) or auto (system where an installed rotation config exists, else "
        "app).", restart=True, editable=False, category="logging")
    mysql: MySQL = MySQL()
    control_database: str = opt("planetgen_control", "The MySQL schema holding admin identities, sessions, API keys "
                                "and the audit log. Never listed or selectable as a galaxy.",
                                env="PLANETGEN_CONTROL_DATABASE", restart=True, editable=False, category="database")
    redis: Redis = Redis()
    ratelimit: RateLimit = RateLimit()
    login_allowlist: List[str] = opt(
        [], "Client addresses or networks (CIDR) never locked out after failed logins; loopback never is either. "
        "Environment: comma- or space-separated.", env="PLANETGEN_LOGIN_ALLOWLIST", restart=True, role="owner",
        category="security")
    admin_cookie_insecure: bool = opt(False, "Send the admin session cookie over plain HTTP. Only for local "
                                      "development without TLS; a production deployment must never set it.",
                                      env="PLANETGEN_ADMIN_COOKIE_INSECURE", restart=True, editable=False,
                                      category="security")
    secret_key: str = opt("", "Signs the CSRF tokens on the pages' forms. Empty makes a random one at startup (a "
                          "form open across a restart then fails once).", env="PLANETGEN_SECRET_KEY", restart=True,
                          secret=True, editable=False, category="security")
    proxy_fix: ProxyFix = ProxyFix()
    tile_cache: TileCache = TileCache()
    page_cache: PageCache = PageCache()
    jobs: Jobs = Jobs()
    wiki: Wiki = Wiki()

    _loose = field_validator("debug", "admin_cookie_insecure", mode="before")(_loose_bool)

    @field_validator("login_allowlist", mode="before")
    @classmethod
    def _split_allowlist(cls, value):
        if isinstance(value, str):
            return value.replace(",", " ").split()
        return value


class SettingsError(ValueError):
    """A settings file or environment variable that doesn't fit the model;
    the message names each field."""


# ---------------------------------------------------------------------
# Walking the model
# ---------------------------------------------------------------------

def _model_of(annotation):
    if isinstance(annotation, type) and issubclass(annotation, BaseModel):
        return annotation
    return None


def iter_fields(model=Settings, prefix=""):
    """
    Yields `(dotted path, FieldInfo, metadata dict)` for every leaf option
    of `model` (sections are walked into). The metadata always has every
    `x-` key (`opt`'s defaults), `description` and the annotation.
    """
    for name, info in model.model_fields.items():
        path = f"{prefix}{name}"
        section = _model_of(info.annotation)
        if section is not None:
            yield from iter_fields(section, f"{path}.")
            continue
        meta = dict(info.json_schema_extra or {})
        meta["description"] = info.description or ""
        meta["annotation"] = info.annotation
        yield path, info, meta


def field_default(info):
    return info.get_default(call_default_factory=True)


def type_name(annotation):
    """`"integer"`, `"text"`, `"list of text"`, ... for the docs table."""
    origin = get_origin(annotation)
    if origin is Literal:
        return " or ".join(f'"{a}"' for a in get_args(annotation))
    if origin in (list, List):
        return "list of text"
    return {bool: "true or false", int: "whole number", float: "number", str: "text"}.get(annotation, str(annotation))


def env_names():
    """Every environment variable the model names, as `{variable: path}`."""
    return {meta["x-env"]: path for path, _info, meta in iter_fields() if meta.get("x-env")}


def env_name(path):
    """The environment variable the model names for option `path`."""
    for found, _info, meta in iter_fields():
        if found == path:
            return meta["x-env"]
    raise KeyError(path)


# ---------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------

def default_web_settings_path():
    """`settings.json`, the web-owned overlay: `PLANETGEN_SETTINGS_FILE`,
    else `/var/lib/planetGen/settings.json` (the checkout's `settings.json`
    on Windows)."""
    env = os.environ.get("PLANETGEN_SETTINGS_FILE")
    if env:
        return env
    if sys.platform.startswith("win"):
        return os.path.join(PROJECT_ROOT, "settings.json")
    return "/var/lib/planetGen/settings.json"


def _read_json_object(path, label):
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
    except FileNotFoundError:
        return {}
    except (OSError, ValueError) as exc:
        raise SettingsError(f"{label} ({path}) can't be read: {exc}") from exc
    if not isinstance(data, dict):
        raise SettingsError(f"{label} ({path}) must contain a JSON object, not {type(data).__name__}")
    return data


def migrate(data):
    """
    Brings a file's dict up to `CONFIG_VERSION` (no migrations exist yet:
    version 1 is the first).

    Raises:
        SettingsError: The file is from a newer release.
    """
    version = data.get("config_version", 1)
    if not isinstance(version, int) or isinstance(version, bool) or version < 1:
        raise SettingsError(f"config_version must be a whole number 1 or above, not {version!r}")
    if version > CONFIG_VERSION:
        raise SettingsError(f"config_version {version} is newer than this release understands ({CONFIG_VERSION}); "
                            "update planetGen")
    return data


def _merge(base, overrides):
    for key, value in overrides.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            _merge(base[key], value)
        else:
            base[key] = value


def _set_path(data, path, value):
    parts = path.split(".")
    for part in parts[:-1]:
        if not isinstance(data.get(part), dict):
            data[part] = {}
        data = data[part]
    data[parts[-1]] = value


def _environment_layer(environ):
    layer = {}
    for path, _info, meta in iter_fields():
        name = meta.get("x-env")
        if not name or name not in environ:
            continue
        value = environ[name]
        if value == "" and not meta.get("x-env-empty"):
            continue
        _set_path(layer, path, value)
    return layer


def _file_only_paths(overlay, prefix=""):
    """The dotted paths in `overlay` the web may not set (`x-editable`
    false), with unknown ones: both are errors in `settings.json`."""
    editable = {path: meta.get("x-editable", True) for path, _i, meta in iter_fields()}
    bad = []

    def walk(node, prefix):
        for key, value in node.items():
            path = f"{prefix}{key}"
            if isinstance(value, dict) and not any(p == path for p in editable):
                walk(value, f"{path}.")
            elif path not in editable:
                bad.append(f"{path} is not a setting")
            elif not editable[path]:
                bad.append(f"{path} can only be set in config.json")
    walk(overlay, prefix)
    return bad


def build(config=None, overlay=None, environ=None, strict=False):
    """
    Validates the layers into a `Settings`: defaults, `config` (the
    `config.json` dict), `overlay` (the `settings.json` dict), then the
    `PLANETGEN_*` variables in `environ`.

    Args:
        strict (bool): Unknown keys in `config` are errors (`config check`,
            ADM.43's form) instead of ignored with a warning.

    Raises:
        SettingsError: Naming every field that doesn't fit.
    """
    environ = os.environ if environ is None else environ
    merged = {}
    problems = []
    for label, layer in (("config.json", config), ("settings.json", overlay)):
        if layer is None:
            continue
        try:
            layer = migrate(dict(layer))
        except SettingsError as exc:
            raise SettingsError(f"{label}: {exc}") from exc
        if label == "settings.json":
            problems.extend(f"settings.json: {message}" for message in _file_only_paths(layer))
        elif strict:
            known = {path for path, _i, _m in iter_fields()}
            problems.extend(f"config.json: {message}" for message in _unknown_paths(layer, known))
        _merge(merged, layer)
    _merge(merged, _environment_layer(environ))
    if problems:
        raise SettingsError("\n".join(problems))
    try:
        return Settings.model_validate(merged)
    except ValidationError as exc:
        lines = []
        for error in exc.errors():
            where = ".".join(str(part) for part in error["loc"])
            lines.append(f"{where}: {error['msg']}")
        raise SettingsError("; ".join(lines) if len(lines) < 4 else "\n".join(lines)) from exc


def _unknown_paths(layer, known, prefix=""):
    sections = {path.rsplit(".", 1)[0] for path in known if "." in path}
    out = []
    for key, value in layer.items():
        path = f"{prefix}{key}"
        if path in known:
            continue
        if isinstance(value, dict) and (path in sections or any(k.startswith(f"{path}.") for k in known)):
            out.extend(_unknown_paths(value, known, f"{path}."))
        else:
            out.append(f"{path} is not a setting")
    return out


_lock = threading.Lock()
_cache = {"key": None, "settings": None}


def _stamp(path):
    try:
        stat = os.stat(path)
    except OSError:
        return None
    return (stat.st_mtime_ns, stat.st_size)


def _environment_stamp(environ):
    return tuple(sorted((name, environ[name]) for name in env_names() if name in environ))


def get_settings(environ=None):
    """
    The validated settings. Re-read only when `config.json`,
    `settings.json` or a variable the model names changes (a stat costs
    about a microsecond), so every option that is read per use stays live.

    Raises:
        SettingsError: A file or variable doesn't fit the model.
    """
    environ = os.environ if environ is None else environ
    config_path = logpaths.CONFIG_PATH
    overlay_path = default_web_settings_path()
    key = (config_path, _stamp(config_path), overlay_path, _stamp(overlay_path), _environment_stamp(environ))
    with _lock:
        if _cache["key"] == key:
            return _cache["settings"]
        config = _read_json_object(config_path, "config.json")
        overlay = _read_json_object(overlay_path, "settings.json") if os.path.isfile(overlay_path) else None
        settings = build(config, overlay, environ)
        _cache["key"], _cache["settings"] = key, settings
        return settings


def reset_cache():
    """Forgets the cached settings (tests)."""
    with _lock:
        _cache["key"] = _cache["settings"] = None


# ---------------------------------------------------------------------
# What the model writes: the docs table, the example file, the schema
# ---------------------------------------------------------------------

def example_dict():
    """`config.json.example`'s content: every default of the options
    `config.json` may hold, in the model's order."""
    return Settings().model_dump(mode="json")


def example_json():
    return json.dumps(example_dict(), indent=2) + "\n"


def schema_json():
    """`config.schema.json`: the model's JSON Schema."""
    return json.dumps(Settings.model_json_schema(), indent=2, sort_keys=True) + "\n"


def _format_default(value):
    if isinstance(value, bool):
        return "true" if value else "false"
    if value == "":
        return "empty"
    if isinstance(value, list):
        return "empty list" if not value else ", ".join(str(v) for v in value)
    return str(value)


def docs_table():
    """The generated option table for `docs/config.md` (Markdown)."""
    rows = ["| Option | Type | Default | Applies | Edited | Environment | What it does |",
            "|---|---|---|---|---|---|---|"]
    roles = {"admin": "any admin", "owner": "the Owner"}
    for path, info, meta in iter_fields():
        applies = "after a restart" if meta["x-restart"] else "at once"
        edited = roles[meta["x-min_role"]] if meta["x-editable"] else "config.json only"
        unit = f" ({meta['x-unit']})" if meta.get("x-unit") else ""
        secret = " **Secret.**" if meta["x-secret"] else ""
        env = f"`{meta['x-env']}`" if meta.get("x-env") else "none"
        rows.append(f"| `{path}` | {type_name(meta['annotation'])}{unit} | {_format_default(field_default(info))} | "
                    f"{applies} | {edited} | {env} | {meta['description']}{secret} |")
    return "\n".join(rows) + "\n"


def settings_hash(settings=None):
    """A short hash of the effective settings, secrets left out (for a
    page that detects a concurrent edit)."""
    settings = settings or get_settings()
    secret_paths = {path for path, _i, meta in iter_fields() if meta["x-secret"]}
    flat = {path: _lookup(settings, path) for path, _i, _m in iter_fields() if path not in secret_paths}
    return hashlib.sha256(json.dumps(flat, sort_keys=True, default=str).encode()).hexdigest()[:16]


def _lookup(settings, path):
    node: Any = settings
    for part in path.split("."):
        node = getattr(node, part)
    return node
