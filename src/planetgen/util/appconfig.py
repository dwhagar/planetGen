# planetgen/util/appconfig.py

"""
Unified deployment configuration loader
========================================

Loads `config.json` -- the single, per-deployment settings file for the
whole planetGen project (generation CLIs, the Flask API, and the `html/`
CGI browser alike), documented in full in
[`config.md`](../../docs/config.md). That file covers what each field
means, why the real `config.json` is gitignored while `config.json.example`
is committed as a template, and how this relates to the `PLANETGEN_*`
environment variables every entry point already reads.

This replaces the old `webconfig.py`/`webconfig.json`, which only ever
covered a handful of unused placeholder fields for the web interface.
`config.json` instead holds every setting a deployment would otherwise
have to repeat across Apache `SetEnv` lines, systemd `EnvironmentFile`s,
and `--mysql-*` CLI flags -- MySQL connection details (including the
control-schema name), the API's rate limits, the admin cookie's `Secure`
flag, the debug log switch, and the site's own display
name/base URL/API endpoint.

Precedence, everywhere a setting has more than one source, is: an
explicit function/CLI argument, then the matching `PLANETGEN_*`
environment variable (so a single shared `config.json` can still be
overridden per-process -- e.g. `planetgen-orbits@.service`'s per-instance
`PLANETGEN_MYSQL_DATABASE=%i`), then `config.json`, then the built-in
default below. Callers (`planetgen.db.store`, `html/api/config.py`,
`html/lib/apiclient.py`) each still read their own
`os.environ.get(VAR, ...)` for the middle two steps; this module only
supplies the `config.json` layer.

Dependency-free (standard library `json`/`os`/`copy` only), matching the
rest of this project's config plumbing.
"""

import copy
import json
import os
import sys

# This file is src/planetgen/util/appconfig.py: four levels up is the repo
# root.
_PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

CONFIG_PATH = os.path.join(_PROJECT_ROOT, "config.json")
"""str: Absolute path to the real config file, at the repo root -- a
sibling of `html/`, `db/`, and `src/`, deliberately outside `html/`'s
Apache `DocumentRoot` (unlike `config.json.example`, which lives at the
repo root too since it holds only placeholder values, not secrets). See
`docs/config.md` for why the location matters."""

DEFAULT_CONFIG = {
    "site_name": "planetGen",
    "base_url": "http://localhost/",
    "api_base_url": "http://127.0.0.1/api",
    "debug": False,
    "log_file": "/var/log/planetgen.log",
    # The always-on activity log's folder (planetgen/admin/activity_log.py).
    # Empty means the platform's standard place (`default_log_dir`).
    "log_dir": "",
    # "auto", "system" (logrotate/newsyslog moves the file) or "app" (the
    # program rotates it itself); see `log_rotation_mode`.
    "log_rotation": "auto",
    "mysql": {
        "host": "127.0.0.1",
        "port": 3306,
        "user": "planetgen",
        "password": "",
        "database": "planetgen",
        "database_prefix": "planetgen",
        # PERF.17: the web's read-only connections stop any statement
        # that runs longer than this many seconds (0 turns it off).
        "statement_timeout_seconds": 10,
    },
    "control_database": "planetgen_control",
    # The Redis server the work queue (PERF.24) and the rate limits
    # (SEC.30) will use; install.sh and update.sh install or check it
    # (OPS.21). Nothing reads it yet.
    "redis": {
        "url": "redis://127.0.0.1:6379/0",
    },
    "ratelimit": {
        "default": "200 per day;50 per hour",
        "storage_uri": "memory://",
        # Per-client-IP limits on the HTML pages (html/web/ratelimits.py)
        # and /api/health. An empty string turns that one limit off.
        "pages": {
            "search": "30 per minute",
            "galaxy": "60 per minute",
            "galaxy_tiles": "600 per minute",
            "health": "60 per minute",
            "other": "300 per minute",
        },
    },
    # Client addresses or networks (CIDR) never locked out after failed
    # logins (html/api/loginguard.py); loopback never is either.
    "login_allowlist": [],
    "admin_cookie_insecure": False,
    "secret_key": "",
    # How many reverse proxies in front of the WSGI server to trust for
    # each X-Forwarded-* header (werkzeug's ProxyFix; see
    # html/api/config.py). All 0 (the default) means the app uses the
    # connection's own address and scheme, which is right under Apache +
    # mod_wsgi.
    "proxy_fix": {
        "x_for": 0,
        "x_proto": 0,
        "x_host": 0,
    },
    "tile_cache": {
        "dir": "",
        "max_mb": 200,
    },
    "jobs": {
        "dir": "",
        "keep": 20,
        "python": "",
    },
    "wiki": {
        # Either, both, or neither backend may be configured at once -- a
        # backend is "configured" (offered as an upload target) purely by
        # having a non-empty base_url plus that backend's own required
        # credential field(s), not by a separate on/off switch (see
        # html/api/config.py's WIKI_CONFIG, which computes this). Both
        # configured at once is exactly what lets a caller choose "the
        # wiki of their choice" per upload (see wikiClient/client.py).
        "wikijs": {
            "base_url": "",
            "api_token": "",
        },
        "mediawiki": {
            "base_url": "",
            "username": "",
            "password": "",
        },
    },
}
"""dict: Fallback values, matching `config.json.example`'s shape, used for
any field/section missing from a deployment's real `config.json` (or when
no such file exists yet at all)."""

def _merge(base, overrides):
    """Recursively merges `overrides` onto `base` in place -- a section
    (`mysql`, `ratelimit`) present in `config.json` only needs to name the
    fields it wants to change; anything it omits keeps its
    `DEFAULT_CONFIG` value rather than disappearing."""
    for key, value in overrides.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            _merge(base[key], value)
        else:
            base[key] = value


def load_config():
    """
    Loads the deployment configuration, deep-merged onto `DEFAULT_CONFIG`
    so every field/section is always present -- callers never need to
    handle a missing file, or a file that only sets a handful of fields,
    themselves.

    Returns:
        dict: `DEFAULT_CONFIG`'s shape, with any values `config.json`
              overrides applied on top.
    """
    merged = copy.deepcopy(DEFAULT_CONFIG)
    if not os.path.isfile(CONFIG_PATH):
        return merged

    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        overrides = json.load(f)
    if not isinstance(overrides, dict):
        # Same failure mode as malformed JSON (json.JSONDecodeError is a
        # ValueError), with a message that says what's actually wrong.
        raise ValueError(f"{CONFIG_PATH} must contain a JSON object, not {type(overrides).__name__}")
    _merge(merged, overrides)
    return merged


_FALSE_STRINGS = ("", "0", "false", "no", "off")


def debug_enabled(config=None):
    """
    Whether this deployment's debug mode is on: the `PLANETGEN_DEBUG`
    environment variable when set (`0`/`false`/`no`/`off`/empty mean off,
    anything else on), else `config.json`'s `"debug"`, which defaults to
    off when missing (a string there is read the same way as the variable). Debug mode turns on the verbose debug log (see
    `planetgen.util.log`) and the web interface's traceback-in-page 500
    responses.

    Args:
        config (dict, optional): An already-loaded `load_config()` result,
            to skip re-reading the file.

    Returns:
        bool
    """
    env = os.environ.get("PLANETGEN_DEBUG")
    if env is not None:
        return env.strip().lower() not in _FALSE_STRINGS
    if config is None:
        config = load_config()
    value = config.get("debug")
    if isinstance(value, str):
        # "false"/"0"/"off"/"no" in config.json mean off, exactly as they
        # do for PLANETGEN_DEBUG -- bool("false") would be True.
        return value.strip().lower() not in _FALSE_STRINGS
    return bool(value)


def log_file_path(config=None):
    """
    Where the debug log goes: `PLANETGEN_LOG_FILE` when set, else
    `config.json`'s `"log_file"`, else `/var/log/planetgen.log`.

    Args:
        config (dict, optional): An already-loaded `load_config()` result.

    Returns:
        str
    """
    env = os.environ.get("PLANETGEN_LOG_FILE")
    if env:
        return env
    if config is None:
        config = load_config()
    return config.get("log_file") or DEFAULT_CONFIG["log_file"]


ACTIVITY_LOG_NAME = "planetgen.log"
"""str: The activity log's file name inside `log_dir_path()`."""

SYSTEM_ROTATION_FILES = ("/etc/logrotate.d/planetgen-log", "/etc/newsyslog.d/planetgen-log.conf")
"""tuple: The rotation configs `examples/apache/setup-debug-log.sh`
installs for the activity log (Linux, macOS). While one exists, `"auto"`
rotation leaves the file to the system tool."""


def default_log_dir(platform=None):
    """
    The platform's standard folder for the activity log: `/var/log/
    planetgen` on Linux (and other Unix systems), `/Library/Logs/planetgen`
    on macOS, and `logs` under the checkout on Windows.

    Args:
        platform (str, optional): A `sys.platform` value (default this
            one's).
    """
    platform = platform or sys.platform
    if platform.startswith("win"):
        return os.path.join(_PROJECT_ROOT, "logs")
    if platform == "darwin":
        return "/Library/Logs/planetgen"
    return "/var/log/planetgen"


def log_dir_path(config=None):
    """
    The activity log's folder: `PLANETGEN_LOG_DIR` when set, else
    `config.json`'s `"log_dir"`, else `default_log_dir()`.
    """
    env = os.environ.get("PLANETGEN_LOG_DIR")
    if env:
        return env
    if config is None:
        config = load_config()
    value = config.get("log_dir")
    if value and not isinstance(value, str):
        raise TypeError(f'"log_dir" must be a string, not {type(value).__name__}')
    return value or default_log_dir()


def activity_log_path(config=None):
    """
    The activity log file: `planetgen.log` in `log_dir_path()`, or
    `planetgen-activity.log` there when that would be the very file the
    debug log uses (`log_file_path()`), so the two never share a file.
    """
    if config is None:
        config = load_config()
    folder = os.path.abspath(log_dir_path(config))
    path = os.path.join(folder, ACTIVITY_LOG_NAME)
    if os.path.normcase(path) == os.path.normcase(os.path.abspath(log_file_path(config))):
        path = os.path.join(folder, "planetgen-activity.log")
    return path


def log_rotation_mode(config=None):
    """
    How the activity log is rotated: `"system"` (logrotate or newsyslog
    moves the file; the program reopens it) or `"app"` (the program
    rotates it itself). `config.json`'s `"log_rotation"` picks one;
    `"auto"` (the default) means `"system"` where one of
    `SYSTEM_ROTATION_FILES` exists, else `"app"` (always on Windows).
    """
    if config is None:
        config = load_config()
    value = config.get("log_rotation")
    value = value.strip().lower() if isinstance(value, str) else "auto"
    if value in ("system", "app"):
        return value
    if sys.platform.startswith("win"):
        return "app"
    return "system" if any(os.path.exists(p) for p in SYSTEM_ROTATION_FILES) else "app"
