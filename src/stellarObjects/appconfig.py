# stellarObjects/appconfig.py

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
default below. Callers (`stellarObjects._db`, `html/api/config.py`,
`html/lib/apiclient.py`, `html/lib/page.py`) each still read their own
`os.environ.get(VAR, ...)` for the middle two steps; this module only
supplies the `config.json` layer.

Dependency-free (standard library `json`/`os`/`copy` only), matching the
rest of this project's config plumbing.
"""

import copy
import json
import os

# stellarObjects/ lives at src/stellarObjects/ (src layout) -- three levels
# up from this file, not two, to reach the actual repo root.
_PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

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
    "mysql": {
        "host": "127.0.0.1",
        "port": 3306,
        "user": "planetgen",
        "password": "",
        "database": "planetgen",
        "database_prefix": "planetgen",
    },
    "control_database": "planetgen_control",
    "ratelimit": {
        "default": "200 per day;50 per hour",
        "storage_uri": "memory://",
    },
    "admin_cookie_insecure": False,
    "tile_cache": {
        "dir": "",
        "max_mb": 200,
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
    _merge(merged, overrides)
    return merged


_FALSE_STRINGS = ("", "0", "false", "no", "off")


def debug_enabled(config=None):
    """
    Whether this deployment's debug mode is on: the `PLANETGEN_DEBUG`
    environment variable when set (`0`/`false`/`no`/`off`/empty mean off,
    anything else on), else `config.json`'s `"debug"`, which defaults to
    off when missing. Debug mode turns on the verbose debug log (see
    `stellarObjects.log`) and the web interface's traceback-in-page 500
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
    return bool(config.get("debug"))


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
