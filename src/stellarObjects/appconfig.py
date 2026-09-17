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
write-capable and control-schema overrides), the API's rate limits, the
admin cookie's `Secure` flag, the CGI browser's debug-page toggle, and the
site's own display name/base URL/API endpoint.

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
    "mysql": {
        "host": "127.0.0.1",
        "port": 3306,
        "user": "planetgen",
        "password": "",
        "database": "planetgen",
        "database_prefix": "planetgen",
    },
    "mysql_write": {
        # host/port/database/database_prefix always come from the `mysql`
        # section above -- the read-only and write-capable accounts share
        # one server/schema, only their credentials differ. Empty string
        # means "inherit the matching `mysql` value", same fallback
        # behavior the old `PLANETGEN_MYSQL_WRITE_*` env vars already had
        # against `PLANETGEN_MYSQL_*`.
        "user": "",
        "password": "",
    },
    "control_database": "planetgen_control",
    "ratelimit": {
        "default": "200 per day;50 per hour",
        "storage_uri": "memory://",
    },
    "admin_cookie_insecure": False,
}
"""dict: Fallback values, matching `config.json.example`'s shape, used for
any field/section missing from a deployment's real `config.json` (or when
no such file exists yet at all)."""

# TODO(wiki publishing): add a `"wiki"` section here (and to
# config.json.example) once the "Upload to Wiki" feature is wired up --
# see docs/TODO.md's "Wiki publishing isn't wired up yet" item and
# src/wikiClient/ (the already-built, standalone WikiClient library --
# one shared interface over a Wiki.js and a MediaWiki backend -- this
# would configure). Something like:
#     "wiki": {
#         "enabled": False,
#         "backend": "wikijs",  # or "mediawiki"
#         "base_url": "",
#         # wikijs backend:
#         "api_token": "",
#         # mediawiki backend:
#         "username": "",
#         "password": "",
#     },
# following this file's own documented precedence (an explicit function
# argument, then a PLANETGEN_WIKI_* env var, then this config.json
# section, then a built-in default) -- see html/api/config.py's own TODO
# for where `PLANETGEN_WIKI_BACKEND`/`PLANETGEN_WIKI_BASE_URL`/etc. would
# be read and layered on top of this section, the same way
# `_write_mysql_config` already layers PLANETGEN_MYSQL_WRITE_* over
# config.json's `mysql_write` section above.


def _merge(base, overrides):
    """Recursively merges `overrides` onto `base` in place -- a section
    (`mysql`, `mysql_write`, `ratelimit`) present in `config.json` only
    needs to name the fields it wants to change; anything it omits keeps
    its `DEFAULT_CONFIG` value rather than disappearing."""
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
