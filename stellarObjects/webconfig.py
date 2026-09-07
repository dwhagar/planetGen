# stellarObjects/webconfig.py

"""
Web interface site configuration loader
========================================

Loads `webconfig.json` -- the site-level configuration file for the
planetGen web interface (`html/`), documented in full in
[`WEBCONFIG.md`](../WEBCONFIG.md) at the repo root. That file covers what
each field means, why the real `webconfig.json` is gitignored while
`webconfig.json.example` is committed as a template, and how this relates
to the `PLANETGEN_DB_DIR`/`PLANETGEN_DEBUG` environment variables `html/`
already reads.

This module is deliberately small and dependency-free (standard library
`json`/`os` only), matching the rest of the web-interface plumbing (see
`html/lib/dbutil.py`). It is not wired into any current feature -- it
exists as scaffolding for code that will want site-level settings
(`site_name`, `base_url`, and eventually a non-SQLite database backend's
credentials) without needing to handle a missing config file itself.
"""

import json
import os

_PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

WEBCONFIG_PATH = os.path.join(_PROJECT_ROOT, "webconfig.json")
"""str: Absolute path to the real site config file, at the repo root --
a sibling of `html/`, `db/`, and `stellarObjects/`, deliberately outside
`html/`'s Apache `DocumentRoot` the same way `db/` already is. See
`WEBCONFIG.md` for why."""

DEFAULT_WEBCONFIG = {
    "site_name": "planetGen",
    "base_url": "http://localhost/",
    "db_username": "",
    "db_password": "",
    "db_name": "",
}
"""dict: Fallback values, matching `webconfig.json.example`'s shape,
used when no `webconfig.json` file exists yet. `db_username`/
`db_password`/`db_name` are unused placeholders reserved for a possible
future non-SQLite database backend -- see `WEBCONFIG.md`."""


def load_webconfig():
    """
    Loads the site configuration, falling back to defaults if the real
    file doesn't exist -- callers never need to handle a missing-file case
    themselves.

    Returns:
        dict: The parsed contents of `webconfig.json`, or a copy of
              `DEFAULT_WEBCONFIG` if that file isn't present.
    """
    if not os.path.isfile(WEBCONFIG_PATH):
        return dict(DEFAULT_WEBCONFIG)

    with open(WEBCONFIG_PATH, "r", encoding="utf-8") as f:
        return json.load(f)
