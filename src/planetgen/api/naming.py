# planetgen/api/naming.py

"""
Names from IDs in the API's answers (GEN.71).

The objects the phoneme codec names (`names/naming_key.py`) keep their
19-digit hex ID in their `name` column. Every JSON answer passes through
`NamingJSONProvider`, which swaps such a `name` for the codec name under
the galaxy's naming key (control database, `galaxy_naming`); the pages
call the API in-process, so they show the same names. Nothing is
rewritten in the database, so changing the key renames at once (the tile
cache follows because `query.galaxy_content_state` mixes the key into its
`base`).
"""

import time

from flask import current_app, g, has_request_context, request
from flask.json.provider import DefaultJSONProvider

from planetgen.api import ids
from planetgen.api.common import get_control_db
from planetgen.names import naming_key

CACHE_SECONDS = 5.0
"""float: How long one process keeps a database's key before asking again."""

_cache = {}
"""dict: database -> (expires at, key or `None`)."""


def forget(database=None):
    """Drops the remembered key of `database` (all, when `None`): an admin
    just changed it."""
    if database is None:
        _cache.clear()
    else:
        _cache.pop(database, None)


def _database():
    return request.args.get("db") or current_app.config["MYSQL_CONFIG"].database


def request_key():
    """The naming key for this request's galaxy database, or `None` before
    one is drawn, when the control database can't be read, or outside a
    request."""
    if not has_request_context():
        return None
    if "naming_key" in g:
        return g.naming_key
    database = _database()
    cached = _cache.get(database)
    if cached is not None and cached[0] > time.monotonic():
        key = cached[1]
    else:
        try:
            key = naming_key.key_of(get_control_db(), database)
        except Exception:  # noqa: BLE001 -- no control database, or a schema older than v9: IDs stay as they are
            key = None
        _cache[database] = (time.monotonic() + CACHE_SECONDS, key)
    g.naming_key = key
    return key


class NamingJSONProvider(DefaultJSONProvider):
    """Flask's JSON provider, naming codec-named objects on the way out."""

    def dumps(self, obj, **kwargs):
        if has_request_context():
            obj = ids.translate_response(obj)
            obj = naming_key.rename_names(obj, request_key)
        return super().dumps(obj, **kwargs)


def install(app):
    """Makes `app` answer with codec names and the rest of the code ask for the key in force."""
    app.json = NamingJSONProvider(app)
    naming_key.set_resolver(request_key)
