#!/usr/bin/env python3
# html/index.py

"""
Landing page: renders the deployment's database straight onto `index.py`
itself (by calling into `browse.handler` directly) instead of showing a
"pick a database" table, or redirecting to `browse.py?db=...`, first.

This project deploys as one branded starmap per vhost now (see
`config.json`'s `site_name`/`api_base_url`, `docs/config.md`) -- a picker
whose only real job was choosing among a list that's realistically always
length 1 just added an extra click/page load in front of every single
visit, with no upside for the common case it actually runs in production.
A redirect would have needed `db` right there in the `Location` URL, put
right back in the address bar the same POST-only `db` reaches every other
page with now avoids (see `lib/page.py`'s `post_link`); calling
`browse.handler` in-process instead means this page never needs a `db` of
its own at all -- `browse.py` still works exactly the same reached
directly (its own bottom-of-module `run(handler)` is guarded by
`if __name__ == "__main__":` specifically so importing it here, to reuse
just that one function, doesn't also execute it against *this* request).
`GET /api/databases` still lists every schema matching this deployment's
prefix (see that endpoint's own docs) for a deployment that genuinely keeps
more than one; each one remains reachable directly by choosing it on
`browse.py` (or typing `db` into any of its own forms), just not from this
landing page -- this only ever renders the first one returned
(name-sorted).
"""

import os
import sys
from urllib.parse import quote

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
# stellarObjects/ lives at src/stellarObjects/ (src layout); this file is at
# src/html/, one level down from src/ -- same pattern lib/page.py already
# uses for its own load_config import, so this page's own "no databases"
# title reflects config.json's site_name too, not a hardcoded "planetGen".
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import browse  # noqa: E402 -- html/'s own dir is already sys.path[0] (this script's own), so a plain sibling import
from apiclient import ApiError, NotFoundError, list_databases  # noqa: E402
from page import render, run  # noqa: E402
from stellarObjects.appconfig import load_config  # noqa: E402


def handler():
    body = """
<section class="panel">
<p>No databases found. Generate one first with
<code>sectorGen.py</code> or <code>systemGen.py</code>
(see the project README).</p>
</section>
"""
    return f"{load_config()['site_name']} Databases", body


try:
    _databases = list_databases()
except (ApiError, NotFoundError):
    # The planetGen API is unreachable, or a misrouted/misconfigured
    # backend 404s instead of returning JSON -- fall through to
    # run(handler), whose own ApiError/NotFoundError handling renders a
    # clear error page instead of a raw traceback (this probe would
    # otherwise crash before run() ever gets a chance to catch anything).
    _databases = None

if _databases:
    _db_name = _databases[0]["name"]
    # The sidenav (`lib/page.py`'s `_sidenav_html`) reads `db` from
    # `nav_params()`, which for this GET request means this script's own
    # (empty) QUERY_STRING -- it has no idea `browse.handler` was just
    # called with one. Setting it here (server-side only; the browser
    # never sees this, unlike a redirect's `Location` URL would) is what
    # gets the sidenav's Search/Galaxy/Sectors/etc. links scoped to this
    # database too, exactly as if this request had arrived as
    # `?db=<name>` in the first place.
    os.environ["QUERY_STRING"] = f"db={quote(_db_name)}"
    render(*browse.handler(_db_name))
else:
    run(handler)
