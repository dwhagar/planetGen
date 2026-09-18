#!/usr/bin/env python3
# html/index.py

"""
Landing page: redirects straight into the deployment's database
(`browse.py`) instead of showing a "pick a database" table first.

This project deploys as one branded starmap per vhost now (see
`config.json`'s `site_name`/`api_base_url`, `docs/config.md`) -- a picker
whose only real job was choosing among a list that's realistically always
length 1 just added an extra click/page load in front of every single
visit, with no upside for the common case it actually runs in production.
`GET /api/databases` still lists every schema matching this deployment's
prefix (see that endpoint's own docs) for a deployment that genuinely keeps
more than one; each one remains reachable directly at
`browse.py?db=<name>`, just no longer from a shared landing page -- this
only ever jumps to the first one returned (name-sorted).
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

from apiclient import ApiError, NotFoundError, list_databases
from page import redirect, run
from stellarObjects.appconfig import load_config


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
    redirect(f"browse.py?db={quote(_databases[0]['name'])}")
else:
    run(handler)
