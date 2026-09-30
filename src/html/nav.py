#!/usr/bin/env python3
# html/nav.py

"""
Moved: the NAV page is now served by the Flask app at `/nav`
(`web/nav_page.py`). This shim 301-redirects an old `nav.py` link,
bookmark or picker form post there, translating the old parameters (GET
query or POST body) to the new scheme:

    from=12                                    -> from=system:12
    from=3&from_kind=phenomenon&from_type=nebula -> from=nebula:3
    (and the same for to), plus from_sector / to_sector as they are.

The `db` parameter is dropped: the Flask pages take the database from
config. Removed in the cleanup PR.
"""

import os
import re
import sys
from urllib.parse import urlencode

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))

from page import nav_params, redirect, start_request_log  # noqa: E402

_KIND = re.compile(r"^[a-z_]+$")


def new_query(params):
    """The `/nav` query (a list of pairs) for the old parameters."""
    query = []
    for prefix in ("from", "to"):
        raw = params.get(prefix, "").strip()
        if not raw.isdigit():
            continue
        kind = "system"
        phenomenon_type = params.get(f"{prefix}_type", "")
        if params.get(f"{prefix}_kind") == "phenomenon" and _KIND.match(phenomenon_type):
            kind = phenomenon_type
        query.append((prefix, f"{kind}:{int(raw)}"))
    for name in ("from_sector", "to_sector"):
        raw = params.get(name, "").strip()
        if raw.isdigit():
            query.append((name, raw))
    return query


if __name__ == "__main__":
    start_request_log()
    query = new_query(nav_params())
    redirect(f"/nav?{urlencode(query, safe=':')}" if query else "/nav", status="301 Moved Permanently")
