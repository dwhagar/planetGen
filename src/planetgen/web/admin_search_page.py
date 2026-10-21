"""
Admin search (ADM.120): `/admin/search`, the plain page behind the shell's search box, and `/admin/palette.json`,
the list the Ctrl+K palette (`static/adminpalette.js`) filters as you type. Both read the registry
(`web/adminregistry.py`), so a page, task or setting is found the day it exists.
"""

from flask import jsonify, request, url_for

from . import adminregistry, bp
from .admin_pages import _render, _require_admin
from .helpers import crumb, page_url

SEARCH_LIMIT = 40


def _address(entry):
    return adminregistry.url(entry, page_url)


@bp.route("/admin/search")
def admin_search():
    """Search every admin page, task and setting by words; works without script."""
    admin, bounce = _require_admin()
    if bounce is not None:
        return bounce
    query = request.args.get("q", "").strip()[:200]
    found = adminregistry.search(query, adminregistry.entries(), limit=SEARCH_LIMIT)
    return _render(
        "admin_search.html", title="Search admin", section="admin",
        breadcrumbs=[crumb("Admin", "admin"), crumb("Search")], query=query,
        results=[{"entry": entry, "url": _address(entry)} for entry in found],
        description="Find an admin page, a task or a setting by what it does.",
    )


@bp.route("/admin/palette.json")
def admin_palette():
    """Every entry for the palette, as JSON, loaded the first time it opens."""
    admin, bounce = _require_admin()
    if bounce is not None:
        return jsonify({"error": "sign in"}), 403
    items = [
        {"key": e.key, "kind": e.kind, "label": e.label, "summary": e.summary, "danger": e.danger,
         "destination": e.destination, "keywords": list(e.keywords), "url": _address(e)}
        for e in adminregistry.entries()
    ]
    response = jsonify({"items": items, "search_url": url_for("web.admin_search")})
    response.headers["Cache-Control"] = "no-store"
    return response


admin_palette.json_only = True  # not a page: tests/test_web_a11y.py skips it
