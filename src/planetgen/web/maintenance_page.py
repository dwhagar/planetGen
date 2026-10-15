# planetgen/web/maintenance_page.py

"""
`/maintenance` (UX.112): what a visitor sees while maintenance mode is on
(`admin/maintenance.py`, `web/maintenance_gate.py`). With it off the page
sends the visitor home, so a browser waiting on the page lands on the site
by itself once the work is done. An admin sees the same page (status 200)
with a link to the admin page's switch.
"""

from flask import redirect, url_for

from planetgen.admin import maintenance

from . import bp
from .maintenance_gate import is_admin, page_response


@bp.route("/maintenance")
def maintenance_page():
    """The "being updated" page, or a redirect home when maintenance is off."""
    state = maintenance.read()
    if not state["enabled"]:
        return redirect(url_for("web.index"), code=302)
    return page_response(state, status=200 if is_admin() else 503)
