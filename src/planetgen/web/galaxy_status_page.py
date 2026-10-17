# planetgen/web/galaxy_status_page.py

"""
`/galaxy-status` (UX.3): JSON saying whether a background run is changing the galaxy, polled every 30 seconds
by `static/galaxybanner.js` on every page so the banner appears and goes without a reload. Public, and nothing
private in it (`lib/galaxystatus.py`).
"""

from flask import jsonify

from . import bp
from .helpers import db_name
from .lib import galaxystatus

@bp.route("/galaxy-status")
def galaxy_status():
    """`{"active": false}`, or `{"active": true, "hours": 4, "percent": 31.5, "stage": "...", "text": "..."}`."""
    response = jsonify(galaxystatus.status(db_name()))
    response.headers["Cache-Control"] = f"public, max-age={int(galaxystatus.CACHE_SECONDS)}"
    return response


galaxy_status.json_only = True  # not a page: tests/test_web_a11y.py skips it
