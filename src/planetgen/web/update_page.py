# planetgen/web/update_page.py

"""
The Update screen (SEC.43), `/admin/update`: shows whether an update is waiting or running, its exit status and the end
of its log, and lets the Owner press Update. The button does not run anything itself: it asks the API
(`POST /api/admin/update`) to make the request file that a root-owned path unit turns into `update.sh` (see
`planetgen/admin/webupdate.py`). The Owner must have typed the password within the last few minutes; the page asks for
it with the same form the Owner settings use. Works without script: while an update is waiting or running the page
reloads itself every ten seconds.
"""

from flask import request, url_for

from planetgen.web.lib import apiclient

from . import bp
from .admin_pages import _api_message, _cookie_header, _flash, _render, _require_admin, _see_other, _take_flash
from .helpers import crumb

REFRESH_SECONDS = 10
"""int: How often the page reloads itself while an update is waiting or running."""


def _page(data, message=None, error=None, status=200):
    live = bool(data) and data["state"] in ("requested", "running")
    response = _render(
        "admin_update.html", flashed=True, title="Update", section="admin_update", data=data, message=message,
        error=error, refresh=REFRESH_SECONDS if live else None,
        breadcrumbs=[crumb("Admin", "admin"), crumb("Update")],
        description="Run update.sh on the server from this page.")
    response.status_code = status
    return response


@bp.route("/admin/update", methods=["GET", "POST"])
def admin_update():
    """The update's state (GET); POST `action=confirm-password` opens the Owner window, `action=update` asks for the
    update."""
    admin, bounce = _require_admin()
    if bounce is not None:
        return bounce
    cookie_header = _cookie_header()
    target = url_for("web.admin_update")
    if request.method == "POST":
        action = request.form.get("action") or ""
        try:
            if action == "confirm-password":
                done = apiclient.auth_confirm_password(cookie_header, request.form.get("current_password", ""))
                return _flash(_see_other(target), message=f"Thank you. You can start an update for {done['minutes']} minutes.")
            if action == "update":
                apiclient.admin_update_request(cookie_header)
                return _flash(_see_other(target), message="The update was requested. It starts in a moment; this page "
                                                         "follows it.")
        except apiclient.ApiError as exc:
            return _flash(_see_other(target), error=_api_message(exc))
        return _flash(_see_other(target), error="Choose Update, or type the password to unlock it.")
    try:
        data = apiclient.admin_update(cookie_header)
    except apiclient.ApiError as exc:
        return _page(None, error=_api_message(exc))
    flashed = _take_flash()
    return _page(data, message=flashed.get("message"), error=flashed.get("error"))
