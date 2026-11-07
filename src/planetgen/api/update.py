# planetgen/api/update.py

"""
The Update button's API (SEC.43), under `/api/admin/update`. Both routes need a signed-in admin past the forced
credential change; an API key is not enough for the request (`session_only`), so a leaked key cannot start an update.

- `GET /api/admin/update`: the state of the machinery (`planetgen.admin.webupdate.view`): `state`
  (`not_installed`, `idle`, `requested`, `running`), `can_request` and `reason`, and the end of `update.log`. `can_ask` says whether this caller may press the button and
  `needs_password` whether only the typed password is missing.
- `POST /api/admin/update` `{}`: asks root to run `update.sh`. Only the Owner whose password was typed within the
  last `adminAuth.RECENT_PASSWORD_MINUTES` may (403 otherwise; `POST /api/auth/confirm-password` opens the window). The
  request is one empty file in the state folder; no value from the request reaches the command, and the body takes no
  field. A press while an update runs is accepted and runs after it; a second one while a request waits is a 409.
"""

from flask import Blueprint, g, jsonify

from planetgen.admin import activity_log, webupdate
from planetgen.admin import auth as adminAuth

from . import settings as config_routes
from .authz import audit, require_admin
from .common import ApiError, require_json_body
from .limiter import limiter
from .schemas import UpdateRequest, parse_body

bp = Blueprint("update", __name__, url_prefix="/api/admin/update")

REQUEST_RATE_LIMIT = "6 per hour"
"""str: The request route's own limit, on top of the app-wide default."""

OWNER_ONLY = "Updating the site is for the Owner, who must have typed the password in the last few minutes."


def _may_ask():
    """`(may, needs_password)` for the caller: the Owner signed in with a session (not an API key), and whether the
    typed password is the only thing missing."""
    owner = g.api_key_id is None and adminAuth.is_owner(g.admin_user)
    return owner and config_routes.writer_role() == "owner", config_routes.owner_needs_password()


@bp.route("", methods=["GET"])
@require_admin(fresh=True, session_only=True)
def get_update():
    """`GET /api/admin/update` -- see the module docstring."""
    may, needs_password = _may_ask()
    return jsonify({**webupdate.view(), "can_ask": may, "needs_password": needs_password,
                    "owner_minutes": adminAuth.RECENT_PASSWORD_MINUTES})


@bp.route("", methods=["POST"])
@limiter.limit(REQUEST_RATE_LIMIT)
@require_admin(fresh=True, session_only=True)
def post_update():
    """`POST /api/admin/update` -- ask root to run `update.sh`; see the module docstring."""
    parse_body(UpdateRequest, require_json_body())
    may, needs_password = _may_ask()
    if not may:
        activity_log.event("AUTHZ", "update.refused", user=g.admin_user["username"], needs_password=needs_password)
        raise ApiError(OWNER_ONLY, status_code=403, extra={"needs_password": needs_password})
    try:
        webupdate.request()
    except webupdate.UpdateRequestError as exc:
        raise ApiError(str(exc), status_code=409) from exc
    audit("update.request", target="site")
    return jsonify({**webupdate.view(), "can_ask": False, "needs_password": False}), 202
