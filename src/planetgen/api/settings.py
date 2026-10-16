# planetgen/api/settings.py

"""
The configuration endpoints behind the Admin settings page (ADM.84), under `/api/admin/settings`. All of them
need a logged-in admin past the forced credential change.

- `GET /api/admin/settings`: every option with its help, limits, default, the value in config.json, where the
  value in use comes from, whether this admin may change it and whether a change needs a restart
  (`planetgen.util.settings_view.describe`), the file's `file_hash`, why the file can't be written (`problem`),
  and the options saved since this process started that still need a restart (`pending_restart`).
- `POST /api/admin/settings/plan` `{"changes": {...}}`: the diff a save would make, or the error at each option;
  nothing is written.
- `POST /api/admin/settings` `{"changes": {...}, "expect_hash": "..."}`: the same check, then the save
  (`planetgen.util.settings_writer.save`): a backup archive, an atomic write and one audit row per option. A
  value that doesn't fit is a 200 with `saved: false` and `errors`; a file changed by someone else since the form
  was drawn is a 409; a file that can't be written is a 409 with the reason.

A secret's value is never in an answer: a change to one is shown as `(changed)`, and a stored one as `is_set`.
Writing the file needs the web server's user to own config.json (`examples/apache/set-permissions.sh`).
"""

from flask import Blueprint, g, jsonify

from planetgen.admin import auth as adminAuth
from planetgen.util import settings_view, settings_writer

from .authz import require_admin
from .common import ApiError, get_control_db, require_json_body
from .schemas import SettingsChanges, parse_body

bp = Blueprint("settings", __name__, url_prefix="/api/admin/settings")

ROLE = "admin"
"""str: The writer role of a signed-in admin (`settings_writer.allowed`); the Owner's comes with accounts (8.2)."""


def _changes(body):
    """The writer's `{path: value or RESET}` from a request's `changes`."""
    out = {}
    for path, value in body.changes.items():
        if isinstance(value, dict):
            if value != {"reset": True}:
                raise ApiError(f"'{path}': a change is a value, or {{\"reset\": true}}")
            value = settings_writer.RESET
        out[path] = value
    return out


def _view(plan):
    """A `Plan` or `Result` as JSON: the changes (secrets masked), the errors, the options that need a restart
    and the options an environment variable overrides."""
    return {
        "ok": plan.ok,
        "changes": [{"path": c.path, "old": c.shown("old"), "new": c.shown("new"), "secret": c.secret,
                     "restart": c.restart, "reset": c.reset} for c in plan.changes],
        "errors": plan.errors,
        "shadowed": plan.shadowed,
        "restart": plan.restart,
        "file_hash": plan.file_hash,
        "saved": getattr(plan, "saved", False),
        "backup": getattr(plan, "backup", None),
        "audit_errors": getattr(plan, "audit_errors", []),
    }


def _audit(action, target, detail):
    """One control-DB audit row (the writer writes the activity-log line itself)."""
    admin = g.admin_user
    adminAuth.record_audit(get_control_db(), admin["id"], admin["username"], action, target=target, detail=detail)


@bp.route("", methods=["GET"])
@require_admin(fresh=True)
def get_settings():
    """`GET /api/admin/settings` -- see the module docstring."""
    view = settings_view.describe(role=ROLE)
    view["pending_restart"] = [{"path": path, "title": title} for path, title in settings_view.pending_restart()]
    return jsonify(view)


@bp.route("/plan", methods=["POST"])
@require_admin(fresh=True)
def plan_settings():
    """`POST /api/admin/settings/plan` -- the diff a save would make; nothing is written."""
    body = parse_body(SettingsChanges, require_json_body())
    return jsonify(_view(settings_writer.plan(_changes(body), role=ROLE)))


@bp.route("", methods=["POST"])
@require_admin(fresh=True)
def save_settings():
    """`POST /api/admin/settings` -- checks and saves; see the module docstring."""
    body = parse_body(SettingsChanges, require_json_body())
    try:
        result = settings_writer.save(_changes(body), actor=g.admin_user["username"], role=ROLE, audit=_audit,
                                      expect_hash=body.expect_hash)
    except settings_writer.StaleSettingsError as exc:
        raise ApiError(str(exc), status_code=409) from exc
    except settings_writer.SettingsWriteError as exc:
        raise ApiError(str(exc), status_code=409) from exc
    return jsonify(_view(result))
