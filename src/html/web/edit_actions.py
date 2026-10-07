# html/web/edit_actions.py

"""
The admin Delete and Regenerate buttons (TODO ADM.8) on the sector,
system and phenomenon pages, and the system page's Change class (ADM.6)
and Change star (ADM.7) forms.

Each button is a `<details>` disclosure (the confirm step, no script
needed) holding a small POST form back to the page it is on, carrying
`csrf_field()`, `edit_action` (`regenerate` or `delete`), `edit_target`
(`<kind>:<id>`) and, for a delete or a planet regenerate, an
"also delete its facilities" checkbox (`drop_facilities`). The page hands
the POST to `handle_post`, which calls the API (`api/edits.py`, or the
system `PATCH`/`DELETE`), flashes what happened (with the re-validation's
warnings) and redirects (303), so a reload never repeats it.
"""

import re

from flask import abort, current_app, flash, get_flashed_messages, redirect, request

import apiclient

from planetgen.admin import activity_log

from .helpers import current_admin, db_name, page_url

FLASH_OK = "edit"
FLASH_ERROR = "edit-error"

ACTIONS = ("regenerate", "delete", "class", "star")
"""tuple: The `edit_action` values the pages accept."""

_BODY_PATHS = {"planet": "/planets", "moon": "/moons", "belt": "/belts"}


def can_edit(admin):
    """Whether the buttons show: an admin past the forced credential
    change (the API refuses edits from the installer's default login)."""
    return admin is not None and not admin["must_change_credentials"]


def messages():
    """`[(is_error, text)]`: what the last edit flashed. Every page shows
    them (`base.html`), so a delete's message still shows on the page it
    redirects to. A request with no Flask session cookie (nearly every
    visitor) has none, and its session isn't touched."""
    if current_app.config.get("SESSION_COOKIE_NAME", "session") not in request.cookies:
        return []
    return [(category == FLASH_ERROR, text)
            for category, text in get_flashed_messages(with_categories=True,
                                                       category_filter=[FLASH_OK, FLASH_ERROR])]


def _api_message(exc):
    message = re.sub(r"^planetGen API error \(\d+\): ", "", str(exc))
    if getattr(exc, "status_code", None) in (401, 403):
        return "Not allowed. Change the admin username and password the installer set first."
    if "drop_facilities" in message:
        return message.split(";")[0] + ". Tick “also delete its facilities” to go ahead."
    return message


def class_options(system_id):
    """`{"recommended": {"planet:<id>": [classes]}, "all": [classes]}` for
    the Change class menus, or None when the API won't say."""
    try:
        return apiclient.admin_edit(request.headers.get("Cookie"), db_name(), "GET",
                                    f"/systems/{system_id}/class-options")
    except apiclient.ApiError:
        return None


def _call(kind, target_id, action, drop_facilities):
    """The API call for one edit; returns its answer."""
    cookie_header = request.headers.get("Cookie")
    body = {"drop_facilities": True} if drop_facilities else None
    db = db_name()
    if action == "class":
        choice = request.form.get("planet_class") or ""
        force, _sep, planet_class = choice.rpartition(":")
        return apiclient.admin_edit(cookie_header, db, "POST", f"{_BODY_PATHS[kind]}/{target_id}/class",
                                    {"class": planet_class, "force": force == "force"})
    if action == "star":
        return apiclient.admin_edit(cookie_header, db, "POST", f"/systems/{target_id}/star",
                                    {"star_type": request.form.get("star_type") or "",
                                     "drop_facilities": drop_facilities})
    if kind == "system":
        if action == "delete":
            return apiclient.admin_edit(cookie_header, db, "DELETE", f"/systems/{target_id}")
        return apiclient.admin_edit(cookie_header, db, "PATCH", f"/systems/{target_id}",
                                    {"regenerate": {}, "drop_facilities": drop_facilities})
    if kind == "sector":
        if action == "delete":
            return apiclient.admin_edit(cookie_header, db, "DELETE", f"/sectors/{target_id}/contents")
        return apiclient.admin_edit(cookie_header, db, "POST", f"/sectors/{target_id}/regenerate")
    if kind in _BODY_PATHS:
        base = f"{_BODY_PATHS[kind]}/{target_id}"
    else:
        base = f"/phenomena/{kind}/{target_id}"
    if action == "delete":
        return apiclient.admin_edit(cookie_header, db, "DELETE", base, body)
    return apiclient.admin_edit(cookie_header, db, "POST", f"{base}/regenerate", body)


def _flash_result(kind, action, result):
    summary = result.get("summary")
    if not summary:
        if kind == "system":
            summary = "System regenerated." if action == "regenerate" else "System deleted."
        elif kind == "sector" and action == "delete":
            summary = (f"Sector deleted with {result.get('systems', 0)} system(s) and "
                       f"{result.get('phenomena', 0)} phenomena.")
        elif kind == "sector":
            summary = "Sector regenerated." if result.get("sector_id") else \
                "Sector deleted; its slot is outside the galaxy's outline, so nothing was generated."
        else:
            summary = "Done."
    flash(summary, FLASH_OK)
    if result.get("moved"):
        flash(f"Moved to keep the orbits stable: {', '.join(result['moved'])}.", FLASH_OK)
    if result.get("reclassified"):
        flash(f"Reclassified after moving: {', '.join(result['reclassified'])}.", FLASH_OK)
    if result.get("removed"):
        flash(f"Removed, with no stable orbit left: {', '.join(result['removed'])}.", FLASH_OK)
    for warning in result.get("warnings") or []:
        flash(f"Still unstable: {warning}", FLASH_ERROR)


def handle_post(allowed, here, after_delete=None):
    """
    Runs one edit POST (the CSRF token was already checked app-wide).

    Args:
        allowed (set): The `(kind, id)` targets this page shows buttons
            for; anything else is refused.
        here (str): The page's own URL, to redirect back to.
        after_delete (dict, optional): `{kind: url}` to go to instead after
            deleting the page's own object (it no longer exists).

    Returns:
        A redirect response.
    """
    admin = current_admin()
    if admin is None:
        activity_log.event("AUTHZ", "admin.required", path=request.path)
        abort(403)
    action = request.form.get("edit_action")
    kind, _sep, raw_id = (request.form.get("edit_target") or "").partition(":")
    try:
        target_id = int(raw_id)
    except ValueError:
        target_id = None
    if action not in ACTIONS or (kind, target_id) not in allowed \
            or (action == "class" and kind not in ("planet", "moon")) or (action == "star" and kind != "system"):
        flash("That isn't something on this page.", FLASH_ERROR)
        return redirect(here, code=303)
    try:
        result = _call(kind, target_id, action, request.form.get("drop_facilities") == "1")
    except apiclient.NotFoundError:
        flash("It was already gone.", FLASH_ERROR)
        return redirect((after_delete or {}).get(kind, here), code=303)
    except apiclient.ApiError as exc:
        flash(_api_message(exc), FLASH_ERROR)
        return redirect(here, code=303)
    _flash_result(kind, action, result)
    if kind == "sector" and action == "regenerate":
        if result.get("sector_id"):
            return redirect(page_url("sector", sector_id=result["sector_id"]), code=303)
        return redirect(page_url("sectors"), code=303)
    if action == "delete" and after_delete and kind in after_delete:
        return redirect(after_delete[kind], code=303)
    return redirect(here, code=303)
