# planetgen/web/settings_page.py

"""
The Admin settings page (ADM.84): `/admin/settings` lists every option of config.json with its help text, default
and where the value in use comes from, and lets an admin change them.

A save is two steps. Posting the form (`action=review`) asks the API for the diff (`/api/admin/settings/plan`):
an option that doesn't fit is shown again with its message and nothing is saved; a clean form shows the old and
new value of each changed option, the ones that need a restart and the ones an environment variable overrides.
Confirming (`action=save`) posts the same changes with the file's hash, so a file changed by someone else in
between is refused instead of overwritten. The save makes a backup archive first (`planetgen.backup`) and writes
an audit row per option (`config.update`).

Only options that differ from what the page showed are sent. An option the admin may not change (command line
only, or the Owner's) is listed read-only with the reason. A secret is never shown: its box is empty, "set" or
"empty" says what is stored, a new value replaces it and "Clear it" removes it. The page works without script.

The review page carries the new values in hidden fields (secrets included, so the confirmed save is the one
reviewed); it is sent `no-store` and protected by the CSRF token like every form.

The options themselves, their help and their limits come from the settings model through the API, so a new
`settings.opt(...)` appears here with no change to this module.
"""

import json
import re

from flask import request, url_for

from planetgen.util import settings_view
from planetgen.web.lib import apiclient

from . import bp
from .admin_pages import _api_message, _cookie_header, _flash, _render, _require_admin, _see_other, _take_flash
from .helpers import crumb
from .limit_sliders import LimitSlider, ladder

FIELD = "opt."
"""str: Prefix of each option's form field; the rest is the dotted option path."""

RESET_FIELD = "reset."
"""str: Prefix of the "use the default" check box."""

CLEAR_FIELD = "clear."
"""str: Prefix of a secret's "clear it" check box."""

RATE_FIELD = "rate."
"""str: Prefix of a rate control's fields: `rate.<option>.count`, `.per` and `.off` (a limit of "N per period")."""

RATE_PERIODS = ("second", "minute", "hour", "day")

_RATE = re.compile(r"^\s*(\d+)\s*(?:per|/)\s*(second|minute|hour|day)s?\s*$", re.IGNORECASE)

DURATION_UNITS = {"seconds": 1, "hours": 3600}


def _text(value):
    """An option's value as the text of its box."""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return str(int(value)) if value.is_integer() and abs(value) < 1e15 else repr(value)
    if value is None:
        return ""
    return str(value)


def _show(value):
    """A value in the review table: a list as its items, a bool as on or off."""
    if isinstance(value, list):
        return ", ".join(str(item) for item in value)
    if isinstance(value, bool):
        return "on" if value else "off"
    return _text(value)


def _number(text, kind):
    """`text` as an int or float, or `None` when it isn't one."""
    try:
        return int(text) if kind == "int" else float(text)
    except ValueError:
        return None


def parse_rate(text):
    """`{"limit", "period"}` of one limit written "30 per minute", `{"limit": "", "period": "minute"}` for no limit,
    or `None` for anything else (several limits, a period the control does not offer), which the page edits as
    text."""
    if not text.strip():
        return {"limit": "", "period": "minute"}
    found = _RATE.match(text)
    return {"limit": found.group(1), "period": found.group(2).lower()} if found else None


def duration_text(value, unit):
    """A duration in words ("6 hours", "15 seconds", "1 day 2 hours"); empty when `value` is not a number."""
    try:
        seconds = float(value) * DURATION_UNITS.get(unit, 1)
    except (TypeError, ValueError):
        return ""
    if seconds < 0:
        return ""
    if seconds == 0:
        return "off"
    parts = []
    for name, size in (("day", 86400), ("hour", 3600), ("minute", 60), ("second", 1)):
        count, seconds = divmod(seconds, size)
        if count >= 1:
            parts.append(f"{int(count)} {name}{'' if int(count) == 1 else 's'}")
        if len(parts) == 2:
            break
    return " ".join(parts) if parts else f"{seconds:g} seconds"


def _lines(text):
    return [line.strip() for line in text.replace("\r", "").split("\n") if line.strip()]


def collect(options, form):
    """
    The changes a submitted form asks for.

    Args:
        options (list): The API's option list (`settings_view.describe`).
        form (Mapping): The submitted fields (`request.form`).

    Returns:
        tuple: `(changes, submitted)`: `changes` is `{option: value or {"reset": true}}` for the options whose
            field differs from what the page showed, `submitted` is `{option: what to show in its field}` for
            every option the form carried, to draw the form again around an error.
    """
    changes, submitted = {}, {}
    for option in options:
        path = option["path"]
        name = FIELD + path
        rate_count = f"{RATE_FIELD}{path}.count"
        if not option["allowed"] or (name not in form and RESET_FIELD + path not in form
                                     and CLEAR_FIELD + path not in form and rate_count not in form):
            continue
        kind = option["kind"]
        in_file = option["source"] == "config.json"
        if option["secret"]:
            submitted[path] = ""
            if form.get(CLEAR_FIELD + path):
                if option["is_set"] or in_file:
                    changes[path] = {"reset": True}
            elif form.get(name):
                changes[path] = form.get(name)
            continue
        if rate_count in form and not form.get(RESET_FIELD + path):
            count = form.get(rate_count, "").strip()
            period = form.get(f"{RATE_FIELD}{path}.per", "minute")
            raw = "" if form.get(f"{RATE_FIELD}{path}.off") or not count else f"{count} per {period}"
            submitted[path] = raw
            if raw != _text(option["value"]):
                changes[path] = raw
            continue
        if form.get(RESET_FIELD + path):
            submitted[path] = _text(option["default"]) if kind != "list" else "\n".join(option["default"] or [])
            if in_file:
                changes[path] = {"reset": True}
            continue
        if kind == "bool":
            wanted = "1" in form.getlist(name)
            submitted[path] = wanted
            if wanted != bool(option["value"]):
                changes[path] = wanted
        elif kind == "list":
            raw = form.get(name, "")
            submitted[path] = raw
            if _lines(raw) != list(option["value"] or []):
                changes[path] = _lines(raw)
        elif kind in ("int", "float"):
            raw = form.get(name, "").strip()
            submitted[path] = raw
            parsed = _number(raw, kind)
            if parsed is None:
                changes[path] = raw                    # the model words the error ("not a valid integer")
            elif parsed != option["value"]:
                changes[path] = parsed
        elif option.get("multiline"):
            raw = form.get(name, "").replace("\r\n", "\n").replace("\r", "\n")
            submitted[path] = raw
            if raw != _text(option["value"]):
                changes[path] = raw
        else:
            raw = form.get(name, "").replace("\r", "").replace("\n", "")
            submitted[path] = raw
            if raw != _text(option["value"]):
                changes[path] = raw
    return changes, submitted


def _shown(option, submitted):
    """What the option's field holds when the page is drawn."""
    if option["path"] in submitted:
        return submitted[option["path"]]
    if option["kind"] == "list":
        return "\n".join(option["value"] or [])
    if option["kind"] == "bool":
        return bool(option["value"])
    return _text(option["value"])


def _titles(view):
    return {option["path"]: option["title"] for option in view["options"]}


def _search_text(option):
    """The lowercase text the search matches: the label, summary, long help, key and group of the option."""
    group = settings_view.GROUPS.get(option["group"], "")
    return " ".join((option["title"], option["summary"], option["description"], option["path"], group)).lower()


def _slider(option):
    """The `LimitSlider` of a `slider` option: logarithmic stops over its `x-slider` range, and 0 when it allows it."""
    low, high = option["slider"]
    stops = ladder(low, high)
    if option["limits"].get("ge") == 0:
        stops = (0.0,) + stops
    return LimitSlider(FIELD + option["path"], option["title"], "log", stops, float(option["default"]))


def _form_page(view, submitted=None, errors=None, message=None, error=None, status=200):
    submitted = submitted or {}
    errors = errors or {}
    show_fixed = request.args.get("show") == "fixed"
    query = request.args.get("q", "").strip().lower()
    for option in (option for group in view["groups"] for option in group["options"]):
        option["shown"] = _shown(option, submitted)
        option["error"] = errors.get(option["path"], "")
        option["default_text"] = ("\n".join(option["default"] or []) if option["kind"] == "list"
                                  else _text(option["default"])) if not option["secret"] else ""
        option["effective_text"] = ("\n".join(option["effective"] or []) if option["kind"] == "list"
                                    else _text(option["effective"])) if not option["secret"] else ""
        option["value_text"] = ("\n".join(option["value"] or []) if option["kind"] == "list"
                                else _text(option["value"])) if not option["secret"] else ""
        option["search_text"] = _search_text(option)
        if option["control"] == "slider" and option["slider"]:
            option["slider_spec"] = _slider(option)
            option["shown"] = _text(option["shown"])
        option["rate"] = parse_rate(str(option["shown"])) if option["control"] == "rate" else None
    shown_count = len(view["options"])
    if query:
        terms = query.split()
        for group in view["groups"]:
            group["options"] = [o for o in group["options"] if all(term in o["search_text"] for term in terms)]
        view = dict(view, groups=[group for group in view["groups"] if group["options"]])
        shown_count = sum(len(group["options"]) for group in view["groups"])
    titles = _titles(view)
    whole_file = errors.get("")
    response = _render(
        "admin_settings.html", flashed=True, title="Settings", section="admin_settings",
        breadcrumbs=[crumb("Admin", "admin"), crumb("Settings")],
        description="Every option of config.json, with its help, default and current value.",
        view=view, errors=[(path, titles.get(path, path), text) for path, text in errors.items() if path],
        whole_file_error=whole_file, message=message, error=error, can_save=not view["problem"],
        file_error=view["file_error"] if view["file_error"] != view["problem"] else "",
        show_fixed=show_fixed, q=request.args.get("q", "").strip(), shown_count=shown_count,
        duration_text=duration_text,
    )
    response.status_code = status
    return response


def _review_page(view, plan, hidden, confirm_url=None, heading="These changes will be saved",
                 button="Save these changes", back_url=None):
    """The diff before a save. `hidden` is the `(name, value)` fields the confirming form posts to `confirm_url`."""
    titles = _titles(view)
    units = {option["path"]: option["unit"] for option in view["options"] if option["unit"]}

    def shown(path, value):
        text = _show(value)
        return f"{text} {units[path]}" if path in units and text != "" and not isinstance(value, bool) else text

    rows = [dict(change, title=titles.get(change["path"], change["path"]), old_text=shown(change["path"], change["old"]),
                 new_text=shown(change["path"], change["new"])) for change in plan["changes"]]
    return _render(
        "admin_settings_review.html", title="Review the changes", section="admin_settings",
        breadcrumbs=[crumb("Admin", "admin"), crumb("Settings", "admin_settings"), crumb("Review")],
        description="The changes a save would make to config.json.",
        rows=rows, shadowed={titles.get(path, path): var for path, var in plan["shadowed"].items()},
        restart=[titles.get(path, path) for path in plan["restart"]],
        hidden=hidden, confirm_url=confirm_url or url_for("web.admin_settings"), heading=heading, button=button,
        back_url=back_url or url_for("web.admin_settings"),
    )


def _saved_message(result, view):
    titles = _titles(view)
    count = len(result["changes"])
    text = f"Saved {count} {'setting' if count == 1 else 'settings'}."
    if result["backup"]:
        text += f" The old file is kept in the backup {result['backup']}."
    restart = [titles.get(path, path) for path in result["restart"]]
    if restart:
        text += " These apply after the web server restarts: " + ", ".join(restart) + "."
    if result["shadowed"]:
        text += (" An environment variable still overrides: "
                 + ", ".join(f"{titles.get(path, path)} ({var})" for path, var in result["shadowed"].items()) + ".")
    if result["audit_errors"]:
        text += " The audit log could not be written for: " + "; ".join(result["audit_errors"]) + "."
    return text


@bp.route("/admin/settings", methods=["GET", "POST"])
def admin_settings():
    """The settings form, its review step and the confirmed save."""
    admin, bounce = _require_admin()
    if bounce is not None:
        return bounce
    cookie_header = _cookie_header()
    try:
        view = apiclient.admin_settings(cookie_header)
    except apiclient.ApiError as exc:
        return _render("admin_settings.html", title="Settings", section="admin_settings",
                       breadcrumbs=[crumb("Admin", "admin"), crumb("Settings")], view=None, errors=[],
                       whole_file_error="", message=None, error=_api_message(exc), can_save=False, file_error="",
                       show_fixed=False, q="", shown_count=0, duration_text=duration_text,
                       description="Every option of config.json, with its help, default and current value.")
    if request.method == "GET":
        flashed = _take_flash()
        return _form_page(view, error=flashed.get("error"), message=flashed.get("message"))

    action = request.form.get("action") or "review"
    if action == "save":
        return _save(cookie_header, view)
    changes, submitted = collect(view["options"], request.form)
    drawn_from = request.form.get("expect_hash")
    if drawn_from and drawn_from != view["file_hash"]:
        return _form_page(view, error="The settings were changed by someone else after this form was opened; "
                                      "the form now shows the current values. Make your change again.", status=409)
    if not changes:
        return _form_page(view, submitted=submitted, message="Nothing to change: every value is as it was.")
    try:
        plan = apiclient.admin_settings_plan(cookie_header, changes)
    except apiclient.ApiError as exc:
        return _form_page(view, submitted=submitted, error=_api_message(exc), status=400)
    if plan["errors"]:
        return _form_page(view, submitted=submitted, errors=plan["errors"], status=422,
                          error="Nothing was saved: fix the settings marked below.")
    if not plan["changes"]:
        return _form_page(view, submitted=submitted, message="Nothing to change: every value is as it was.")
    return _review_page(view, plan, [("action", "save"), ("expect_hash", view["file_hash"]),
                                     ("changes", json.dumps(changes, sort_keys=True))])


def _save(cookie_header, view):
    """The confirmed save: the reviewed changes, with the hash of the file they were reviewed against."""
    try:
        changes = json.loads(request.form.get("changes") or "")
    except ValueError:
        changes = None
    if not isinstance(changes, dict) or not changes:
        return _form_page(view, error="That review is no longer valid; make the changes again.", status=400)
    expect = request.form.get("expect_hash") or None
    try:
        result = apiclient.admin_settings_save(cookie_header, changes, expect_hash=expect)
    except apiclient.ApiError as exc:
        return _form_page(view, error=_api_message(exc), status=409 if getattr(exc, "status_code", 0) == 409 else 400)
    if result["errors"]:
        return _form_page(view, errors=result["errors"], status=422,
                          error="Nothing was saved: fix the settings marked below.")
    if not result["saved"]:
        return _form_page(view, message="Nothing to change: every value is as it was.")
    return _flash(_see_other(url_for("web.admin_settings")), message=_saved_message(result, view))
