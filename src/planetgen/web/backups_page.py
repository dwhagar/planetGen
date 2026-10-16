# planetgen/web/backups_page.py

"""
The backup archives screen and the settings export and import (ADM.85), under the Settings tab.

- `/admin/settings/backups`: the archives (when, who, why, which parts, how big) and Back up now.
- `/admin/settings/backups/<archive>`: one archive: whether it verifies, the diff of its settings part against
  today's `config.json`, and a form that restores the whole archive or only the parts ticked. The diff page is
  the confirmation: the button says exactly what it does, and the restore itself takes a safety archive of what
  it replaces and writes the audit rows.
- `/admin/settings/export` downloads the settings an admin may change as JSON, without secrets or
  command-line-only options. `/admin/settings/import` takes such a file, checks it the way an edit is checked
  (the same validation, the same diff on the review page, a bad value changes nothing) and
  `/admin/settings/import/save` saves the reviewed document.

An archive holds every secret in `config.json`, so downloading one is the Owner's (the link appears for the Owner
from 8.2); an admin sees why it is not offered.
"""

import json

from flask import make_response, request, url_for

from planetgen.web.lib import apiclient

from . import bp
from .admin_pages import _api_message, _cookie_header, _flash, _no_store, _render, _require_admin, _see_other, _take_flash
from .helpers import crumb
from .settings_page import _review_page, _saved_message, _show, _titles

MAX_IMPORT_BYTES = 1024 * 1024
"""int: The largest settings file the import accepts."""


def _backups_page(data, message=None, error=None, import_errors=None, status=200):
    response = _render(
        "admin_backups.html", flashed=True, title="Backups", section="admin_settings",
        breadcrumbs=[crumb("Admin", "admin"), crumb("Settings", "admin_settings"), crumb("Backups")],
        description="The backup archives, Back up now, and the settings export and import.",
        data=data, message=message, error=error, import_errors=import_errors or [],
    )
    response.status_code = status
    return response


@bp.route("/admin/settings/backups", methods=["GET", "POST"])
def admin_backups():
    """The archive list; POST is Back up now."""
    admin, bounce = _require_admin()
    if bounce is not None:
        return bounce
    cookie_header = _cookie_header()
    if request.method == "POST":
        try:
            archive = apiclient.admin_backup_take(cookie_header)
        except apiclient.ApiError as exc:
            return _flash(_see_other(url_for("web.admin_backups")), error=_api_message(exc))
        parts = ", ".join(part["label"] for part in archive["parts"]) or "nothing"
        return _flash(_see_other(url_for("web.admin_backups")),
                      message=f"Backed up {parts} in {archive['name']}.")
    try:
        data = apiclient.admin_backups(cookie_header)
    except apiclient.ApiError as exc:
        return _backups_page(None, error=_api_message(exc))
    flashed = _take_flash()
    return _backups_page(data, message=flashed.get("message"), error=flashed.get("error"))


def _detail_page(data, name, message=None, error=None, status=200):
    diff = data["diff"] or {}
    titles = {}
    try:
        titles = _titles(apiclient.admin_settings(_cookie_header()))
    except apiclient.ApiError:
        pass
    rows = [dict(change, title=titles.get(change["path"], change["path"]), old_text=_show(change["old"]),
                 new_text=_show(change["new"])) for change in diff.get("changes", [])]
    response = _render(
        "admin_backup.html", flashed=True, title=name, section="admin_settings",
        breadcrumbs=[crumb("Admin", "admin"), crumb("Settings", "admin_settings"),
                     crumb("Backups", "admin_backups"), crumb("Archive")],
        description="One backup archive and what restoring it would change.",
        data=data, rows=rows, diff=diff, message=message, error=error,
        restart=[titles.get(path, path) for path in diff.get("restart", [])],
        shadowed={titles.get(path, path): var for path, var in (diff.get("shadowed") or {}).items()},
    )
    response.status_code = status
    return response


@bp.route("/admin/settings/backups/<archive>")
def admin_backup(archive):
    """One archive, its check and the settings diff."""
    admin, bounce = _require_admin()
    if bounce is not None:
        return bounce
    try:
        data = apiclient.admin_backup(_cookie_header(), archive)
    except apiclient.NotFoundError:
        return _flash(_see_other(url_for("web.admin_backups")), error="There is no such backup.")
    except apiclient.ApiError as exc:
        if getattr(exc, "status_code", 0) == 404:
            return _flash(_see_other(url_for("web.admin_backups")), error="There is no such backup.")
        return _flash(_see_other(url_for("web.admin_backups")), error=_api_message(exc))
    flashed = _take_flash()
    return _detail_page(data, archive, message=flashed.get("message"), error=flashed.get("error"))


@bp.route("/admin/settings/backups/<archive>/restore", methods=["POST"])
def admin_backup_restore(archive):
    """Restores the ticked parts (all when none are named in the form's `all` field)."""
    admin, bounce = _require_admin()
    if bounce is not None:
        return bounce
    cookie_header = _cookie_header()
    parts = request.form.getlist("part")
    target = url_for("web.admin_backup", archive=archive)
    if not parts:
        return _flash(_see_other(target), error="Tick at least one part to restore.")
    try:
        result = apiclient.admin_backup_restore(cookie_header, archive, parts=parts,
                                                expect_hash=request.form.get("expect_hash") or None)
    except apiclient.ApiError as exc:
        return _flash(_see_other(target), error=_api_message(exc))
    text = f"Restored {', '.join(result['restored']) or 'nothing'} from {archive}."
    settings_result = result["results"].get("settings")
    if settings_result is not None:
        count = len(settings_result["changes"])
        text += f" {count} {'setting' if count == 1 else 'settings'} changed."
        if settings_result["restart"]:
            text += " These apply after the web server restarts: " + ", ".join(settings_result["restart"]) + "."
    if result["skipped"]:
        text += " This version of planetGen cannot restore, and so skipped: " + ", ".join(result["skipped"]) + "."
    if result["safety"]:
        text += f" What was replaced is kept in the backup {result['safety']}."
    return _flash(_see_other(target), message=text)


@bp.route("/admin/settings/export")
def admin_settings_export():
    """The settings an admin may change as a download."""
    admin, bounce = _require_admin()
    if bounce is not None:
        return bounce
    try:
        document = apiclient.admin_settings_export(_cookie_header())
    except apiclient.ApiError as exc:
        return _flash(_see_other(url_for("web.admin_backups")), error=_api_message(exc))
    response = make_response(json.dumps(document, indent=2) + "\n")
    response.headers["Content-Type"] = "application/json"
    response.headers["Content-Disposition"] = 'attachment; filename="planetgen-settings.json"'
    return _no_store(response)


admin_settings_export.json_only = True   # a download, not a page: the accessibility sweep leaves it out


@bp.route("/admin/settings/import", methods=["POST"])
def admin_settings_import():
    """Checks an uploaded export and shows the diff, or the errors."""
    admin, bounce = _require_admin()
    if bounce is not None:
        return bounce
    cookie_header = _cookie_header()
    upload = request.files.get("file")
    try:
        data = apiclient.admin_backups(cookie_header)
    except apiclient.ApiError as exc:
        return _backups_page(None, error=_api_message(exc))
    if upload is None or not upload.filename:
        return _backups_page(data, error="Choose a settings file to import.", status=400)
    raw = upload.read(MAX_IMPORT_BYTES + 1)
    if len(raw) > MAX_IMPORT_BYTES:
        return _backups_page(data, error="That file is too big to be a settings export.", status=400)
    try:
        document = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, ValueError):
        return _backups_page(data, error="That file is not valid JSON; nothing was imported.", status=400)
    if not isinstance(document, dict):
        return _backups_page(data, error="A settings file is a JSON object; nothing was imported.", status=400)
    try:
        plan = apiclient.admin_settings_import(cookie_header, document, check_only=True)
        view = apiclient.admin_settings(cookie_header)
    except apiclient.ApiError as exc:
        return _backups_page(data, error=_api_message(exc), status=400)
    titles = _titles(view)
    if plan["errors"]:
        problems = [(titles.get(path, path) if path else "The file", text) for path, text in plan["errors"].items()]
        return _backups_page(data, error="Nothing was imported: the file has problems.", import_errors=problems,
                             status=422)
    if not plan["changes"]:
        return _backups_page(data, message="Nothing to import: every value in the file is as it is.")
    return _review_page(
        view, plan, [("expect_hash", plan["file_hash"]), ("document", json.dumps(document, sort_keys=True))],
        confirm_url=url_for("web.admin_settings_import_save"), heading="This import will change these settings",
        button="Import these settings", back_url=url_for("web.admin_backups"))


@bp.route("/admin/settings/import/save", methods=["POST"])
def admin_settings_import_save():
    """Saves the reviewed import."""
    admin, bounce = _require_admin()
    if bounce is not None:
        return bounce
    cookie_header = _cookie_header()
    back = url_for("web.admin_backups")
    try:
        document = json.loads(request.form.get("document") or "")
    except ValueError:
        document = None
    if not isinstance(document, dict):
        return _flash(_see_other(back), error="That review is no longer valid; import the file again.")
    try:
        result = apiclient.admin_settings_import(cookie_header, document,
                                                 expect_hash=request.form.get("expect_hash") or None)
        view = apiclient.admin_settings(cookie_header)
    except apiclient.ApiError as exc:
        return _flash(_see_other(back), error=_api_message(exc))
    if result["errors"]:
        return _flash(_see_other(back), error="Nothing was imported: " + "; ".join(result["errors"].values()))
    if not result["saved"]:
        return _flash(_see_other(back), message="Nothing to import: every value in the file is as it is.")
    return _flash(_see_other(url_for("web.admin_settings")),
                  message=_saved_message(result, view).replace("Saved", "Imported", 1))
