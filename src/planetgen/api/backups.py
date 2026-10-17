# planetgen/api/backups.py

"""
The backup endpoints behind the Admin backups screen (ADM.85), under `/api/admin/backups`, and the settings export
and import under `/api/admin/settings`. All need a logged-in admin past the forced credential change.

- `GET /api/admin/backups`: the archives newest first, each with its parts, size, who took it and why, and
  whether it holds every part this build knows; `problem` says why the folder can't be used (missing or not
  writable by the web user); `can_download` is whether this caller may download an archive (the Owner only).
- `POST /api/admin/backups` `{}`: Back up now, every part. Audited as `backup.take`.
- `GET /api/admin/backups/<name>`: one archive, whether it verifies, and the diff of its settings part against
  today's `config.json` (what a restore would change; secrets are `(changed)`). A restore that would be refused
  (the backup no longer fits the model) says why in `diff.errors`.
- `POST /api/admin/backups/<name>/restore` `{"parts": [...], "expect_hash": "..."}`: restores the whole archive
  (`parts` null) or the named parts after verifying it. A safety archive of the parts it replaces is taken first;
  each changed option gets a `config.restore` audit row and the restore a `backup.restore` row. Below the command
  line only the options an admin may change are put back from the settings part, so a web restore can never swap
  the database password or the secret key.
- `GET /api/admin/backups/<name>/download`: the archive file. It holds secrets, so it is the Owner's (403 for an
  admin until accounts, 8.2).
- `GET /api/admin/settings/export`: the settings an admin may change as a JSON file (the file's own layout, with
  `config_version`; only options set in config.json), without secrets or command-line-only options. `POST /api/admin/settings/import` `{"document": {...}, "expect_hash"}` checks an
  export the way an edit is checked and, unless `check_only`, saves it as `config.import`.
"""

import json

from flask import Blueprint, Response, g, jsonify

from planetgen import backup
from planetgen.util import settings_writer

from . import settings as config_routes
from .authz import audit, require_admin
from .common import ApiError, get_control_db, require_json_body
from .schemas import BackupRestore, SettingsImport, parse_body

bp = Blueprint("backups", __name__, url_prefix="/api/admin")

OWNER_ONLY = "Downloading a backup is for the Owner: it holds every secret in config.json."


def _labels():
    return {name: part.label for name, part in backup.registered().items()}


def _archive_view(archive):
    labels = _labels()
    return {
        "name": archive.name, "taken_at": archive.taken_at, "taken_by": archive.taken_by, "reason": archive.reason,
        "version": archive.version, "size": archive.size, "complete": archive.complete,
        "readable": archive.readable,
        "parts": [{"name": p.name, "label": labels.get(p.name, p.name), "size": p.size, "files": len(p.files),
                   "known": p.known} for p in archive.parts],
    }


def _can_download():
    return config_routes.ROLE in ("owner", "host")


def _find(name):
    try:
        for archive in backup.list_archives():
            if archive.name == name:
                return archive
    except backup.BackupError as exc:
        raise ApiError(str(exc), status_code=409) from exc
    raise ApiError(f"There is no backup named {name!r}.", status_code=404)


@bp.route("/backups", methods=["GET"])
@require_admin(fresh=True)
def list_backups():
    """`GET /api/admin/backups` -- see the module docstring."""
    try:
        archives = [_archive_view(a) for a in backup.list_archives()]
        error = ""
    except backup.BackupError as exc:
        archives, error = [], str(exc)
    return jsonify({
        "problem": backup.write_problem() or error or "", "archives": archives, "kept": backup.kept(),
        "parts": [{"name": n, "label": label} for n, label in _labels().items()],
        "can_download": _can_download(), "file_hash": settings_writer.file_hash(),
    })


@bp.route("/backups", methods=["POST"])
@require_admin(fresh=True)
def take_backup():
    """`POST /api/admin/backups` -- Back up now: one archive of every part."""
    problem = backup.write_problem()
    if problem:
        raise ApiError(problem, status_code=409)
    try:
        archive = backup.take(who=g.admin_user["username"], reason="manual")
    except backup.BackupError as exc:
        raise ApiError(str(exc), status_code=409) from exc
    audit("backup.take", target=f"backup:{archive.name}", detail=json.dumps({"parts": archive.part_names}))
    return jsonify(_archive_view(archive))


@bp.route("/backups/<name>", methods=["GET"])
@require_admin(fresh=True)
def show_backup(name):
    """`GET /api/admin/backups/<name>` -- one archive, whether it verifies, and the settings diff."""
    archive = _find(name)
    checked = backup.verify(name)
    diff = None
    if checked.ok and settings_writer.PART in archive.part_names:
        try:
            diff = config_routes._view(settings_writer.plan_restore(name, role=config_routes.ROLE))
        except settings_writer.SettingsWriteError as exc:
            diff = {"ok": False, "changes": [], "errors": {"": str(exc)}, "shadowed": {}, "restart": [],
                    "file_hash": settings_writer.file_hash(), "saved": False, "backup": None, "audit_errors": []}
    return jsonify({"archive": _archive_view(archive), "verified": checked.ok, "problems": checked.problems,
                    "diff": diff, "can_download": _can_download(), "file_hash": settings_writer.file_hash()})


@bp.route("/backups/<name>/restore", methods=["POST"])
@require_admin(fresh=True)
def restore_backup(name):
    """`POST /api/admin/backups/<name>/restore` -- the whole archive or the named parts."""
    body = parse_body(BackupRestore, require_json_body())
    archive = _find(name)
    wanted = body.parts if body.parts is not None else archive.part_names
    if not wanted:
        raise ApiError("Choose at least one part to restore.")
    if settings_writer.PART in wanted and body.expect_hash is not None \
            and body.expect_hash != settings_writer.file_hash():
        raise ApiError("The settings were changed by someone else since this page was opened; reload it and "
                       "try again.", status_code=409)
    context = backup.RestoreContext(who=g.admin_user["username"], role=config_routes.ROLE,
                                    audit=config_routes._audit)
    try:
        report = backup.restore(name, body.parts, context)
    except backup.BackupError as exc:
        raise ApiError(str(exc), status_code=409) from exc
    results = {part: config_routes._view(result) for part, result in report.results.items()
               if hasattr(result, "changes")}
    return jsonify({"archive": name, "restored": report.restored, "skipped": report.skipped,
                    "safety": report.safety, "results": results})


@bp.route("/backups/<name>/download", methods=["GET"])
@require_admin(fresh=True)
def download_backup(name):
    """`GET /api/admin/backups/<name>/download` -- the archive; the Owner's."""
    if not _can_download():
        raise ApiError(OWNER_ONLY, status_code=403)
    archive = _find(name)
    with open(archive.path, "rb") as handle:
        data = handle.read()
    audit("backup.download", target=f"backup:{name}")
    return Response(data, mimetype="application/x-bzip2",
                    headers={"Content-Disposition": f'attachment; filename="{name}"'})


@bp.route("/settings/export", methods=["GET"])
@require_admin(fresh=True)
def export_settings():
    """`GET /api/admin/settings/export` -- the settings an admin may change, without secrets."""
    try:
        document = settings_writer.export_document(role=config_routes.ROLE)
    except settings_writer.SettingsWriteError as exc:
        raise ApiError(str(exc), status_code=409) from exc
    return jsonify(document)


@bp.route("/settings/import", methods=["POST"])
@require_admin(fresh=True)
def import_settings():
    """`POST /api/admin/settings/import` -- checks (and, unless `check_only`, saves) an export."""
    body = parse_body(SettingsImport, require_json_body())
    changes, errors = settings_writer.import_changes(body.document, role=config_routes.ROLE)
    if errors:
        return jsonify({"ok": False, "changes": [], "errors": errors, "shadowed": {}, "restart": [],
                        "file_hash": settings_writer.file_hash(), "saved": False, "backup": None,
                        "audit_errors": [], "request": {}})
    try:
        if body.check_only:
            result = settings_writer.plan(changes, role=config_routes.ROLE)
        else:
            result = settings_writer.save(changes, actor=g.admin_user["username"], role=config_routes.ROLE,
                                          audit=config_routes._audit, expect_hash=body.expect_hash,
                                          action="config.import")
    except settings_writer.SettingsWriteError as exc:
        raise ApiError(str(exc), status_code=409) from exc
    view = config_routes._view(result)
    view["request"] = changes
    return jsonify(view)
