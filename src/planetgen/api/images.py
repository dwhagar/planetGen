# planetgen/api/images.py

"""
The site's icons and link-preview picture (ADM.44), under `/api/admin/images`. All of them need a logged-in admin past
the forced credential change.

- `GET /api/admin/images`: for each kind (`favicon_svg`, `favicon_png`, `apple_touch_icon`, `default_image`) what it
  accepts, the stored file now in use (`name`, `url`, `mime`, `width`, `height`; `null` for the built-in one) and the
  built-in preview picture's address.
- `POST /api/admin/images/<kind>` `{"data": "<base64>"}`: checks the picture (`planetgen.web.assets.check`), stores
  it under its content name and saves the setting that names it (`icons.favicon_svg`, ...) through the settings
  writer, so the change is backed up and audited like any other. A picture that cannot be used is a 422 saying why.
- `DELETE /api/admin/images/<kind>`: resets the setting (the built-in picture or no icon applies again) and removes
  the stored file when no other setting names it.

The file is served at `/site-assets/<name>`.
"""

import base64
import binascii
import os

from flask import Blueprint, g, jsonify

from planetgen.admin import auth as adminAuth
from planetgen.util import settings as settingsModule
from planetgen.util import settings_writer
from planetgen.web import assets

from .authz import require_admin
from .common import ApiError, get_control_db, require_json_body
from .schemas import ImageUpload, parse_body

bp = Blueprint("images", __name__, url_prefix="/api/admin/images")

ROLE = "admin"


def _audit(action, target, detail):
    admin = g.admin_user
    adminAuth.record_audit(get_control_db(), admin["id"], admin["username"], action, target=target, detail=detail)


def _kind(kind):
    spec = assets.KINDS.get(kind)
    if spec is None:
        raise ApiError(f"Unknown kind {kind!r}; use one of {', '.join(assets.KINDS)}.", status_code=404)
    return spec


def _current(spec):
    section, field = spec["setting"]
    return getattr(getattr(settingsModule.get_settings(), section), field)


def _describe(kind):
    spec = assets.KINDS[kind]
    name = _current(spec)
    stored = bool(name) and os.path.exists(assets.path_of(name) or "")
    width, height = assets.dimensions(name) if stored else (None, None)
    return {
        "kind": kind, "title": spec["title"], "help": spec["help"], "formats": list(spec["formats"]),
        "setting": ".".join(spec["setting"]), "name": name if stored else None,
        "url": f"/site-assets/{name}" if stored else None, "mime": assets.mime_of(name) if stored else None,
        "width": width, "height": height,
    }


@bp.route("", methods=["GET"])
@require_admin(fresh=True)
def list_images():
    """`GET /api/admin/images` -- see the module docstring."""
    return jsonify({"images": [_describe(kind) for kind in assets.KINDS], "built_in_preview": assets.BUILT_IN_OG["path"],
                    "max_bytes": assets.MAX_BYTES})


@bp.route("/<kind>", methods=["POST"])
@require_admin(fresh=True)
def upload_image(kind):
    """`POST /api/admin/images/<kind>` -- see the module docstring."""
    spec = _kind(kind)
    body = parse_body(ImageUpload, require_json_body())
    try:
        data = base64.b64decode(body.data, validate=True)
    except (binascii.Error, ValueError) as exc:
        raise ApiError("'data' must be the picture as base64 text", status_code=422) from exc
    try:
        name = assets.store(kind, data)
    except assets.AssetError as exc:
        raise ApiError(str(exc), status_code=422) from exc
    path = ".".join(spec["setting"])
    previous = _current(spec)
    try:
        result = settings_writer.save({path: name}, actor=g.admin_user["username"], role=ROLE, audit=_audit)
    except settings_writer.SettingsWriteError as exc:
        raise ApiError(str(exc), status_code=409) from exc
    if not result.ok:
        raise ApiError("; ".join(f"{p}: {m}" for p, m in result.errors.items()), status_code=422)
    _forget(previous, name)
    return jsonify({"image": _describe(kind), "saved": True})


@bp.route("/<kind>", methods=["DELETE"])
@require_admin(fresh=True)
def remove_image(kind):
    """`DELETE /api/admin/images/<kind>` -- see the module docstring."""
    spec = _kind(kind)
    previous = _current(spec)
    try:
        settings_writer.save({".".join(spec["setting"]): settings_writer.RESET}, actor=g.admin_user["username"],
                             role=ROLE, audit=_audit)
    except settings_writer.SettingsWriteError as exc:
        raise ApiError(str(exc), status_code=409) from exc
    _forget(previous, None)
    return jsonify({"image": _describe(kind), "saved": True})


def _forget(previous, keep):
    """Removes the stored file `previous` held when nothing names it any more."""
    if previous and previous != keep and not assets.users_of(previous):
        try:
            assets.remove(previous)
        except assets.AssetError:
            pass   # a file left behind costs a few KB; the change itself was saved
