# planetgen/web/images_page.py

"""
The Images screen (ADM.44), under the Settings tab: `/admin/settings/images` shows the site icon, the home-screen
icon and the link-preview picture in use, and lets an admin upload a new one or go back to the built-in.

Each picture is checked when it is uploaded (`planetgen.web.assets.check`: the format, the size, square where an
icon must be, a sandboxed SVG), stored under a name made from its content and named by the setting that uses it
(`icons.favicon_svg`, `icons.favicon_png`, `icons.apple_touch_icon`, `og.default_image`), so the change is backed
up and audited like any settings save. The page works without script: the browser posts the file as a form, and
the page hands it to the API as base64.
"""

import base64

from flask import request, url_for

from planetgen.web import assets
from planetgen.web.lib import apiclient

from . import bp
from .admin_pages import _api_message, _cookie_header, _flash, _render, _require_admin, _see_other, _take_flash
from .helpers import crumb


@bp.route("/admin/settings/images", methods=["GET", "POST"])
def admin_images():
    """The pictures in use; POST uploads (`kind` and `file`) or, with `remove`, goes back to the built-in one."""
    admin, bounce = _require_admin()
    if bounce is not None:
        return bounce
    cookie_header = _cookie_header()
    target = url_for("web.admin_images")
    if request.method == "POST":
        kind = request.form.get("kind", "")
        if kind not in assets.KINDS:
            return _flash(_see_other(target), error="Choose which picture to change.")
        title = assets.KINDS[kind]["title"]
        try:
            if request.form.get("remove"):
                apiclient.admin_image_remove(cookie_header, kind)
                return _flash(_see_other(target), message=f"{title}: back to the built-in one.")
            upload = request.files.get("file")
            data = upload.read(assets.MAX_BYTES + 1) if upload else b""
            if not data:
                return _flash(_see_other(target), error=f"{title}: choose a file to upload.")
            if len(data) > assets.MAX_BYTES:
                return _flash(_see_other(target), error=f"{title}: the file is larger than {assets.MAX_BYTES // 1024} KB.")
            apiclient.admin_image_upload(cookie_header, kind, base64.b64encode(data).decode("ascii"))
        except apiclient.ApiError as exc:
            return _flash(_see_other(target), error=f"{title}: {_api_message(exc)}")
        return _flash(_see_other(target), message=f"{title}: saved.")
    try:
        data = apiclient.admin_images(cookie_header)
    except apiclient.ApiError as exc:
        return _render("admin_images.html", title="Images", section="admin_settings", data=None,
                       error=_api_message(exc), message=None, flashed=True,
                       breadcrumbs=[crumb("Admin", "admin"), crumb("Settings", "admin_settings"), crumb("Images")])
    flashed = _take_flash()
    return _render(
        "admin_images.html", flashed=True, title="Images", section="admin_settings",
        breadcrumbs=[crumb("Admin", "admin"), crumb("Settings", "admin_settings"), crumb("Images")],
        description="The site icon, the home-screen icon and the link-preview picture.",
        data=data, message=flashed.get("message"), error=flashed.get("error"),
    )
