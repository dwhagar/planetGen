# planetgen/web/assets.py

"""
The site's own icons and preview picture (ADM.44)
=================================================

An admin uploads a favicon, an Apple touch icon or the default link-preview picture on the Images page; this module
checks the file, stores it under a name made from its content, and says which stored file each setting names. The
file name is `<first 16 hex digits of the SHA-256>.<png|ico|jpg|svg>`, so a served file never changes (the response
says `immutable`) and two uploads of the same picture are one file.

Files live in `web.assets_dir`, else `assets` in the state folder, and are served at `/site-assets/<name>`
(`seo_pages.py`).

`KINDS` says what each upload kind accepts. An SVG is parsed (no DOCTYPE or entities, no scripts, event handlers or
external references) and is served with `Content-Security-Policy: default-src 'none'; sandbox`, so even one that
slipped through can run nothing.

This module imports Flask only to build URLs for `icon_links`.
"""

import hashlib
import io
import os
import re
import threading
import xml.etree.ElementTree as ET

from planetgen.backup import core as backup_core
from planetgen.util import settings as settingsModule

MAX_BYTES = 1024 * 1024
"""int: The largest upload (the API takes a request body of 2 MB, and a picture goes in it as base64)."""

BUILT_IN_OG = {"path": "/static/og-default.png", "mime": "image/png", "width": 1200, "height": 630}
"""dict: The picture a link preview shows when none was uploaded (`src/html/static/og-default.png`)."""

MIME = {"png": "image/png", "ico": "image/x-icon", "jpg": "image/jpeg", "svg": "image/svg+xml"}

KINDS = {
    "favicon_svg": {"setting": ("icons", "favicon_svg"), "title": "Icon (SVG)", "formats": ("svg",),
                    "help": "A square SVG. Browsers draw it at any size and in both themes."},
    "favicon_png": {"setting": ("icons", "favicon_png"), "title": "Icon (PNG or ICO)", "formats": ("png", "ico"),
                    "min": (32, 32), "square": True,
                    "help": "A square PNG (at least 32 by 32) or an ICO, for browsers that do not take an SVG."},
    "apple_touch_icon": {"setting": ("icons", "apple_touch_icon"), "title": "Home-screen icon", "formats": ("png",),
                         "min": (120, 120), "square": True,
                         "help": "A square PNG, 180 by 180 is best; iPhones and iPads use it for a bookmark."},
    "default_image": {"setting": ("og", "default_image"), "title": "Link preview picture",
                      "formats": ("png", "jpg"), "min": (600, 315),
                      "help": "A PNG or JPEG, 1200 by 630 is best (at least 600 by 315), shown when a link is shared."},
}
"""dict: Each upload kind: the setting that names the file, what it accepts."""

_NAME = re.compile(r"^[0-9a-f]{16}\.(png|ico|jpg|svg)$")
_lock = threading.Lock()
_info_cache = {}


class AssetError(ValueError):
    """An upload that can't be used; the message says why, for the person who uploaded it."""


def assets_dir():
    """The folder the uploaded files live in."""
    return settingsModule.get_settings().web.assets_dir or os.path.join(settingsModule.state_dir(), "assets")


def valid_name(name):
    return bool(_NAME.match(name or ""))


def path_of(name):
    """The file for a stored `name`, or `None` for a name that is not one of ours (never a path out of the folder)."""
    return os.path.join(assets_dir(), name) if valid_name(name) else None


def mime_of(name):
    return MIME[name.rsplit(".", 1)[1]]


# ---------------------------------------------------------------------
# Checking an upload
# ---------------------------------------------------------------------

_SVG_BAD_TAGS = {"script", "foreignobject", "iframe", "object", "embed", "audio", "video", "link", "meta", "style"}


def _check_svg(data):
    if b"<!DOCTYPE" in data.upper() or b"<!ENTITY" in data.upper():
        raise AssetError("The SVG has a DOCTYPE or entity definition; save it as plain SVG.")
    try:
        root = ET.fromstring(data)
    except ET.ParseError as exc:
        raise AssetError(f"That is not a readable SVG ({exc}).") from exc
    if root.tag.rsplit("}", 1)[-1].lower() != "svg":
        raise AssetError("That is not an SVG: the top element must be <svg>.")
    for element in root.iter():
        tag = element.tag.rsplit("}", 1)[-1].lower() if isinstance(element.tag, str) else ""
        if tag in _SVG_BAD_TAGS:
            raise AssetError(f"The SVG contains a <{tag}> element, which is not allowed in an icon.")
        for key, value in element.attrib.items():
            attribute = key.rsplit("}", 1)[-1].lower()
            if attribute.startswith("on"):
                raise AssetError("The SVG has an event handler, which is not allowed in an icon.")
            if attribute == "href" and not str(value).startswith("#"):
                raise AssetError("The SVG refers to something outside itself, which is not allowed in an icon.")
            if "url(" in str(value).lower() and "url(#" not in str(value).lower():
                raise AssetError("The SVG refers to something outside itself, which is not allowed in an icon.")
    return None


def _image_info(data):
    """`(format, width, height)` of a PNG, JPEG or ICO, checked by decoding it."""
    from PIL import Image, UnidentifiedImageError
    try:
        with Image.open(io.BytesIO(data)) as image:
            kind = {"PNG": "png", "JPEG": "jpg", "ICO": "ico"}.get(image.format)
            if kind is None:
                raise AssetError(f"That is a {image.format or 'unknown'} file; use a PNG, JPEG or ICO.")
            image.load()
            width, height = image.size
            if kind == "ico":
                width, height = max(image.info.get("sizes", {(width, height)}))
            return kind, width, height
    except AssetError:
        raise
    except (UnidentifiedImageError, OSError, ValueError, SyntaxError) as exc:
        raise AssetError("That file is not a readable picture.") from exc


def check(kind, data):
    """
    Checks `data` (bytes) as an upload of `kind`.

    Returns:
        tuple: `(extension, width, height)`; the size is `(None, None)` for an SVG.

    Raises:
        AssetError: What is wrong, in a sentence.
    """
    spec = KINDS.get(kind)
    if spec is None:
        raise AssetError(f"Unknown kind {kind!r}; use one of {', '.join(KINDS)}.")
    if not data:
        raise AssetError("The file is empty.")
    if len(data) > MAX_BYTES:
        raise AssetError(f"The file is {len(data) // 1024:,} KB; the most allowed is {MAX_BYTES // 1024:,} KB.")
    if "svg" in spec["formats"] and data.lstrip()[:1] == b"<":
        _check_svg(data)
        return "svg", None, None
    if spec["formats"] == ("svg",):
        raise AssetError("This must be an SVG file.")
    if data.lstrip()[:1] == b"<":
        raise AssetError(f"{spec['title']} takes {' or '.join(f.upper() for f in spec['formats'])}, not an SVG.")
    extension, width, height = _image_info(data)
    if extension not in spec["formats"]:
        raise AssetError(f"{spec['title']} takes {' or '.join(f.upper() for f in spec['formats'])}, not {extension.upper()}.")
    minimum = spec.get("min")
    if minimum and (width < minimum[0] or height < minimum[1]):
        raise AssetError(f"The picture is {width} by {height}; it must be at least {minimum[0]} by {minimum[1]}.")
    if spec.get("square") and width != height:
        raise AssetError(f"The picture is {width} by {height}; it must be square.")
    return extension, width, height


def store(kind, data):
    """
    Checks `data` and writes it under its content name.

    Returns:
        str: The stored file's name (what the setting holds).

    Raises:
        AssetError: The upload is not usable, or the folder can't be written (the message says how to fix it).
    """
    extension, _width, _height = check(kind, data)
    name = f"{hashlib.sha256(data).hexdigest()[:16]}.{extension}"
    path = path_of(name)
    with _lock:
        if not os.path.exists(path):
            try:
                backup_core.atomic_write(path, data, mode=0o644)
            except backup_core.BackupError as exc:
                raise AssetError(str(exc)) from exc
    return name


def read(name):
    """`(bytes, mime)` of a stored file, or `None` when there is none by that name."""
    path = path_of(name)
    if path is None:
        return None
    try:
        with open(path, "rb") as handle:
            return handle.read(), mime_of(name)
    except OSError:
        return None


def users_of(name, settings=None):
    """The upload kinds whose setting names the stored file `name`."""
    settings = settings or settingsModule.get_settings()
    return [kind for kind, spec in KINDS.items() if getattr(getattr(settings, spec["setting"][0]), spec["setting"][1]) == name]


def remove(name):
    path = path_of(name)
    if path is None:
        raise AssetError("That is not an uploaded file's name.")
    try:
        os.unlink(path)
    except FileNotFoundError:
        pass
    except OSError as exc:
        raise AssetError(f"{name} can't be removed ({exc.strerror or exc}).") from exc


# ---------------------------------------------------------------------
# What the pages link to
# ---------------------------------------------------------------------

def dimensions(name):
    """`(width, height)` of a stored picture, or `(None, None)` (an SVG, or unreadable); remembered per file."""
    if name.endswith(".svg"):
        return None, None
    with _lock:
        if name in _info_cache:
            return _info_cache[name]
    found = read(name)
    result = (None, None)
    if found:
        try:
            _kind, width, height = _image_info(found[0])
            result = (width, height)
        except AssetError:
            pass
    with _lock:
        _info_cache[name] = result
    return result


def asset_url(name):
    from flask import url_for
    return url_for("web.site_asset", name=name)


def og_image(settings=None):
    """The picture link previews show: `{path, mime, width, height}`, the built-in one when none was uploaded (or
    the uploaded file has gone from the folder)."""
    settings = settings or settingsModule.get_settings()
    name = settings.og.default_image
    if name and path_of(name) and os.path.exists(path_of(name)):
        width, height = dimensions(name)
        return {"path": asset_url(name), "mime": mime_of(name), "width": width, "height": height}
    return dict(BUILT_IN_OG)


def icon_links(settings=None):
    """The `<link>` tags for the site's icons, as `[{rel, href, type, sizes}]`; an icon not uploaded falls back to the
    built-in SVG (the SVG one) or is left out (the others)."""
    from planetgen.web import static_url
    settings = settings or settingsModule.get_settings()
    links = []
    svg = settings.icons.favicon_svg
    if svg and os.path.exists(path_of(svg) or ""):
        links.append({"rel": "icon", "href": asset_url(svg), "type": "image/svg+xml", "sizes": "any"})
    else:
        links.append({"rel": "icon", "href": static_url("favicon.svg"), "type": "image/svg+xml", "sizes": "any"})
    for field, rel in (("favicon_png", "icon"), ("apple_touch_icon", "apple-touch-icon")):
        name = getattr(settings.icons, field)
        if name and os.path.exists(path_of(name) or ""):
            width, height = dimensions(name)
            sizes = f"{width}x{height}" if width else ("any" if name.endswith(".ico") else None)
            links.append({"rel": rel, "href": asset_url(name), "type": mime_of(name), "sizes": sizes})
    return links
