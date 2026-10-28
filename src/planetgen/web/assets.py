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

An account's picture (`prepare_avatar`, SEC.36) is the one place a GIF, animated or not, is accepted. It has its own
format list and its own caps, and what is stored is a fresh encoding made here, never the uploaded bytes.

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

MIME = {"png": "image/png", "ico": "image/x-icon", "jpg": "image/jpeg", "svg": "image/svg+xml", "gif": "image/gif"}

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

_NAME = re.compile(r"^[0-9a-f]{16}\.(png|ico|jpg|svg|gif)$")
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


PICTURE_FORMATS = {"PNG": "png", "JPEG": "jpg", "ICO": "ico"}
"""dict: The picture formats an upload may be (Pillow's name -> the file kind). Pillow is only asked to try these (SEC.35)."""


def _image_info(data):
    """`(format, width, height)` of a PNG, JPEG or ICO, checked by decoding it. Pillow tries no other decoder on the
    upload (TIFF, GIF, PPM and the rest are not read at all), so a file of another format is "not a readable picture"."""
    from PIL import Image, UnidentifiedImageError
    try:
        with Image.open(io.BytesIO(data), formats=list(PICTURE_FORMATS)) as image:
            kind = PICTURE_FORMATS[image.format]
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


# ---------------------------------------------------------------------
# An account's picture (SEC.36)
# ---------------------------------------------------------------------

AVATAR_FORMATS = {"PNG": "png", "JPEG": "jpg", "ICO": "ico", "GIF": "gif"}
"""dict: The formats an account picture may be uploaded as (Pillow's name -> the kind). Pillow is only asked to try
these. An ICO is stored as a PNG; a PNG or GIF that is animated stays animated."""

AVATAR_MAX_PIXELS = 1024 * 1024
"""int: The most width times height one frame of an account picture may have (1024 by 1024)."""

AVATAR_MAX_FRAMES = 100
"""int: The most frames an animated account picture may have."""

AVATAR_MAX_TOTAL_PIXELS = 16_000_000
"""int: The most width times height times frames decoded for one upload (about 64 MB held as RGBA)."""

PILLOW_MAX_PIXELS = 25_000_000
"""int: What Pillow's own decompression-bomb limit is set to for this process (its default is 89 million)."""


def _pillow():
    """Pillow's `Image`, with the process-wide bomb limit lowered to `PILLOW_MAX_PIXELS` (never raised)."""
    from PIL import Image
    if Image.MAX_IMAGE_PIXELS is None or Image.MAX_IMAGE_PIXELS > PILLOW_MAX_PIXELS:
        Image.MAX_IMAGE_PIXELS = PILLOW_MAX_PIXELS
    return Image


def _avatar_frames(image, kind):
    """`(frames, durations, loop)` read from an opened picture, frame by frame: the caps are tested before each frame
    is decoded, so a picture over them is refused having decoded at most the frames that fit. Every frame is a copy
    with no `info` (comments, XMP, ICC profile, EXIF and application extensions stay behind)."""
    from PIL import ImageOps
    _pillow()
    width, height = image.size
    if width < 1 or height < 1 or width * height > AVATAR_MAX_PIXELS:
        raise AssetError(f"The picture is {width} by {height}; a picture may have at most {AVATAR_MAX_PIXELS:,} pixels "
                         f"(1024 by 1024).")
    mode = "RGB" if kind == "jpg" else "RGBA"
    frames, durations, total = [], [], 0
    loop = None
    for index in range(AVATAR_MAX_FRAMES + 1):
        try:
            image.seek(index)
        except EOFError:
            break
        if index >= AVATAR_MAX_FRAMES:
            raise AssetError(f"The picture has more than {AVATAR_MAX_FRAMES} frames; the most allowed is {AVATAR_MAX_FRAMES}.")
        frame_width, frame_height = image.size
        total += frame_width * frame_height
        if total > AVATAR_MAX_TOTAL_PIXELS:
            raise AssetError(f"The picture is too large to unpack: {width} by {height} with {index + 1} or more frames "
                             f"is over {AVATAR_MAX_TOTAL_PIXELS:,} pixels in all.")
        if index == 0:
            loop = image.info.get("loop")
        durations.append(int(image.info.get("duration") or 0))
        source = ImageOps.exif_transpose(image) if kind == "jpg" else image   # the turn is kept, the EXIF is not
        frame = source.convert(mode)
        frame.info = {}
        frames.append(frame)
    if not frames:
        raise AssetError("That file is not a readable picture.")
    return frames, durations, loop


def prepare_avatar(data):
    """
    Checks an account picture and makes the file to store.

    The picture must be a PNG, JPEG, ICO or GIF of at most `MAX_BYTES`, `AVATAR_MAX_PIXELS` per frame,
    `AVATAR_MAX_FRAMES` frames and `AVATAR_MAX_TOTAL_PIXELS` pixels over all the frames. The size comes from the header
    before anything is decoded, Pillow is told which formats to try, a decompression-bomb warning is an error, and the
    frames are read one by one under the caps. A GIF or PNG with several frames keeps them, their durations and its loop
    count; nothing else the upload carried (comments, XMP, ICC profiles, EXIF, application extensions) is kept, because
    the stored bytes are Pillow's new encoding of the pixels, not the uploaded ones.

    Returns:
        tuple: `(extension, bytes, width, height, frames)` -- `extension` is `png`, `jpg` or `gif`.

    Raises:
        AssetError: What is wrong, in a sentence.
    """
    import warnings
    from PIL import UnidentifiedImageError
    Image = _pillow()
    if not data:
        raise AssetError("The file is empty.")
    if len(data) > MAX_BYTES:
        raise AssetError(f"The file is {len(data) // 1024:,} KB; the most allowed is {MAX_BYTES // 1024:,} KB.")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", Image.DecompressionBombWarning)
            with Image.open(io.BytesIO(data), formats=list(AVATAR_FORMATS)) as image:
                kind = AVATAR_FORMATS[image.format]
                frames, durations, loop = _avatar_frames(image, kind)
                if kind == "ico":
                    kind = "png"
            width, height = frames[0].size
            out = io.BytesIO()
            options = {}
            if len(frames) > 1 and kind in ("gif", "png"):
                options = {"save_all": True, "append_images": frames[1:], "duration": durations, "disposal": 2}
                if loop is not None:
                    options["loop"] = loop
            elif kind == "jpg":
                options = {"quality": 90}
            frames[0].save(out, {"png": "PNG", "jpg": "JPEG", "gif": "GIF"}[kind], **options)
    except AssetError:
        raise
    except (Image.DecompressionBombWarning, Image.DecompressionBombError) as exc:
        raise AssetError("The picture is far larger than a picture may be (a decompression bomb).") from exc
    except (UnidentifiedImageError, OSError, ValueError, SyntaxError) as exc:
        raise AssetError("That file is not a readable picture (use a PNG, JPEG, ICO or GIF).") from exc
    stored = out.getvalue()
    if len(stored) > MAX_BYTES:
        raise AssetError(f"The picture is {len(stored) // 1024:,} KB once cleaned; the most allowed is {MAX_BYTES // 1024:,} KB.")
    return kind, stored, width, height, len(frames)


def store_avatar(data):
    """
    Checks an account picture (`prepare_avatar`) and writes the cleaned file under its content name.

    Returns:
        dict: `{"name", "width", "height", "frames", "mime"}`.

    Raises:
        AssetError: The picture is not usable, or the folder can't be written.
    """
    extension, stored, width, height, frames = prepare_avatar(data)
    name = f"{hashlib.sha256(stored).hexdigest()[:16]}.{extension}"
    path = path_of(name)
    with _lock:
        if not os.path.exists(path):
            try:
                backup_core.atomic_write(path, stored, mode=0o644)
            except backup_core.BackupError as exc:
                raise AssetError(str(exc)) from exc
    return {"name": name, "width": width, "height": height, "frames": frames, "mime": MIME[extension]}


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
