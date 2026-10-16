# planetgen/web/seo_pages.py

"""
The files search engines, browsers and link previews ask for (ADM.44):

- `/robots.txt`, `/sitemap.xml` (an index), `/sitemap-pages.xml` and, when the system, sector and phenomenon pages
  may be indexed, `/sitemap-systems-<n>.xml` and `/sitemap-sectors-<n>.xml` (`seo.sitemap_max_urls` addresses a
  file, built from the database in pages of 500 and kept for `seo.sitemap_cache_hours`);
- `/.well-known/security.txt` and `/humans.txt` (the text from the settings; 404 while empty);
- `/site.webmanifest` (while `icons.manifest_enabled`);
- `/site-assets/<name>`: an uploaded icon or picture (`assets.py`), cached for a year since its name is its content;
- `/og/<kind>/<id>.png`: the link-preview picture of one system or sector, drawn by a renderer registered with
  `register_card_renderer` (the card drawer arrives with VIEW.3); with none, or with `og.generated_cards` off,
  it sends the visitor to the default picture.

What each page's `<head>` says is in `seo.py`.
"""

import hashlib
import os
import threading

from flask import Response, abort, redirect, request, url_for

from planetgen.util import settings as settingsModule
from planetgen.web import assets, seo
from planetgen.web.lib import apiclient

from . import bp
from .helpers import db_name

PAGE_SIZE = 500
"""int: The API's largest page, which the sitemaps are built from."""

CARD_RENDERERS = {}
"""dict: Kind (`"system"`, `"sector"`) -> `renderer(object_id) -> png bytes or None`."""


def register_card_renderer(kind, renderer):
    """Installs the function that draws the preview card of a `kind` of object (`"system"` or `"sector"`)."""
    CARD_RENDERERS[kind] = renderer


bp.after_request(seo.response_headers)


def _file(view):
    """Marks a view as a generated file, not a page (`json_only`: the accessibility sweep skips it)."""
    view.json_only = True
    return view


def _text(body, mimetype="text/plain", cache=3600):
    response = Response(body, mimetype=mimetype)
    response.headers["Cache-Control"] = f"public, max-age={cache}"
    return response


@bp.route("/robots.txt")
@_file
def robots_txt():
    """The generated robots.txt (or the override)."""
    return _text(seo.robots_txt())


def _sitemaps_on():
    settings = settingsModule.get_settings()
    if not settings.seo.sitemap_enabled or settings.seo.indexing == "off" or seo.base_url() is None:
        abort(404)
    return settings


def _chunks(total, size):
    return max(1, -(-int(total) // size))


@bp.route("/sitemap.xml")
@_file
def sitemap_index():
    """The sitemap index: the top pages, then the detail pages' files when those may be indexed."""
    settings = _sitemaps_on()

    def build():
        clear_boundaries()
        files = [seo.absolute(url_for("web.sitemap_pages"))]
        if seo.detail_indexed(settings):
            db = db_name()
            size = settings.seo.sitemap_max_urls
            for kind, total in (("systems", apiclient.get_systems(db, limit=1)["total"]),
                                ("sectors", apiclient.get_sectors(db, limit=1)["total"])):
                files += [seo.absolute(url_for("web.sitemap_detail", kind=kind, n=n))
                          for n in range(1, _chunks(total, size) + 1)]
        return seo.sitemap_index_xml(files)

    return _text(seo.cached("index", build), "application/xml")


@bp.route("/sitemap-pages.xml")
@_file
def sitemap_pages():
    """The pages that need no id: home, the lists, classes, species, polities, the galaxy."""
    _sitemaps_on()
    return _text(seo.cached("pages", lambda: seo.sitemap_urls_xml(
        [seo.absolute(path) for path in seo.top_page_paths()])), "application/xml")


_boundaries = {}
_boundary_lock = threading.Lock()


def clear_boundaries():
    """Forgets where each system file ended (rebuilt with the index, so the files tile the list)."""
    with _boundary_lock:
        _boundaries.clear()


def clear_caches():
    """Forgets every built sitemap (tests, and after a settings save)."""
    seo.clear_cache()
    clear_boundaries()


def _system_chunk(db, n, size):
    """The addresses of systems `(n - 1) * size` to `n * size`, walked with `after` from the nearest chunk boundary
    already known (each boundary is the ID the chunk before ended on)."""
    with _boundary_lock:
        known = {k: v for k, v in _boundaries.get(db, {}).items() if k < n}
    start = max(known, default=0)
    after = known.get(start)
    if start and after is None:
        return []                                    # the walk already ran off the end
    urls = []
    for chunk in range(start + 1, n + 1):
        urls, last = [], None
        while len(urls) < size:
            page = apiclient.get_systems(db, limit=min(PAGE_SIZE, size - len(urls)), after=after)["items"]
            if not page:
                break
            urls += [seo.absolute(url_for("web.system", system_id=row["id"])) for row in page]
            last = after = page[-1]["id"]
        with _boundary_lock:
            _boundaries.setdefault(db, {})[chunk] = last
        if last is None:
            return urls if chunk == n else []
    return urls


def _sector_chunk(db, n, size):
    urls, offset = [], (n - 1) * size
    while len(urls) < size:
        page = apiclient.get_sectors(db, limit=min(PAGE_SIZE, size - len(urls)), offset=offset + len(urls))["items"]
        if not page:
            break
        urls += [seo.absolute(url_for("web.sector", sector_id=row["id"])) for row in page]
    return urls


@bp.route("/sitemap-<any(systems,sectors):kind>-<int:n>.xml")
@_file
def sitemap_detail(kind, n):
    """One file of system or sector addresses; 404 unless the detail pages may be indexed."""
    settings = _sitemaps_on()
    if not seo.detail_indexed(settings) or n < 1:
        abort(404)
    db, size = db_name(), settings.seo.sitemap_max_urls
    build = (lambda: _system_chunk(db, n, size)) if kind == "systems" else (lambda: _sector_chunk(db, n, size))
    urls = seo.cached((kind, db, n, size), lambda: seo.sitemap_urls_xml(build()))
    return _text(urls, "application/xml")


@bp.route("/.well-known/security.txt")
@_file
def security_txt():
    text = settingsModule.get_settings().seo.security_txt.strip()
    if not text:
        abort(404)
    return _text(text + "\n")


@bp.route("/humans.txt")
@_file
def humans_txt():
    text = settingsModule.get_settings().seo.humans_txt.strip()
    if not text:
        abort(404)
    return _text(text + "\n")


@bp.route("/site.webmanifest")
@_file
def site_webmanifest():
    import json
    settings = settingsModule.get_settings()
    if not settings.icons.manifest_enabled:
        abort(404)
    return _text(json.dumps(seo.manifest(settings), indent=2) + "\n", "application/manifest+json")


@bp.route("/site-assets/<name>")
@_file
def site_asset(name):
    """An uploaded icon or picture."""
    found = assets.read(name)
    if found is None:
        abort(404)
    data, mime = found
    response = Response(data, mimetype=mime)
    response.headers["Cache-Control"] = "public, max-age=31536000, immutable"
    response.headers["X-Content-Type-Options"] = "nosniff"
    response.headers["ETag"] = f'"{name.split(".")[0]}"'
    if mime == "image/svg+xml":
        response.headers["Content-Security-Policy"] = "default-src 'none'; sandbox"
    return response


def card_url(endpoint, view_args):
    """The preview card address of the page `endpoint` shows, or `None` (cards off, no renderer, no such card)."""
    settings = settingsModule.get_settings()
    kind = {"web.system": "system", "web.sector": "sector"}.get(endpoint)
    if not settings.og.generated_cards or kind not in CARD_RENDERERS:
        return None
    object_id = view_args.get("system_id") or view_args.get("sector_id")
    return url_for("web.og_card", kind=kind, object_id=object_id) if object_id else None


def _card_folder():
    return settingsModule.get_settings().og.card_cache_dir or os.path.join(settingsModule.state_dir(), "og-cards")


def _prune(folder, limit_bytes):
    """Removes the oldest drawn cards until the folder is within `limit_bytes`."""
    try:
        entries = [(e.stat().st_mtime, e.stat().st_size, e.path) for e in os.scandir(folder) if e.name.endswith(".png")]
    except OSError:
        return
    total = sum(size for _m, size, _p in entries)
    for _mtime, size, path in sorted(entries):
        if total <= limit_bytes:
            break
        try:
            os.unlink(path)
            total -= size
        except OSError:
            pass


@bp.route("/og/<any(system,sector):kind>/<uid:object_id>.png")
@_file
def og_card(kind, object_id):
    """The drawn preview card of one object; the default picture when none can be drawn."""
    settings = settingsModule.get_settings()
    renderer = CARD_RENDERERS.get(kind)
    fallback = redirect(assets.og_image(settings)["path"], 302)
    if not settings.og.generated_cards or renderer is None:
        return fallback
    folder = _card_folder()
    key = hashlib.sha256(f"{db_name()}|{kind}|{object_id}".encode()).hexdigest()[:32]
    path = os.path.join(folder, f"{key}.png")
    data = None
    try:
        with open(path, "rb") as handle:
            data = handle.read()
    except OSError:
        pass
    if data is None:
        try:
            data = renderer(object_id)
        except Exception:  # noqa: BLE001 -- a card that can't be drawn is the default picture, never a 500
            data = None
        if not data:
            return fallback
        try:
            os.makedirs(folder, exist_ok=True)
            with open(path, "wb") as handle:
                handle.write(data)
            _prune(folder, int(settings.og.card_cache_mb * 1024 * 1024))
        except OSError:
            pass
    response = Response(data, mimetype="image/png")
    response.headers["Cache-Control"] = "public, max-age=86400"
    return response
