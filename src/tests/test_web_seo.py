# tests/test_web_seo.py

"""
Web, Open Graph and search-engine settings (ADM.44): the page `<head>` (`planetgen.web.seo.head`), the indexing
matrix, robots.txt, the sitemaps, security.txt, humans.txt, the manifest, uploaded icons and pictures
(`planetgen.web.assets`, `/api/admin/images`, `/admin/settings/images`) and the preview card route, through the real
app and control database over a temporary `config.json`.

Run with: pytest tests/test_web_seo.py
"""

import base64
import io
import json
import re
import xml.etree.ElementTree as ET

import pytest
from PIL import Image

from planetgen.util import settings
from planetgen.web import assets, csrf, seo, seo_pages
from planetgen.web.app import create_app
from planetgen.api.config import Config

from tests.test_web_settings import START, _audit_rows, config, site  # noqa: F401 -- fixtures and helpers

SITEMAP = "{http://www.sitemaps.org/schemas/sitemap/0.9}"


def _configure(config, **changes):
    document = dict(START)
    document.update(changes)
    config.write_text(json.dumps(document, indent=2) + "\n")
    settings.reset_cache()
    seo.clear_cache()
    seo_pages.clear_caches()


def _png(width=1200, height=630, color=(10, 20, 30)):
    buffer = io.BytesIO()
    Image.new("RGB", (width, height), color).save(buffer, "PNG")
    return buffer.getvalue()


def _jpeg(width=1200, height=630):
    buffer = io.BytesIO()
    Image.new("RGB", (width, height), (200, 10, 10)).save(buffer, "JPEG")
    return buffer.getvalue()


def _other_format(name):
    """A picture in a format Pillow reads but the upload does not take (SEC.35)."""
    buffer = io.BytesIO()
    Image.new("RGB", (8, 8), (10, 20, 30)).save(buffer, name)
    return buffer.getvalue()


SVG = b'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 10 10"><circle cx="5" cy="5" r="4" fill="#fc0"/></svg>'


def _upload(site, kind, data):
    return site.post(f"/api/admin/images/{kind}", json={"data": base64.b64encode(data).decode()})


@pytest.fixture(autouse=True)
def _assets_folder(tmp_path, monkeypatch):
    monkeypatch.setenv("PLANETGEN_ASSETS_DIR", str(tmp_path / "assets"))
    settings.reset_cache()
    yield
    settings.reset_cache()
    seo.clear_cache()
    seo_pages.clear_caches()


def _tags(page, name, attr="name"):
    return re.findall(rf'<meta {attr}="{re.escape(name)}" content="([^"]*)"', page)


# --- The indexing matrix ----------------------------------------------------------------

def test_every_web_page_is_classified_on_purpose():
    app = create_app(Config)
    endpoints = {r.endpoint for r in app.url_map.iter_rules() if r.endpoint.startswith("web.") and "GET" in r.methods}
    unclassified = sorted(e for e in endpoints if not seo.is_classified(e))
    assert unclassified == [], f"add these to the indexing matrix in web/seo.py: {unclassified}"


@pytest.mark.parametrize("endpoint, klass", [
    ("web.index", "top"), ("web.classes", "top"), ("web.sectors", "list"), ("web.system", "detail"),
    ("web.sector", "detail"), ("web.search", "utility"), ("web.galaxy_nodes", "utility"), ("web.login", "private"),
    ("web.admin_settings", "private"), ("web.generate_job_log", "private"), ("web.not_a_page_yet", "utility"),
    ("web.preferences", "private"),
])
def test_the_classes(endpoint, klass):
    assert seo.classify(endpoint) == klass


def test_what_may_be_indexed_follows_the_settings(config):
    _configure(config)
    assert seo.indexing_of("web.index") == (True, True)
    assert seo.indexing_of("web.system") == (False, True)          # detail pages stay out by default
    assert seo.indexing_of("web.search") == (False, True)
    assert seo.indexing_of("web.admin") == (False, False)
    _configure(config, seo={"detail_pages": "index"})
    assert seo.indexing_of("web.system") == (True, True)
    _configure(config, seo={"indexing": "all"})
    assert seo.indexing_of("web.sector") == (True, True)
    _configure(config, seo={"indexing": "off"})
    assert seo.indexing_of("web.index") == (False, True)


# --- The head ---------------------------------------------------------------------------

def test_without_a_base_url_there_is_no_absolute_link(site, config):
    _configure(config)
    page = site.get("/classes").get_data(as_text=True)
    assert 'rel="canonical"' not in page and "og:url" not in page and "og:image" not in page
    assert "<title>Classes - Mine</title>" in page or "- Mine</title>" in page
    assert '<html lang="en">' in page


def test_the_head_of_a_public_page(site, config):
    _configure(config, base_url="https://galaxy.example.org", site={"tagline": "Every star", "language": "de"},
               web={"title_template": "{site_name} | {title}"}, seo={"keywords": ["stars", "planets"]},
               og={"image_alt": "A galaxy", "twitter_site": "@pg"})
    page = site.get("/classes?utm_source=x&sort=name&q=g").get_data(as_text=True)
    assert '<html lang="de">' in page
    assert re.search(r"<title>Mine \| [^<]+</title>", page)
    assert '<link rel="canonical" href="https://galaxy.example.org/classes?q=g">' in page   # tracking and sort dropped
    assert '<meta name="robots"' in page and "noindex,follow" in page                       # a query: not indexed
    assert _tags(page, "description") == ["Every star"] or _tags(page, "description")
    assert _tags(page, "keywords") == ["stars, planets"]
    assert _tags(page, "og:site_name", "property") == ["Mine"]
    assert _tags(page, "og:locale", "property") == ["de_DE"]
    assert _tags(page, "og:image", "property") == ["https://galaxy.example.org/static/og-default.png"]
    assert _tags(page, "og:image:width", "property") == ["1200"]
    assert _tags(page, "og:image:alt", "property") == ["A galaxy"]
    assert _tags(page, "twitter:card") == ["summary_large_image"] and _tags(page, "twitter:site") == ["@pg"]
    plain = site.get("/classes").get_data(as_text=True)
    assert '<meta name="robots"' not in plain and 'href="https://galaxy.example.org/classes"' in plain


def test_open_graph_can_be_switched_off(site, config):
    _configure(config, base_url="https://galaxy.example.org/", og={"enabled": False})
    page = site.get("/classes").get_data(as_text=True)
    assert "og:title" not in page and "twitter:card" not in page and 'rel="canonical"' in page


def test_private_pages_say_noindex_nofollow_and_carry_no_preview(site, config):
    _configure(config, base_url="https://galaxy.example.org/")
    for client, path in ((site.application.test_client(), "/login"), (site, "/admin/settings"),
                         (site.application.test_client(), "/preferences")):
        response = client.get(path)
        page = response.get_data(as_text=True)
        assert "noindex,nofollow" in page and "og:title" not in page and 'rel="canonical"' not in page
        assert response.headers["X-Robots-Tag"] == "noindex, nofollow"


def test_data_routes_carry_the_header_instead_of_a_meta_tag(site, config):
    _configure(config)
    assert site.get("/search").headers["X-Robots-Tag"] == "noindex"
    assert "X-Robots-Tag" not in site.get("/classes").headers
    assert "X-Robots-Tag" not in site.get("/robots.txt").headers


def test_indexing_off_says_noindex_everywhere(site, config):
    _configure(config, seo={"indexing": "off"})
    assert "noindex,follow" in site.get("/classes").get_data(as_text=True)


def test_the_theme_colours_and_json_ld(site, config):
    _configure(config, base_url="https://galaxy.example.org/", web={"theme_color_light": "#ffffff"})
    page = site.get("/").get_data(as_text=True)
    assert 'content="#ffffff" media="(prefers-color-scheme: light)"' in page
    block = re.search(r'<script type="application/ld\+json">(.*?)</script>', page, re.S)
    assert block and json.loads(block.group(1)) == {"@context": "https://schema.org", "@type": "WebSite",
                                                    "name": "Mine", "url": "https://galaxy.example.org/"}
    assert "ld+json" not in site.get("/classes").get_data(as_text=True)
    _configure(config, base_url="https://galaxy.example.org/", seo={"json_ld": False})
    assert "ld+json" not in site.get("/").get_data(as_text=True)


def test_the_tagline_sits_under_the_home_heading(site, config):
    _configure(config, site={"tagline": "Every star we know"})
    assert '<p class="tagline">Every star we know</p>' in site.get("/").get_data(as_text=True)


def test_the_privacy_note_is_on_the_account_page(site, config):
    _configure(config, seo={"privacy_note": "We keep a login cookie."})
    assert "We keep a login cookie." in site.get("/account").get_data(as_text=True)


# --- robots.txt -------------------------------------------------------------------------

def test_robots_txt_disallows_the_tools_and_names_the_sitemap(site, config):
    _configure(config, base_url="https://galaxy.example.org/")
    response = site.get("/robots.txt")
    text = response.get_data(as_text=True)
    assert response.mimetype == "text/plain"
    lines = text.splitlines()
    assert lines[0] == "User-agent: *"
    for path in ("/admin", "/search", "/nav", "/api/", "/galaxy/nodes", "/login", "/*?*page=", "/system/*/scene"):
        assert f"Disallow: {path}" in lines, path
    assert "Disallow: /system/" not in lines and "Disallow: /system" not in lines and "Disallow: /classes" not in text
    assert lines[-1] == "Sitemap: https://galaxy.example.org/sitemap.xml"


def test_robots_txt_when_indexing_is_off_and_with_ai_crawlers_blocked(site, config):
    _configure(config, seo={"indexing": "off", "block_ai_crawlers": True})
    text = site.get("/robots.txt").get_data(as_text=True)
    assert "Disallow: /\n" in text and "Sitemap:" not in text
    for agent in ("GPTBot", "ClaudeBot", "CCBot", "Google-Extended", "Applebot-Extended"):
        assert f"User-agent: {agent}\nDisallow: /" in text


def test_robots_extra_is_appended_and_override_replaces(site, config):
    _configure(config, seo={"robots_extra": "User-agent: Foo\nDisallow: /bar"})
    text = site.get("/robots.txt").get_data(as_text=True)
    assert text.startswith("User-agent: *") and text.rstrip().endswith("Disallow: /bar")
    _configure(config, seo={"robots_override": "User-agent: *\nAllow: /"})
    assert site.get("/robots.txt").get_data(as_text=True) == "User-agent: *\nAllow: /\n"


def test_a_robots_file_with_a_bad_line_is_refused_at_the_option(config):
    assert settings.robots_problems("User-agent: *\nDisallow: /x") == []
    assert settings.robots_problems("hello")
    document = dict(START, seo={"robots_extra": "not a directive"})
    config.write_text(json.dumps(document))
    settings.reset_cache()
    errors = settings.field_errors(document)[1]
    assert "seo.robots_extra" in [path for path, _message in errors]


def test_an_override_that_hides_the_site_is_a_warning(config):
    _configure(config)
    assert "disallows everything" in seo.robots_warning("User-agent: *\nDisallow: /\n")
    assert seo.robots_warning("User-agent: *\nDisallow: /admin\n") is None
    _configure(config, seo={"indexing": "off"})
    assert seo.robots_warning("User-agent: *\nDisallow: /\n") is None


# --- Sitemaps ---------------------------------------------------------------------------

def _urls(response):
    return [e.text for e in ET.fromstring(response.get_data()).iter(SITEMAP + "loc")]


def test_there_is_no_sitemap_without_a_base_url_or_when_off(site, config):
    _configure(config)
    assert site.get("/sitemap.xml").status_code == 404
    _configure(config, base_url="https://galaxy.example.org/", seo={"sitemap_enabled": False})
    assert site.get("/sitemap.xml").status_code == 404 and site.get("/sitemap-pages.xml").status_code == 404
    _configure(config, base_url="https://galaxy.example.org/", seo={"indexing": "off"})
    assert site.get("/sitemap.xml").status_code == 404


def test_the_index_lists_the_pages_file_and_the_pages_file_the_top_pages(site, config):
    _configure(config, base_url="https://galaxy.example.org/")
    index = site.get("/sitemap.xml")
    assert index.mimetype == "application/xml"
    assert _urls(index) == ["https://galaxy.example.org/sitemap-pages.xml"]
    pages = _urls(site.get("/sitemap-pages.xml"))
    assert pages[0] == "https://galaxy.example.org/"
    for path in ("classes", "galaxy", "sectors", "systems", "phenomena", "species", "polities"):
        assert f"https://galaxy.example.org/{path}" in pages, path
    assert not any("/admin" in u or "/search" in u or "/nav" in u for u in pages)
    assert site.get("/sitemap-systems-1.xml").status_code == 404       # detail pages are not indexed


def _fake_api(monkeypatch, systems, sectors):
    calls = []

    def get_systems(db, limit=None, after=None, **kwargs):
        calls.append(("systems", limit, after))
        start = 0 if after is None else next(i for i, s in enumerate(systems) if s["id"] == after) + 1
        return {"items": systems[start:start + (limit or 50)], "total": len(systems)}

    def get_sectors(db, limit=None, offset=0, **kwargs):
        calls.append(("sectors", limit, offset))
        return {"items": sectors[offset:offset + (limit or 50)], "total": len(sectors)}

    monkeypatch.setattr(seo_pages.apiclient, "get_systems", get_systems)
    monkeypatch.setattr(seo_pages.apiclient, "get_sectors", get_sectors)
    return calls


def test_detail_sitemaps_split_at_the_limit_and_are_cached(site, config, monkeypatch):
    systems = [{"id": f"FE8100000A2B-000000{i}-000"} for i in range(5)]
    sectors = [{"id": f"{i:X}"} for i in range(1, 4)]
    calls = _fake_api(monkeypatch, systems, sectors)
    _configure(config, base_url="https://galaxy.example.org/", seo={"detail_pages": "index", "sitemap_max_urls": 2})
    index = _urls(site.get("/sitemap.xml"))
    assert index == ["https://galaxy.example.org/sitemap-pages.xml"] + [
        f"https://galaxy.example.org/sitemap-systems-{n}.xml" for n in (1, 2, 3)] + [
        f"https://galaxy.example.org/sitemap-sectors-{n}.xml" for n in (1, 2)]
    got = [u for n in (1, 2, 3) for u in _urls(site.get(f"/sitemap-systems-{n}.xml"))]
    assert got == [f"https://galaxy.example.org/system/{s['id']}" for s in systems]
    assert [len(_urls(site.get(f"/sitemap-sectors-{n}.xml"))) for n in (1, 2)] == [2, 1]
    assert _urls(site.get("/sitemap-systems-4.xml")) == []
    before = len(calls)
    site.get("/sitemap-systems-2.xml")
    assert len(calls) == before                                          # served from the cache


def test_a_system_file_asked_first_walks_from_the_start(site, config, monkeypatch):
    systems = [{"id": f"FE8100000A2B-000000{i}-000"} for i in range(5)]
    _fake_api(monkeypatch, systems, [])
    _configure(config, base_url="https://galaxy.example.org/", seo={"indexing": "all", "sitemap_max_urls": 2})
    assert _urls(site.get("/sitemap-systems-3.xml")) == [f"https://galaxy.example.org/system/{systems[4]['id']}"]
    assert len(_urls(site.get("/sitemap-systems-2.xml"))) == 2


# --- The small files --------------------------------------------------------------------

def test_security_and_humans_txt_are_404_until_there_is_text(site, config):
    _configure(config)
    assert site.get("/.well-known/security.txt").status_code == 404 and site.get("/humans.txt").status_code == 404
    _configure(config, seo={"security_txt": "Contact: mailto:a@example.org\nExpires: 2027-01-01T00:00:00Z",
                            "humans_txt": "Built by people."})
    security = site.get("/.well-known/security.txt")
    assert security.mimetype == "text/plain" and security.get_data(as_text=True).startswith("Contact: mailto:a@")
    assert site.get("/humans.txt").get_data(as_text=True) == "Built by people.\n"


def test_the_manifest_names_the_site_and_its_icons(site, config):
    _configure(config)
    assert site.get("/site.webmanifest").status_code == 404
    assert 'rel="manifest"' not in site.get("/classes").get_data(as_text=True)
    _configure(config, icons={"manifest_enabled": True}, site={"language": "en-GB"})
    response = site.get("/site.webmanifest")
    body = json.loads(response.get_data())
    assert response.mimetype == "application/manifest+json"
    assert body["name"] == "Mine" and body["theme_color"] == "#f6f7fb" and body["lang"] == "en-GB"
    assert 'rel="manifest" href="/site.webmanifest"' in site.get("/classes").get_data(as_text=True)


# --- Uploaded pictures ------------------------------------------------------------------

def test_the_images_api_lists_every_kind(site, config):
    _configure(config)
    body = site.get("/api/admin/images").get_json()
    assert [i["kind"] for i in body["images"]] == ["favicon_svg", "favicon_png", "apple_touch_icon", "default_image"]
    assert all(i["name"] is None and i["url"] is None for i in body["images"])
    assert body["built_in_preview"] == "/static/og-default.png"
    assert site.application.test_client().get("/api/admin/images").status_code == 401


def test_an_uploaded_preview_picture_is_stored_named_by_content_and_used(site, config):
    _configure(config, base_url="https://galaxy.example.org/")
    data = _jpeg()
    response = _upload(site, "default_image", data)
    assert response.status_code == 200
    image = response.get_json()["image"]
    assert re.fullmatch(r"[0-9a-f]{16}\.jpg", image["name"]) and image["url"] == f"/site-assets/{image['name']}"
    assert (image["width"], image["height"], image["mime"]) == (1200, 630, "image/jpeg")
    served = site.get(image["url"])
    assert served.get_data() == data and served.mimetype == "image/jpeg"
    assert "immutable" in served.headers["Cache-Control"] and served.headers["X-Content-Type-Options"] == "nosniff"
    assert json.loads(config.read_text())["og"]["default_image"] == image["name"]
    page = site.get("/classes").get_data(as_text=True)
    assert _tags(page, "og:image", "property") == [f"https://galaxy.example.org/site-assets/{image['name']}"]
    assert ("config.update", "config:og.default_image") in [(a, t) for a, t, _d in _audit_rows(site)]
    again = _upload(site, "default_image", data).get_json()["image"]
    assert again["name"] == image["name"]                                # same picture, same file


def test_replacing_a_picture_removes_the_old_file(site, config):
    _configure(config)
    first = _upload(site, "default_image", _png(1200, 630)).get_json()["image"]["name"]
    second = _upload(site, "default_image", _png(1200, 631)).get_json()["image"]["name"]
    assert first != second and site.get(f"/site-assets/{first}").status_code == 404
    assert site.get(f"/site-assets/{second}").status_code == 200


def test_removing_a_picture_goes_back_to_the_built_in_one(site, config):
    _configure(config)
    name = _upload(site, "default_image", _png()).get_json()["image"]["name"]
    removed = site.delete("/api/admin/images/default_image")
    assert removed.status_code == 200 and removed.get_json()["image"]["name"] is None
    assert "default_image" not in json.loads(config.read_text()).get("og", {})
    assert site.get(f"/site-assets/{name}").status_code == 404


def test_icons_appear_in_the_head(site, config):
    _configure(config)
    svg = _upload(site, "favicon_svg", SVG).get_json()["image"]
    png = _upload(site, "favicon_png", _png(64, 64)).get_json()["image"]
    touch = _upload(site, "apple_touch_icon", _png(180, 180)).get_json()["image"]
    assert site.get(svg["url"]).headers["Content-Security-Policy"] == "default-src 'none'; sandbox"
    page = site.get("/classes").get_data(as_text=True)
    assert f'<link rel="icon" type="image/svg+xml" href="{svg["url"]}" sizes="any">' in page
    assert f'<link rel="icon" type="image/png" href="{png["url"]}" sizes="64x64">' in page
    assert f'<link rel="apple-touch-icon" type="image/png" href="{touch["url"]}" sizes="180x180">' in page
    assert "/static/favicon.svg" not in page


@pytest.mark.parametrize("kind, data, why", [
    ("default_image", b"", "empty"),
    ("default_image", b"not a picture", "not a readable picture"),
    pytest.param("default_image", _other_format("TIFF"), "not a readable picture", id="tiff-is-not-read"),
    pytest.param("default_image", _other_format("GIF"), "not a readable picture", id="gif-is-not-read"),
    pytest.param("default_image", _other_format("PPM"), "not a readable picture", id="ppm-is-not-read"),
    ("default_image", _png(300, 100), "at least 600 by 315"),
    ("default_image", SVG, "PNG or JPG"),
    ("apple_touch_icon", _png(180, 120), "must be square"),
    ("favicon_png", _png(16, 16), "at least 32 by 32"),
    ("favicon_svg", _png(), "must be an SVG"),
    ("favicon_svg", b'<svg xmlns="http://www.w3.org/2000/svg"><script>alert(1)</script></svg>', "<script>"),
    ("favicon_svg", b'<svg xmlns="http://www.w3.org/2000/svg" onload="x()"></svg>', "event handler"),
    ("favicon_svg", b'<svg xmlns="http://www.w3.org/2000/svg"><image href="http://evil/x.png"/></svg>', "outside"),
    ("favicon_svg", b'<!DOCTYPE svg [<!ENTITY a "b">]><svg xmlns="http://www.w3.org/2000/svg"/>', "DOCTYPE"),
    ("favicon_svg", b"<html></html>", "not an SVG"),
    ("favicon_svg", b'<svg xmlns="http://www.w3.org/2000/svg"' + b" " * 1_100_000 + b"/>", "most allowed"),
])
def test_a_picture_that_cannot_be_used_is_refused_with_the_reason(site, config, kind, data, why):
    _configure(config)
    response = _upload(site, kind, data)
    assert response.status_code == 422, response.get_data(as_text=True)
    assert why in response.get_json()["error"] or why in json.dumps(response.get_json())
    assert not (assets.assets_dir() and __import__("os").path.exists(assets.assets_dir())
                and __import__("os").listdir(assets.assets_dir()))


def test_bad_requests_to_the_images_api(site, config):
    _configure(config)
    assert site.post("/api/admin/images/nothing", json={"data": "AAAA"}).status_code == 404
    assert site.post("/api/admin/images/default_image", json={"data": "***"}).status_code == 422
    assert site.post("/api/admin/images/default_image", json={"data": "AAAA", "extra": 1}).status_code == 400
    assert site.get("/site-assets/../../etc/passwd").status_code == 404
    assert site.get("/site-assets/0123456789abcdef.png").status_code == 404


def test_the_images_page_uploads_and_removes(site, config):
    _configure(config)
    page = site.get("/admin/settings/images").get_data(as_text=True)
    assert "Link preview picture" in page and "the built-in picture" in page
    resp = site.post("/admin/settings/images", data={
        csrf.FIELD_NAME: site.csrf, "kind": "default_image", "file": (io.BytesIO(_png()), "card.png")},
        content_type="multipart/form-data")
    assert resp.status_code == 303
    page = site.get("/admin/settings/images").get_data(as_text=True)
    assert "Link preview picture: saved." in page and "Use the built-in one" in page
    resp = site.post("/admin/settings/images", data={csrf.FIELD_NAME: site.csrf, "kind": "default_image", "remove": "1"})
    assert resp.status_code == 303
    assert "back to the built-in one" in site.get("/admin/settings/images").get_data(as_text=True)


def test_the_images_page_shows_why_an_upload_was_refused(site, config):
    _configure(config)
    site.post("/admin/settings/images", data={
        csrf.FIELD_NAME: site.csrf, "kind": "apple_touch_icon", "file": (io.BytesIO(_png(180, 150)), "a.png")},
        content_type="multipart/form-data")
    assert "must be square" in site.get("/admin/settings/images").get_data(as_text=True)


# --- Preview cards ----------------------------------------------------------------------

def test_a_card_is_the_default_picture_until_a_renderer_is_installed(site, config, monkeypatch):
    _configure(config, base_url="https://galaxy.example.org/", og={"generated_cards": True})
    monkeypatch.setattr(seo_pages, "CARD_RENDERERS", {})
    response = site.get("/og/system/FE8100000A2B-0000001-000.png")
    assert response.status_code == 302 and response.headers["Location"] == "/static/og-default.png"


def test_a_registered_renderer_draws_once_and_is_cached(site, config, monkeypatch, tmp_path):
    _configure(config, base_url="https://galaxy.example.org/", og={"generated_cards": True, "card_cache_dir": str(tmp_path / "cards")})
    drawn = []
    monkeypatch.setattr(seo_pages, "CARD_RENDERERS", {"system": lambda object_id: drawn.append(object_id) or _png(1200, 630)})
    path = "/og/system/FE8100000A2B-0000001-000.png"
    first = site.get(path)
    assert first.status_code == 200 and first.mimetype == "image/png" and first.get_data() == _png(1200, 630)
    assert site.get(path).status_code == 200 and len(drawn) == 1
    with site.application.test_request_context():
        assert seo_pages.card_url("web.system", {"system_id": "FE8100000A2B-0000001-000"}) == path
        assert seo_pages.card_url("web.sector", {"sector_id": "1"}) is None       # no sector renderer


def test_a_renderer_that_fails_falls_back_to_the_default_picture(site, config, monkeypatch, tmp_path):
    _configure(config, og={"generated_cards": True, "card_cache_dir": str(tmp_path / "cards")})

    def broken(object_id):
        raise RuntimeError("no font")

    monkeypatch.setattr(seo_pages, "CARD_RENDERERS", {"system": broken})
    assert site.get("/og/system/FE8100000A2B-0000001-000.png").status_code == 302


def test_old_cards_are_pruned_to_the_size_limit(tmp_path):
    for i in range(4):
        path = tmp_path / f"{i}.png"
        path.write_bytes(b"x" * 1000)
        import os
        os.utime(path, (1000 + i, 1000 + i))
    seo_pages._prune(str(tmp_path), 2500)
    assert sorted(p.name for p in tmp_path.iterdir()) == ["2.png", "3.png"]


# --- The options ------------------------------------------------------------------------

def test_the_new_options_are_validated(config):
    def errors(**section):
        return settings.field_errors(dict(START, **section))[1]

    assert errors(web={"theme_color_light": "red"}) and errors(web={"title_template": "no title here"})
    assert errors(web={"title_template": "{title} {nope}"})
    assert errors(site={"language": "english please"}) and errors(og={"twitter_site": "pg"})
    assert errors(icons={"favicon_svg": "../etc/passwd"}) and errors(og={"default_image": "x.png"})
    assert errors(seo={"sitemap_max_urls": 50001}) and errors(seo={"sitemap_cache_hours": 0})
    assert errors(site={"tagline": "x" * 121}) and errors(base_url="ftp://x/") and errors(base_url="http://u:p@x/")
    assert not errors(base_url="https://galaxy.example.org", og={"twitter_site": "@pg_1"}, icons={"favicon_svg": "0123456789abcdef.svg"})
    assert settings.field_errors(dict(START, base_url="https://galaxy.example.org"))[0].base_url == "https://galaxy.example.org/"


def test_the_robots_boxes_are_multiline_on_the_settings_page(site, config):
    page = site.get("/admin/settings").get_data(as_text=True)
    assert '<textarea' in page and "Images (icons and link preview picture)" in page


def test_the_error_page_for_an_undecodable_query_string_still_has_a_head(site, config):
    _configure(config, base_url="https://galaxy.example.org/")
    response = site.get("/classes", environ_overrides={"QUERY_STRING": "q=\xff\xfe"})
    assert response.status_code == 400 and "<title>" in response.get_data(as_text=True)
