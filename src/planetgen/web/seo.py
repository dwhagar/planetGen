# planetgen/web/seo.py

"""
Web, Open Graph and search-engine settings at request time (ADM.44, `docs/design/interface/settings-seo-and-accounts.md` A8).

Everything here reads `get_settings()` per request, so a change saved on the settings page applies to the next page
view with no restart.

- `head(title, description)` is what `base.html`'s `<head>` is built from: the title, language, robots meta,
  canonical link, icons, theme colours, Open Graph and Twitter tags and the home page's JSON-LD.
- `classify(endpoint)` puts every page into one class of the **indexing matrix**: `top` (index), `list` (a list page:
  indexable when it has no query), `detail` (generated system, sector and phenomenon pages: `noindex,follow` unless
  `seo.detail_pages` or `seo.indexing` says index), `utility` (data and tool routes: `X-Robots-Tag: noindex`,
  disallowed in robots.txt) and `private` (admin, login, account: `noindex,nofollow`, disallowed). A page not listed
  below is `utility`, the safe class, and a test fails until it is classified on purpose.
- `robots_txt`, `sitemap_*`, `manifest` build the generated files `seo_pages.py` serves.

The absolute address of the site is `base_url` (never the request's Host header); while it is still the default
`http://localhost/` no canonical link, `og:url` or `og:image` is written, since those must be absolute.

This module imports Flask only for the request.
"""

import json
import re
import threading
import time
from urllib.parse import urlencode, urlsplit

from flask import current_app, request, url_for

from planetgen.util import settings as settingsModule
from planetgen.web import assets

TOP = frozenset({
    "web.index", "web.help_index", "web.classes", "web.class_type_page", "web.class_page", "web.habitability_levels", "web.property_help_page", "web.life_help_page", "web.life_help_topic", "web.life_help_class", "web.galaxy",
    "web.species", "web.species_page", "web.polities", "web.polity_page",
})
LISTS = frozenset({"web.sectors", "web.systems", "web.phenomena", "web.lists", "web.object_list"})
DETAIL = frozenset({
    "web.system", "web.sector", "web.sector_at", "web.phenomenon", "web.uncharted_phenomenon", "web.uncharted_star",
    "web.object_page",
})
OTHER_PUBLIC_UTILITY = frozenset({
    "web.search", "web.nav", "web.nav_chart", "web.nearby", "web.table_data",
    "web.galaxy_locate", "web.galaxy_status", "web.galaxy_system_nearest", "web.galaxy_regions", "web.galaxy_near", "web.galaxy_nebula_shape", "web.galaxy_nebula_surroundings", "web.galaxy_stage", "web.galaxy_territories",
    "web.galaxy_nodes", "web.galaxy_uncharted_scene", "web.galaxy_window", "web.sector_on_galaxy_map",
    "web.sector_scene", "web.system_on_galaxy_map", "web.system_scene", "web.old_page", "web.object_tree",
})
PRIVATE_PREFIXES = ("web.admin", "web.generate", "web.download_galaxy_settings")
PRIVATE = frozenset({"web.login", "web.logout", "web.account", "web.maintenance_page", "web.preferences"})
GENERATED_FILES = frozenset({
    "web.robots_txt", "web.sitemap_index", "web.sitemap_pages", "web.sitemap_detail", "web.security_txt",
    "web.humans_txt", "web.site_webmanifest", "web.site_asset", "web.og_card",
})
"""frozenset: The files `seo_pages.py` serves; not pages, so not classified."""


def classify(endpoint):
    """`"top"`, `"list"`, `"detail"`, `"utility"` or `"private"` for a `web.*` endpoint."""
    if endpoint in TOP:
        return "top"
    if endpoint in LISTS:
        return "list"
    if endpoint in DETAIL:
        return "detail"
    if endpoint in PRIVATE or endpoint.startswith(PRIVATE_PREFIXES):
        return "private"
    return "utility"


def is_classified(endpoint):
    """Whether `endpoint` was put in a class on purpose (the drift test)."""
    return (endpoint in TOP or endpoint in LISTS or endpoint in DETAIL or endpoint in OTHER_PUBLIC_UTILITY
            or endpoint in PRIVATE or endpoint in GENERATED_FILES or endpoint.startswith(PRIVATE_PREFIXES))


_TRACKING = re.compile(r"^(utm_.*|fbclid|gclid|msclkid|ref|(.*_)?(sort|order|dir|desc))$")


def canonical_query(args):
    """The query parameters that identify a page: the tracking and sort parameters are dropped."""
    return [(key, value) for key, value in args.items(multi=True) if not _TRACKING.match(key) and value != ""]


def base_url():
    """The configured address with its trailing slash, or `None` while it is still the unconfigured default."""
    value = settingsModule.get_settings().base_url
    host = (urlsplit(value).hostname or "").lower()
    if host in ("localhost", "127.0.0.1", "::1", ""):
        return None
    return value


def absolute(path):
    """`path` (starting with `/`) on the configured address, or `None` without one."""
    base = base_url()
    return base + path.lstrip("/") if base else None


def indexing_of(endpoint, settings=None):
    """`(index, follow)` for a page of this endpoint, before the query is looked at."""
    settings = settings or settingsModule.get_settings()
    klass = classify(endpoint)
    if klass == "private":
        return False, False
    if settings.seo.indexing == "off" or klass == "utility":
        return False, True
    if klass == "detail":
        return settings.seo.detail_pages == "index" or settings.seo.indexing == "all", True
    return True, True


def _locale(language):
    if "-" in language:
        lang, region = language.split("-", 1)
        return f"{lang.lower()}_{region.upper()}"
    return {"en": "en_US", "de": "de_DE", "fr": "fr_FR", "es": "es_ES", "it": "it_IT", "pt": "pt_PT",
            "nl": "nl_NL", "sv": "sv_SE", "pl": "pl_PL", "ja": "ja_JP"}.get(language.lower())


def _title(settings, title):
    if title is None:
        return settings.site_name
    return settings.web.title_template.format(title=title, site_name=settings.site_name)


def head(title, description=None):
    """
    What the page's `<head>` needs, as a dict: `title`, `lang`, `description`, `keywords`, `robots` (the meta
    content, or `None` for the default), `canonical`, `theme_light`, `theme_dark`, `icons` (`[{rel, href, type,
    sizes}]`), `manifest`, `og` and `twitter` (`[(name, content)]`) and `json_ld` (a script body or `None`).
    """
    settings = settingsModule.get_settings()
    endpoint = request.endpoint or ""
    index, follow = indexing_of(endpoint, settings)
    try:
        query = canonical_query(request.args)
    except (UnicodeDecodeError, ValueError):
        query = []                                      # the error page for an undecodable query string
    klass = classify(endpoint)
    if klass in ("list", "top") and query and index:
        index = False                                   # ?page=N, filters: not indexed, but followed
    robots = None if (index and follow) else f"{'index' if index else 'noindex'},{'follow' if follow else 'nofollow'}"
    canonical = None
    base = base_url()
    if base and klass != "private" and request.endpoint:
        canonical = base + request.path.lstrip("/") + (("?" + urlencode(query)) if query else "")
    text = description or settings.seo.description or settings.site.tagline or None
    og, twitter = [], []
    if settings.og.enabled and klass != "private":
        og_site = settings.og.site_name or settings.site_name
        og = [("og:title", title or settings.site_name), ("og:type", "website"), ("og:site_name", og_site)]
        if canonical:
            og.append(("og:url", canonical))
        if text:
            og.append(("og:description", text))
        locale = _locale(settings.site.language)
        if locale:
            og.append(("og:locale", locale))
        image = assets.og_image(settings)
        if image and base:
            og += [("og:image", base + image["path"].lstrip("/")), ("og:image:type", image["mime"])]
            if image["width"]:
                og += [("og:image:width", str(image["width"])), ("og:image:height", str(image["height"]))]
            if settings.og.image_alt:
                og.append(("og:image:alt", settings.og.image_alt))
        twitter = [("twitter:card", settings.og.twitter_card)]
        if settings.og.twitter_site:
            twitter.append(("twitter:site", settings.og.twitter_site))
        if image and base and settings.og.image_alt:
            twitter.append(("twitter:image:alt", settings.og.image_alt))
    json_ld = None
    if settings.seo.json_ld and endpoint == "web.index" and base:
        document = {"@context": "https://schema.org", "@type": "WebSite", "name": settings.site_name, "url": base}
        json_ld = json.dumps(document).replace("<", "\\u003c")
    return {
        "title": _title(settings, title), "lang": settings.site.language, "description": text,
        "keywords": ", ".join(settings.seo.keywords) or None, "robots": robots, "canonical": canonical,
        "theme_light": settings.web.theme_color_light, "theme_dark": settings.web.theme_color_dark,
        "icons": assets.icon_links(settings), "manifest": url_for("web.site_webmanifest")
        if settings.icons.manifest_enabled else None,
        "og": og, "twitter": twitter, "json_ld": json_ld,
    }


def response_headers(response):
    """`after_request` for the `web` blueprint: `X-Robots-Tag: noindex` on data and tool routes (they have no
    `<head>` to carry a meta tag), `noindex, nofollow` on private ones."""
    endpoint = request.endpoint or ""
    if endpoint.startswith("web.") and endpoint not in GENERATED_FILES:
        klass = classify(endpoint)
        if klass == "private":
            response.headers.setdefault("X-Robots-Tag", "noindex, nofollow")
        elif klass == "utility":
            response.headers.setdefault("X-Robots-Tag", "noindex")
    return response


# ---------------------------------------------------------------------
# robots.txt
# ---------------------------------------------------------------------

AI_CRAWLERS = ("GPTBot", "ClaudeBot", "CCBot", "Google-Extended", "Applebot-Extended")


def _disallow_patterns():
    """`Disallow` paths for every utility and private page, from the routes themselves (a `<converter>` is `*`)."""
    patterns = set()
    for rule in current_app.url_map.iter_rules():
        endpoint = rule.endpoint
        if not endpoint.startswith("web.") or "GET" not in rule.methods or endpoint in GENERATED_FILES:
            continue
        if classify(endpoint) not in ("utility", "private") or endpoint == "web.old_page":
            continue
        patterns.add(re.sub(r"<[^>]+>", "*", rule.rule))
    patterns.add("/api/")
    patterns.add("/*?*page=")
    kept = []
    for pattern in sorted(patterns):
        if not any(pattern != other and not other.endswith("*") and pattern.startswith(other.rstrip("/") + "/")
                   for other in patterns if "*" not in other):
            kept.append(pattern)
    return kept


def robots_txt():
    """The robots.txt body. `seo.robots_override` replaces it whole; `seo.robots_extra` is appended."""
    settings = settingsModule.get_settings()
    if settings.seo.robots_override.strip():
        return settings.seo.robots_override.rstrip("\n") + "\n"
    lines = ["User-agent: *"]
    if settings.seo.indexing == "off":
        lines.append("Disallow: /")
    else:
        lines += [f"Disallow: {pattern}" for pattern in _disallow_patterns()]
    if settings.seo.block_ai_crawlers:
        for agent in AI_CRAWLERS:
            lines += ["", f"User-agent: {agent}", "Disallow: /"]
    sitemap = absolute("/sitemap.xml")
    if settings.seo.sitemap_enabled and settings.seo.indexing != "off" and sitemap:
        lines += ["", f"Sitemap: {sitemap}"]
    if settings.seo.robots_extra.strip():
        lines += ["", settings.seo.robots_extra.rstrip("\n")]
    return "\n".join(lines) + "\n"


def robots_warning(text, settings=None):
    """A sentence when the robots.txt `text` disallows `/` while `seo.indexing` is not `off` (a mistake that hides the
    site), else `None`."""
    settings = settings or settingsModule.get_settings()
    if settings.seo.indexing == "off":
        return None
    current = None
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        field, _colon, value = line.partition(":")
        field, value = field.strip().lower(), value.strip()
        if field == "user-agent":
            current = value
        elif field == "disallow" and value == "/" and current == "*":
            return "disallows everything while search engines are allowed to index the site"
    return None


# ---------------------------------------------------------------------
# Sitemaps
# ---------------------------------------------------------------------

_cache = {}
_cache_lock = threading.Lock()


def cached(key, build):
    """`build()` reused for `seo.sitemap_cache_hours`."""
    settings = settingsModule.get_settings()
    ttl = settings.seo.sitemap_cache_hours * 3600
    with _cache_lock:
        found = _cache.get(key)
        if found and time.monotonic() - found[0] < ttl:
            return found[1]
    value = build()
    with _cache_lock:
        _cache[key] = (time.monotonic(), value)
    return value


def clear_cache():
    """Forgets every built sitemap (tests, and after a settings save)."""
    with _cache_lock:
        _cache.clear()


def sitemap_urls_xml(urls):
    """A `<urlset>` of absolute addresses (`loc` only: Google ignores changefreq and priority)."""
    from xml.sax.saxutils import escape
    body = "".join(f"<url><loc>{escape(url)}</loc></url>" for url in urls)
    return ('<?xml version="1.0" encoding="UTF-8"?>\n'
            f'<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">{body}</urlset>\n')


def sitemap_index_xml(urls):
    from xml.sax.saxutils import escape
    body = "".join(f"<sitemap><loc>{escape(url)}</loc></sitemap>" for url in urls)
    return ('<?xml version="1.0" encoding="UTF-8"?>\n'
            f'<sitemapindex xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">{body}</sitemapindex>\n')


def top_page_paths():
    """The paths of the pages every sitemap lists: the ones that need no id, the home page first."""
    paths = set()
    for rule in current_app.url_map.iter_rules():
        if rule.endpoint in TOP | LISTS and not rule.arguments and "GET" in rule.methods:
            paths.add(rule.rule)
    return sorted(paths, key=lambda path: (path != "/", path))


def detail_indexed(settings=None):
    settings = settings or settingsModule.get_settings()
    return settings.seo.indexing != "off" and (settings.seo.detail_pages == "index" or settings.seo.indexing == "all")


def manifest(settings=None):
    """`/site.webmanifest` as a dict."""
    settings = settings or settingsModule.get_settings()
    icons = []
    for link in assets.icon_links(settings):
        if link["rel"] in ("icon", "apple-touch-icon") and link.get("sizes") not in (None, "any"):
            icons.append({"src": link["href"], "sizes": link["sizes"], "type": link["type"]})
    return {"name": settings.site_name, "short_name": settings.site_name[:12], "start_url": "/",
            "display": "standalone", "background_color": settings.web.theme_color_dark,
            "theme_color": settings.web.theme_color_light, "icons": icons, "lang": settings.site.language}
