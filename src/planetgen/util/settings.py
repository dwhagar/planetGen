# planetgen/util/settings.py

"""
The settings model (ADM.42)
===========================

One Pydantic model, `Settings`, describes every `config.json` option: its
type, default, help text, unit, whether it is a secret, whether changing it
needs a restart, which environment variable overrides it and whether the
web may edit it. The loader (`get_settings`), the generated tables in
`docs/config.md`, `config.json.example` and `config.schema.json`
(`planetgen.cli.config`) all read the model, and tests fail when they
drift or when code reads a `PLANETGEN_*` variable the model doesn't name.

One file, one model (Boss, 2026-10-16: "1 file for all program configuration
settings"): `config.json` in the install folder, which the installers give to the
web server's user, holds every option, and the configuration page and the command
line (`planetgen.util.settings_writer`) edit it. Load order, lowest to highest: the
model's defaults, `config.json`, the `PLANETGEN_*` environment variables, then an
explicit argument at the call site.

Field metadata sits in `json_schema_extra` under `x-` keys, written with
`opt` below:

- `x-env`: the environment variable that overrides it
- `x-env-empty`: an environment variable that is set but empty still
  overrides (a blank MySQL password), instead of counting as unset
- `x-restart`: a change only applies after a restart
- `x-secret`: never shown in a page or log
- `x-editable`: `false` means the web may not change it (the database, secrets
  that open the host, `jobs.python`); only the command line, run on the server, can
- `x-min_role`: `admin` (the default) or `owner`: who may save it from the web
- `x-multiline`: a text option the form edits in a text area (robots.txt text, security.txt)
- `x-unit`, `x-category`: the form's suffix and section

The log options (`debug`, `log_file`, `log_dir`, `log_rotation`) are read
by `planetgen.util.logpaths`, which needs only the standard library (the
installers run it as root before any virtual environment exists); their
defaults come from there, and a test checks the two agree.

Annotations here are `Optional[...]`/`List[...]`, as in most of the package.
"""

import copy
import hashlib
import json
import os
import re
import threading
from typing import Any, List, Literal, get_args, get_origin
from urllib.parse import urlsplit

from pydantic import BaseModel, ConfigDict, Field, ValidationError, field_validator

from planetgen.util import logpaths

CONFIG_VERSION = 1
"""int: The `config_version` this code writes and reads. A file without one
is version 1; a newer one is refused."""

ENV_PREFIX = "PLANETGEN_"

_FALSE_STRINGS = ("", "0", "false", "no", "off")


def opt(default, description, env=None, restart=False, secret=False, editable=True, role="admin",
        unit=None, category="general", env_empty=False, multiline=False, **constraints):
    """A model field with its `x-` metadata (see the module docstring)."""
    extra = {"x-restart": restart, "x-secret": secret, "x-editable": editable, "x-min_role": role,
             "x-category": category, "x-multiline": multiline}
    if env:
        extra["x-env"] = env
    if env_empty:
        extra["x-env-empty"] = True
    if unit:
        extra["x-unit"] = unit
    return Field(default, description=description, json_schema_extra=extra, **constraints)


def _loose_bool(value):
    """`false`, `0`, `no`, `off` and empty text mean off, any other text on
    (what `PLANETGEN_DEBUG` always did)."""
    if isinstance(value, str):
        return value.strip().lower() not in _FALSE_STRINGS
    return value


class _Section(BaseModel):
    model_config = ConfigDict(extra="ignore", validate_default=True)


class MySQL(_Section):
    host: str = opt("127.0.0.1", "The MySQL or MariaDB server's host name or address.", env="PLANETGEN_MYSQL_HOST",
                    restart=True, editable=False, category="database", env_empty=True)
    port: int = opt(3306, "The server's port.", env="PLANETGEN_MYSQL_PORT", restart=True, editable=False,
                    category="database", ge=1, le=65535, env_empty=True)
    user: str = opt("planetgen", "The account every entry point connects with; give it whatever grants the most "
                    "demanding caller needs.", env="PLANETGEN_MYSQL_USER", restart=True, editable=False,
                    category="database", env_empty=True)
    password: str = opt("", "That account's password.", env="PLANETGEN_MYSQL_PASSWORD", restart=True, secret=True,
                        editable=False, category="database", env_empty=True)
    database: str = opt("planetgen", "The galaxy database the pages show and the generators fill.",
                        env="PLANETGEN_MYSQL_DATABASE", restart=True, editable=False, category="database",
                        env_empty=True)
    database_prefix: str = opt("planetgen", "The schema-name prefix the database picker lists, for a deployment "
                               "with several galaxies on one server.", env="PLANETGEN_MYSQL_DATABASE_PREFIX",
                               restart=True, editable=False, category="database")
    statement_timeout_seconds: float = opt(
        10, "The longest any statement on the web interface's and API's read-only connections may run (0 turns "
        "the limit off). A query it stops gives a 504 page, or a QUERY_TIMEOUT error from the API.",
        unit="seconds", category="database", ge=0)


class Redis(_Section):
    url: str = opt("redis://127.0.0.1:6379/0", "The Redis server the work queue, the rate limits and the login "
                   "lockouts use.", env="PLANETGEN_REDIS_URL", restart=True, editable=False, category="queue")


class RatePages(_Section):
    search: str = opt("30 per minute", "Limit per client address on /search and /galaxy/locate. Empty turns it off.",
                      restart=True, category="rate limits")
    galaxy: str = opt("60 per minute", "Limit on /galaxy, the Galaxy Map page. Empty turns it off.", restart=True,
                      category="rate limits")
    galaxy_tiles: str = opt("600 per minute", "Limit on /galaxy/tiles, /galaxy/stage, /galaxy/window and /galaxy/territories, "
                            "fetched as the map's camera moves. Empty turns it off.", restart=True,
                            category="rate limits")
    health: str = opt("60 per minute", "Limit on /api/health. Empty turns it off.", restart=True,
                      category="rate limits")
    other: str = opt("300 per minute", "Limit on every other page, counted together. Empty turns it off.",
                     restart=True, category="rate limits")


class RateLimit(_Section):
    default: str = opt("200 per day;50 per hour", "Flask-Limiter's default limit for the API.",
                       env="PLANETGEN_RATELIMIT_DEFAULT", restart=True, category="rate limits")
    storage_uri: str = opt("", "Where limits count; empty means the Redis server in redis.url, memory:// counts "
                           "per process (right with exactly one worker).", env="PLANETGEN_RATELIMIT_STORAGE_URI",
                           restart=True, category="rate limits")
    pages: RatePages = RatePages()


class ProxyFix(_Section):
    @field_validator("x_for", "x_proto", "x_host", mode="before")
    @classmethod
    def _not_a_bool(cls, value):
        if isinstance(value, bool):
            raise ValueError("must be a whole number of proxies (0 = off), not true or false")
        return value

    x_for: int = opt(0, "Reverse proxies to trust for X-Forwarded-For, the client address (0 = off).",
                     env="PLANETGEN_PROXY_FIX_X_FOR", restart=True, editable=False, category="proxy", ge=0)
    x_proto: int = opt(0, "Reverse proxies to trust for X-Forwarded-Proto, the scheme (0 = off).",
                       env="PLANETGEN_PROXY_FIX_X_PROTO", restart=True, editable=False, category="proxy", ge=0)
    x_host: int = opt(0, "Reverse proxies to trust for X-Forwarded-Host (0 = off).",
                      env="PLANETGEN_PROXY_FIX_X_HOST", restart=True, editable=False, category="proxy", ge=0)


class TileCache(_Section):
    dir: str = opt("", "Where the Galaxy Map's tiles are cached on disk; empty means /var/cache/planetgen/tiles.",
                   env="PLANETGEN_TILE_CACHE_DIR", editable=False, category="cache")
    max_mb: float = opt(200, "Roughly how big the tile cache may grow before its oldest files are pruned (0 turns "
                        "the disk cache off).", env="PLANETGEN_TILE_CACHE_MAX_MB", unit="MB", category="cache", ge=0)
    warm_levels: int = opt(8, "How many levels of the Galaxy Map's tile octree `planetgen warm` builds in full, "
                           "coarsest first (level 8 is tiles 2,048 pc across; every level is eight times the tiles "
                           "of the one above; 0 builds none).", unit="levels", category="cache", ge=0, le=15)
    warm_zoom_levels: int = opt(14, "Finest tile level of the zoom run `planetgen warm` builds from the opening "
                                "view in to the galaxy's centre, on top of the levels above (level 14 is tiles 32 pc "
                                "across and the finest the map's cut uses; at most 15; 0 builds no run).",
                                unit="levels", category="cache", ge=0, le=15)
    prefetch_per_visit: int = opt(128, "Most tiles one request for tiles queues to be built in the background: the "
                                  "next few zoom steps along the camera's line, then the neighbours of the tiles it "
                                  "asked for and the next level inside them (0 turns the background build off).", unit="tiles", category="cache", ge=0)
    prefetch_per_minute: int = opt(600, "Most tiles the background build queues in a minute, all requests together.",
                                   unit="tiles", category="cache", ge=0)
    prefetch_max_queued: int = opt(8, "The background build queues nothing while this many of its jobs (16 tiles "
                                   "each) or other bulk work already wait.", unit="jobs", category="cache", ge=1)


class PageCache(_Section):
    enabled: bool = opt(True, "Keep the API's public answers in memory, per WSGI process.",
                        env="PLANETGEN_PAGE_CACHE", restart=True, category="cache")
    max_entries: int = opt(2000, "Most answers kept.", restart=True, category="cache", ge=1)
    max_mb: float = opt(64, "Most memory the answers may take.", restart=True, unit="MB", category="cache", gt=0)
    stamp_seconds: float = opt(15, "How often the galaxy's content stamp is checked.", restart=True, unit="seconds",
                               category="cache", ge=0)
    max_age_seconds: float = opt(300, "Nothing is kept longer than this.", restart=True, unit="seconds",
                                 category="cache", gt=0)
    shared_seconds: float = opt(21600, "How long an answer is kept in Redis, where every web process and the warm-up "
                                "share it. Each answer is served only while the data it read has not changed, so "
                                "this only bounds the memory.", restart=True, unit="seconds", category="cache", gt=0)
    warm_max_objects: int = opt(5000, "Most sector pages, system pages and phenomenon pages the warm-up builds "
                                "(each, nearest the galaxy centre first); the rest are built on their first visit, "
                                "and the warm-up says how many it left out.", unit="objects", category="cache", ge=0)
    stored_counts: bool = opt(True, "Table pages show stored totals and filter-menu counts, refreshed in the "
                              "background, instead of counting whole tables in the request (off: exact counts "
                              "at any price).", env="PLANETGEN_COUNT_CACHE", category="cache")

    _loose = field_validator("enabled", "stored_counts", mode="before")(_loose_bool)


class Jobs(_Section):
    dir: str = opt("", "Where each background job's command lines, status and output are kept; empty means "
                   "/var/lib/planetGen/jobs.", env="PLANETGEN_JOBS_DIR", editable=False, category="jobs")
    keep: int = opt(20, "How many finished jobs are kept.", category="jobs", ge=1)
    python: str = opt("", "The interpreter that runs planetgen and planetgen.cli.reset for the jobs; empty means "
                      "the web app's own Python. Editing it from the web would be code execution, so it is file-only.",
                      env="PLANETGEN_PYTHON", editable=False, category="jobs")


class WikiJs(_Section):
    base_url: str = opt("", "The Wiki.js instance's root URL. Empty means it isn't offered as an upload target.",
                        env="PLANETGEN_WIKIJS_BASE_URL", role="owner", category="wiki")
    api_token: str = opt("", "A Wiki.js Personal API Token (Admin, API Access).", env="PLANETGEN_WIKIJS_API_TOKEN",
                         secret=True, role="owner", category="wiki")


class MediaWiki(_Section):
    base_url: str = opt("", "The MediaWiki API entry point's directory (everything up to, not including, api.php). "
                        "Empty means it isn't offered as an upload target.", env="PLANETGEN_MEDIAWIKI_BASE_URL",
                        role="owner", category="wiki")
    username: str = opt("", "A Bot Password user name, in User@BotName form.", env="PLANETGEN_MEDIAWIKI_USERNAME",
                        role="owner", category="wiki")
    password: str = opt("", "That Bot Password.", env="PLANETGEN_MEDIAWIKI_PASSWORD", secret=True, role="owner",
                        category="wiki")


class Wiki(_Section):
    wikijs: WikiJs = WikiJs()
    mediawiki: MediaWiki = MediaWiki()


class SettingsFile(_Section):
    backups_kept: int = opt(20, "How many backup archives are kept. One is made just before every settings save, "
                            "restore or import replaces the file, and by Back up now; the oldest is deleted once "
                            "there are more than this, but never the newest archive that holds every part.",
                            unit="backups", category="backups", ge=1, le=200)


_HEX_COLOR = re.compile(r"^#[0-9a-fA-F]{6}$")
_ASSET_NAME = re.compile(r"^([0-9a-f]{16}\.(png|ico|jpg|svg))?$")
_TWITTER_HANDLE = re.compile(r"^(@[A-Za-z0-9_]{1,15})?$")
_LANGUAGE = re.compile(r"^[a-z]{2,3}(-[A-Za-z]{2,4})?$")

ROBOTS_DIRECTIVES = ("user-agent", "allow", "disallow", "sitemap", "crawl-delay", "host")
"""tuple[str]: The lines a robots.txt text may hold besides comments and blank lines (RFC 9309 and the common extras)."""


def robots_problems(text):
    """`[message]` for each line of a robots.txt text that is not a `Field: value` line with a known field, a comment
    or blank; empty when the text is fine. The first problems only, so one bad paste doesn't list a thousand."""
    problems = []
    for number, raw in enumerate(text.splitlines(), 1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        field, colon, value = line.partition(":")
        if not colon or field.strip().lower() not in ROBOTS_DIRECTIVES:
            problems.append(f"line {number} is not a robots.txt line (Field: value; known fields are "
                            f"{', '.join(ROBOTS_DIRECTIVES)})")
            if len(problems) >= 5:
                break
    return problems


class Site(_Section):
    tagline: str = opt("", "A short line shown under the home page's heading and used as the description of pages "
                       "that give none of their own.", category="site", max_length=120)
    language: str = opt("en", "The pages' language as a code (en, de, en-GB): the `lang` of every page and the Open "
                        "Graph locale.", category="site")

    @field_validator("language")
    @classmethod
    def _language_code(cls, value):
        if not _LANGUAGE.match(value):
            raise ValueError("must be a language code such as en or en-GB")
        return value


class Web(_Section):
    title_template: str = opt("{title} - {site_name}", "How a page's title is written in the browser tab and in "
                              "search results; {title} is the page's own title and {site_name} the site name.",
                              category="web", max_length=120)
    theme_color_light: str = opt("#f6f7fb", "The browser's toolbar colour on phones for the light theme, as #rrggbb.",
                                 category="web")
    theme_color_dark: str = opt("#14151e", "The same for the dark theme.", category="web")
    assets_dir: str = opt("", "Where uploaded icons and images are kept; empty means the assets folder in the state "
                          "folder.", env="PLANETGEN_ASSETS_DIR", editable=False, category="web")

    @field_validator("title_template")
    @classmethod
    def _title_has_title(cls, value):
        try:
            shown = value.format(title="X", site_name="Y")
        except (KeyError, IndexError, ValueError) as exc:
            raise ValueError(f"may only use {{title}} and {{site_name}} ({exc})") from exc
        if "X" not in shown:
            raise ValueError("must contain {title}")
        return value

    @field_validator("theme_color_light", "theme_color_dark")
    @classmethod
    def _color(cls, value):
        if not _HEX_COLOR.match(value):
            raise ValueError("must be a colour written #rrggbb")
        return value


class Icons(_Section):
    favicon_svg: str = opt("", "The site icon as an uploaded SVG; empty uses the built-in one. Set it by uploading on "
                           "the Images page.", category="icons")
    favicon_png: str = opt("", "A PNG or ICO icon for browsers that do not take an SVG; empty offers none. Set it by "
                           "uploading on the Images page.", category="icons")
    apple_touch_icon: str = opt("", "The 180 by 180 PNG iPhones and iPads use for a home-screen bookmark. Set it by "
                                "uploading on the Images page.", category="icons")
    manifest_enabled: bool = opt(False, "Serve /site.webmanifest, which names the site, its icon and its theme colour "
                                 "for browsers that can install it.", category="icons")

    _loose = field_validator("manifest_enabled", mode="before")(_loose_bool)

    @field_validator("favicon_svg", "favicon_png", "apple_touch_icon")
    @classmethod
    def _asset(cls, value):
        if not _ASSET_NAME.match(value):
            raise ValueError("is the name of an uploaded file; upload it on the Images page")
        return value


class Seo(_Section):
    indexing: Literal["off", "top_pages", "all"] = opt(
        "top_pages", "Whether search engines may index the site: off (robots.txt disallows everything and every page "
        "says noindex), top_pages (the home, class, species, polity, galaxy and list pages), or all (also the "
        "system, sector and phenomenon pages).", category="seo")
    detail_pages: Literal["noindex", "index"] = opt(
        "noindex", "Whether the generated system, sector and phenomenon pages may be indexed when indexing is "
        "top_pages. Millions of generated pages look like scaled content to search engines, so they stay out by "
        "default; a small galaxy can switch this to index.", category="seo")
    description: str = opt("", "The description search engines show for pages that give none of their own (about 155 "
                           "characters are shown); empty uses the tagline.", category="seo", max_length=300)
    keywords: List[str] = opt([], "Keywords for the keywords meta tag. Google has ignored it since 2009; some "
                              "smaller engines still read it.", category="seo")
    json_ld: bool = opt(True, "Add a WebSite block (schema.org) to the home page.", category="seo")
    robots_extra: str = opt("", "Lines added to the end of the generated robots.txt.", category="seo", multiline=True,
                            max_length=64 * 1024)
    robots_override: str = opt("", "Replaces the generated robots.txt entirely when not empty. A mistake here can "
                               "hide the whole site from search engines.", category="seo", multiline=True,
                               max_length=500 * 1024)
    block_ai_crawlers: bool = opt(False, "Disallow the crawlers that collect text to train AI models (GPTBot, "
                                  "ClaudeBot, CCBot, Google-Extended, Applebot-Extended) in robots.txt.",
                                  category="seo")
    sitemap_enabled: bool = opt(True, "Serve /sitemap.xml and list it in robots.txt.", category="seo")
    sitemap_max_urls: int = opt(50000, "Most addresses in one sitemap file (the limit is 50,000).", unit="addresses",
                                category="seo", ge=1, le=50000)
    sitemap_cache_hours: int = opt(24, "How long a built sitemap is reused.", unit="hours", category="seo", ge=1,
                                   le=168)
    security_txt: str = opt("", "The text of /.well-known/security.txt (RFC 9116 needs a Contact and an Expires "
                            "line); empty serves none.", category="seo", multiline=True, max_length=16 * 1024)
    humans_txt: str = opt("", "The text of /humans.txt; empty serves none.", category="seo", multiline=True,
                          max_length=16 * 1024)
    privacy_note: str = opt("", "A short note shown on the account page about what the site keeps.", category="seo",
                            multiline=True, max_length=4 * 1024)

    _loose = field_validator("json_ld", "block_ai_crawlers", "sitemap_enabled", mode="before")(_loose_bool)

    @field_validator("robots_extra", "robots_override")
    @classmethod
    def _robots(cls, value):
        problems = robots_problems(value)
        if problems:
            raise ValueError("; ".join(problems))
        return value

    @field_validator("keywords")
    @classmethod
    def _keywords(cls, value):
        cleaned = [word.strip() for word in value if word.strip()]
        if len(", ".join(cleaned)) > 300:
            raise ValueError("the keywords together may be at most 300 characters")
        return cleaned


class OpenGraph(_Section):
    enabled: bool = opt(True, "Add Open Graph and Twitter tags to the public pages, so a link shows a title, a "
                        "description and a picture when shared.", category="og")
    site_name: str = opt("", "The site name in link previews; empty uses the site name.", category="og",
                         max_length=100)
    image_alt: str = opt("", "A description of the default preview picture for people who cannot see it.",
                         category="og", max_length=200)
    default_image: str = opt("", "The picture link previews show (a PNG or JPEG about 1200 by 630); empty uses the "
                             "built-in one. Set it by uploading on the Images page.", category="og")
    twitter_card: Literal["summary", "summary_large_image"] = opt(
        "summary_large_image", "How X (Twitter) draws a preview: a large picture or a small one.", category="og")
    twitter_site: str = opt("", "The site's X (Twitter) account, as @name; empty names none.", category="og")
    generated_cards: bool = opt(False, "Draw a preview picture per system and sector instead of one for the whole "
                                "site. Off until a card drawer is installed; every page then uses the "
                                "default picture.", category="og")
    card_cache_mb: float = opt(100, "How much disk the drawn preview pictures may use.", unit="MB", category="og",
                               ge=0)
    card_cache_dir: str = opt("", "Where drawn preview pictures are cached; empty means the state folder.",
                              editable=False, category="og")

    _loose = field_validator("enabled", "generated_cards", mode="before")(_loose_bool)

    @field_validator("default_image")
    @classmethod
    def _image(cls, value):
        if not _ASSET_NAME.match(value):
            raise ValueError("is the name of an uploaded file; upload it on the Images page")
        return value

    @field_validator("twitter_site")
    @classmethod
    def _handle(cls, value):
        if not _TWITTER_HANDLE.match(value):
            raise ValueError("must be an @name of up to 15 letters, digits or underscores")
        return value


class Settings(_Section):
    """Every `config.json` option."""

    config_version: int = opt(CONFIG_VERSION, "The settings file format version; a file without one is version 1.",
                              editable=False, category="file", ge=1)
    site_name: str = opt("planetGen", "Display name for this deployment, shown in the header, the home page "
                         "heading, every page title and the footer.", category="site", min_length=1)
    base_url: str = opt("http://localhost/", "The address this site is reached at, as an absolute http or https "
                        "address with a trailing slash. Canonical links, link previews, the sitemap and emailed links "
                        "are built from it, never from the request's Host header.", role="owner", category="site")

    @field_validator("base_url")
    @classmethod
    def _base_url(cls, value):
        parts = urlsplit(value)
        if parts.scheme not in ("http", "https") or not parts.netloc or parts.query or parts.fragment \
                or parts.username or parts.password:
            raise ValueError("must be an absolute http or https address such as https://example.org/ (no query, "
                             "fragment or password)")
        return value if value.endswith("/") else value + "/"

    api_base_url: str = opt("http://127.0.0.1/api", "Base URL of the Flask API's /api mount point, used by "
                            "web/lib/apiclient.py only outside the Flask app.", env="PLANETGEN_API_BASE_URL",
                            restart=True, editable=False, category="site")
    debug: bool = opt(False, "Write a verbose debug log to log_file: every decision the generator makes, every "
                      "random roll, every SQL statement and request. The log grows fast.", env="PLANETGEN_DEBUG",
                      restart=True, editable=False, category="logging")
    log_file: str = opt(logpaths.DEFAULT_LOG_FILE, "Where the debug log goes.",
                        env="PLANETGEN_LOG_FILE", restart=True, editable=False, category="logging")
    log_dir: str = opt("", "The folder of the always-on activity log; empty means /var/log/planetgen on Linux, "
                       "/Library/Logs/planetgen on macOS.",
                       env="PLANETGEN_LOG_DIR", restart=True, editable=False, category="logging")
    log_rotation: Literal["auto", "system", "app"] = opt(
        logpaths.DEFAULT_LOG_ROTATION, "How the activity log is rotated: system (logrotate or newsyslog), app (the "
        "program, 100 MB per file, 30 old copies) or auto (system where an installed rotation config exists, else "
        "app).", restart=True, editable=False, category="logging")
    settings: SettingsFile = SettingsFile()
    site: Site = Site()
    web: Web = Web()
    icons: Icons = Icons()
    seo: Seo = Seo()
    og: OpenGraph = OpenGraph()
    mysql: MySQL = MySQL()
    control_database: str = opt("planetgen_control", "The MySQL schema holding admin identities, sessions, API keys "
                                "and the audit log. Never listed or selectable as a galaxy.",
                                env="PLANETGEN_CONTROL_DATABASE", restart=True, editable=False, category="database")
    redis: Redis = Redis()
    ratelimit: RateLimit = RateLimit()
    login_allowlist: List[str] = opt(
        [], "Client addresses or networks (CIDR) never locked out after failed logins; loopback never is either. "
        "Environment: comma- or space-separated.", env="PLANETGEN_LOGIN_ALLOWLIST", restart=True, role="owner",
        category="security")
    admin_cookie_insecure: bool = opt(False, "Send the admin session cookie over plain HTTP. Only for local "
                                      "development without TLS; a production deployment must never set it.",
                                      env="PLANETGEN_ADMIN_COOKIE_INSECURE", restart=True, editable=False,
                                      category="security")
    secret_key: str = opt("", "Signs the CSRF tokens on the pages' forms. Empty makes a random one at startup (a "
                          "form open across a restart then fails once).", env="PLANETGEN_SECRET_KEY", restart=True,
                          secret=True, editable=False, category="security")
    proxy_fix: ProxyFix = ProxyFix()
    tile_cache: TileCache = TileCache()
    page_cache: PageCache = PageCache()
    jobs: Jobs = Jobs()
    wiki: Wiki = Wiki()

    _loose = field_validator("debug", "admin_cookie_insecure", mode="before")(_loose_bool)

    @field_validator("login_allowlist", mode="before")
    @classmethod
    def _split_allowlist(cls, value):
        if isinstance(value, str):
            return value.replace(",", " ").split()
        return value


class SettingsError(ValueError):
    """A settings file or environment variable that doesn't fit the model;
    the message names each field."""


# ---------------------------------------------------------------------
# Walking the model
# ---------------------------------------------------------------------

def _model_of(annotation):
    if isinstance(annotation, type) and issubclass(annotation, BaseModel):
        return annotation
    return None


def iter_fields(model=Settings, prefix=""):
    """
    Yields `(dotted path, FieldInfo, metadata dict)` for every leaf option
    of `model` (sections are walked into). The metadata always has every
    `x-` key (`opt`'s defaults), `description` and the annotation.
    """
    for name, info in model.model_fields.items():
        path = f"{prefix}{name}"
        section = _model_of(info.annotation)
        if section is not None:
            yield from iter_fields(section, f"{path}.")
            continue
        meta = dict(info.json_schema_extra or {})
        meta["description"] = info.description or ""
        meta["annotation"] = info.annotation
        yield path, info, meta


def field_default(info):
    return info.get_default(call_default_factory=True)


def type_name(annotation):
    """`"integer"`, `"text"`, `"list of text"`, ... for the docs table."""
    origin = get_origin(annotation)
    if origin is Literal:
        return " or ".join(f'"{a}"' for a in get_args(annotation))
    if origin in (list, List):
        return "list of text"
    return {bool: "true or false", int: "whole number", float: "number", str: "text"}.get(annotation, str(annotation))


def env_names():
    """Every environment variable the model names, as `{variable: path}`."""
    return {meta["x-env"]: path for path, _info, meta in iter_fields() if meta.get("x-env")}


def env_name(path):
    """The environment variable the model names for option `path`."""
    for found, _info, meta in iter_fields():
        if found == path:
            return meta["x-env"]
    raise KeyError(path)


# ---------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------

DEFAULT_STATE_DIR = "/var/lib/planetGen/state"
"""str: The web-owned state folder (the installers give it to Apache's user; `planetgen.admin.maintenance`
keeps its switch there too). The install folder above it is root's, so a file beside `config.json` could be
read by the web app but never replaced."""


def state_dir():
    """The web-owned state folder: `PLANETGEN_STATE_DIR`, else `/var/lib/planetGen/state`."""
    return os.environ.get("PLANETGEN_STATE_DIR") or DEFAULT_STATE_DIR


def _read_json_object(path, label):
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
    except FileNotFoundError:
        return {}
    except (OSError, ValueError) as exc:
        raise SettingsError(f"{label} ({path}) can't be read: {exc}") from exc
    if not isinstance(data, dict):
        raise SettingsError(f"{label} ({path}) must contain a JSON object, not {type(data).__name__}")
    return data


def read_config():
    """The `config.json` dict (empty when there is no file); raises `SettingsError` for one that can't be read."""
    return _read_json_object(logpaths.CONFIG_PATH, "config.json")


def migrate(data):
    """
    Brings a file's dict up to `CONFIG_VERSION` (no migrations exist yet:
    version 1 is the first).

    Raises:
        SettingsError: The file is from a newer release.
    """
    version = data.get("config_version", 1)
    if not isinstance(version, int) or isinstance(version, bool) or version < 1:
        raise SettingsError(f"config_version must be a whole number 1 or above, not {version!r}")
    if version > CONFIG_VERSION:
        raise SettingsError(f"config_version {version} is newer than this release understands ({CONFIG_VERSION}); "
                            "update planetGen")
    return data


def _merge(base, overrides):
    for key, value in overrides.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            _merge(base[key], value)
        else:
            base[key] = copy.deepcopy(value)  # never share a layer's nested dict with the merge


def _set_path(data, path, value):
    parts = path.split(".")
    for part in parts[:-1]:
        if not isinstance(data.get(part), dict):
            data[part] = {}
        data = data[part]
    data[parts[-1]] = value


def _environment_layer(environ):
    layer = {}
    for path, _info, meta in iter_fields():
        name = meta.get("x-env")
        if not name or name not in environ:
            continue
        value = environ[name]
        if value == "" and not meta.get("x-env-empty"):
            continue
        _set_path(layer, path, value)
    return layer


def build(config=None, environ=None, strict=False):
    """
    Validates the layers into a `Settings`: defaults, `config` (the
    `config.json` dict), then the `PLANETGEN_*` variables in `environ`.

    Args:
        strict (bool): Unknown keys in `config` are errors (`config check`)
            instead of ignored.

    Raises:
        SettingsError: Naming every field that doesn't fit.
    """
    environ = os.environ if environ is None else environ
    merged = {}
    problems = []
    if config is not None:
        try:
            layer = migrate(dict(config))
        except SettingsError as exc:
            raise SettingsError(f"config.json: {exc}") from exc
        if strict:
            known = {path for path, _i, _m in iter_fields()}
            problems.extend(f"config.json: {message}" for message in _unknown_paths(layer, known))
        _merge(merged, layer)
    _merge(merged, _environment_layer(environ))
    if problems:
        raise SettingsError("\n".join(problems))
    try:
        return Settings.model_validate(merged)
    except ValidationError as exc:
        lines = [f"{path}: {message}" for path, message in _error_pairs(exc)]
        raise SettingsError("; ".join(lines) if len(lines) < 4 else "\n".join(lines)) from exc


def _error_pairs(exc):
    """`[(dotted path, message)]` of a `ValidationError`, the path cut back to the option it names (a list
    option's element error is the list's)."""
    leaves = {path for path, _i, _m in iter_fields()}
    pairs = []
    for error in exc.errors():
        parts = [str(part) for part in error["loc"]]
        path = next((".".join(parts[:n]) for n in range(len(parts), 0, -1) if ".".join(parts[:n]) in leaves),
                    ".".join(parts))
        pairs.append((path, error["msg"]))
    return pairs


def field_errors(config=None, environ=None):
    """
    Validates `config` like `build` but without raising and without the unknown-key check: `(settings, [])` when
    it fits, else `(None, [(dotted option path, message), ...])`. The environment is left out unless `environ`
    is given. The writer uses it to put each message at its option.
    """
    environ = {} if environ is None else environ
    merged = {}
    if config:
        _merge(merged, migrate(dict(config)))
    _merge(merged, _environment_layer(environ))
    try:
        return Settings.model_validate(merged), []
    except ValidationError as exc:
        return None, _error_pairs(exc)


def _unknown_paths(layer, known, prefix=""):
    sections = {path.rsplit(".", 1)[0] for path in known if "." in path}
    out = []
    for key, value in layer.items():
        path = f"{prefix}{key}"
        if path in known:
            continue
        if isinstance(value, dict) and (path in sections or any(k.startswith(f"{path}.") for k in known)):
            out.extend(_unknown_paths(value, known, f"{path}."))
        else:
            out.append(f"{path} is not a setting")
    return out


_lock = threading.Lock()
_cache = {"key": None, "settings": None}


def _stamp(path):
    try:
        stat = os.stat(path)
    except OSError:
        return None
    return (stat.st_mtime_ns, stat.st_size)


def _environment_stamp(environ):
    return tuple(sorted((name, environ[name]) for name in env_names() if name in environ))


def get_settings(environ=None):
    """
    The validated settings. Re-read only when `config.json` or a variable
    the model names changes (a stat costs about a microsecond), so every
    option that is read per use stays live.

    Raises:
        SettingsError: A file or variable doesn't fit the model.
    """
    environ = os.environ if environ is None else environ
    config_path = logpaths.CONFIG_PATH
    key = (config_path, _stamp(config_path), _environment_stamp(environ))
    with _lock:
        if _cache["key"] == key:
            return _cache["settings"]
        settings = build(_read_json_object(config_path, "config.json"), environ)
        _cache["key"], _cache["settings"] = key, settings
        return settings


def reset_cache():
    """Forgets the cached settings (tests)."""
    with _lock:
        _cache["key"] = _cache["settings"] = None


# ---------------------------------------------------------------------
# What the model writes: the docs table, the example file, the schema
# ---------------------------------------------------------------------

def example_dict():
    """`config.json.example`'s content: every default of the options
    `config.json` may hold, in the model's order."""
    return Settings().model_dump(mode="json")


def example_json():
    return json.dumps(example_dict(), indent=2) + "\n"


def schema_json():
    """`config.schema.json`: the model's JSON Schema."""
    return json.dumps(Settings.model_json_schema(), indent=2, sort_keys=True) + "\n"


def _format_default(value):
    if isinstance(value, bool):
        return "true" if value else "false"
    if value == "":
        return "empty"
    if isinstance(value, list):
        return "empty list" if not value else ", ".join(str(v) for v in value)
    return str(value)


def docs_table():
    """The generated option table for `docs/config.md` (Markdown)."""
    rows = ["| Option | Type | Default | Applies | Edited | Environment | What it does |",
            "|---|---|---|---|---|---|---|"]
    roles = {"admin": "any admin", "owner": "the Owner"}
    for path, info, meta in iter_fields():
        applies = "after a restart" if meta["x-restart"] else "at once"
        edited = roles[meta["x-min_role"]] if meta["x-editable"] else "command line only"
        unit = f" ({meta['x-unit']})" if meta.get("x-unit") else ""
        secret = " **Secret.**" if meta["x-secret"] else ""
        env = f"`{meta['x-env']}`" if meta.get("x-env") else "none"
        rows.append(f"| `{path}` | {type_name(meta['annotation'])}{unit} | {_format_default(field_default(info))} | "
                    f"{applies} | {edited} | {env} | {meta['description']}{secret} |")
    return "\n".join(rows) + "\n"


def settings_hash(settings=None):
    """A short hash of the effective settings, secrets left out (for a
    page that detects a concurrent edit)."""
    settings = settings or get_settings()
    secret_paths = {path for path, _i, meta in iter_fields() if meta["x-secret"]}
    flat = {path: _lookup(settings, path) for path, _i, _m in iter_fields() if path not in secret_paths}
    return hashlib.sha256(json.dumps(flat, sort_keys=True, default=str).encode()).hexdigest()[:16]


def _lookup(settings, path):
    node: Any = settings
    for part in path.split("."):
        node = getattr(node, part)
    return node
