# Settings model, configuration page, SEO and user accounts

How `config.json` becomes one validated settings model with a web configuration page, how the site's icons, Open Graph tags, robots.txt and sitemap are configured, and how user accounts, invites, password reset, SMTP and the Owner role are built on the control database that already holds the admin accounts. Part A covers the settings model, the Admin page and SEO; Part B covers accounts and security. Boss's uploaded "Web UX and Job Management Guide.md" and "Web UX Development Notes.md" were read for anything on auth or settings and have none; their queue sections describe an SQLite design that [library-migration.md](library-migration.md) replaced with RQ on Redis.

Informs: ADM.42, ADM.43, ADM.44, ADM.18, USR.1, USR.2, USR.3, USR.4, USR.5, USR.6, USR.7, USR.8, API.6, VIEW.3 (card renderer hook)

Status: research, 2026-10-09; decisions marked "Boss" are his, everything else is a recommendation

Evidence tags: [S] seen in a search result or package page, [C] computed in the research session (Python 3.9.25, lock-file versions; timings are from a slow, noisy VM, so only ratios carry over), [R] recalled and unconfirmed. The research environment could read search-result text only, not the papers or the vendors' own pages.

## Decisions already taken

- Configuration page (Boss, GitHub issue #515, 2026-10-08 00:39Z): "A page under admin is needed to configure everything EXCEPT for the database information. All other settings should be accessed and changeable from there." (ADM.43)
- Web settings (Boss, issue #743, 2026-10-09 02:55Z): icons, Open Graph fields for Discord and other link previews, SEO description and keywords, and editable web documents such as robots.txt (ADM.44).
- Accounts (Boss, 2026-10-01): an email loop for setting and resetting passwords; invite only, by unique link, with a use count and an expiry of 1, 4, 6, 12 or 24 hours, 3 or 7 days, 1 month, 1 year or never; infinite uses with never expiring asks for confirmation. One "Owner" account, made at install, cannot be overridden. Owner transfer needs specific approval, double password confirmation and an email loop, and the target accepts by email. SMTP settings live in the admin config (USR.1 to USR.6).
- API keys (Boss, 2026-10-02 01:31Z): "all keys are attached to the account that created them"; admins create user-level keys that read but do not upload (API.6). All signed-in users may generate a one-off system (Boss, 2026-10-02 05:16Z, USR.8).
- Library swaps leave no shims (Boss, 2026-10-07), so `util/appconfig.py` is deleted, not wrapped.
- SEC.29, SEC.30 and SEC.31 are done (PRs #612, #643, #457 per [todo-number-map.md](todo-number-map.md)).

## Part A: settings model, configuration page, SEO and Open Graph

### A1. What exists today and the library choice

Verified against the code on this branch:

- `util/appconfig.py` holds `DEFAULT_CONFIG` (a dict) and `load_config()`, which re-reads and deep-merges `config.json` on every call with no schema, no unknown-key check and no type checks. A typo such as `"ratelimt"` is accepted silently, and `debug` needs a hand-written string-to-bool shim. A Pydantic bool does that (`"false"` and `"0"` give `False`, `"maybe"` is an error) [C].
- About 14 modules import `appconfig`: `api/config.py` freezes most options at import, the rest read per use.
- `docs/config.md` has a hand-written field table of about 25 rows. Keep its prose, generate the table (A8). `base_url` is documented as "not yet read by any page"; ADM.44 makes it load-bearing, so it needs validation.
- `base.html` hard-codes the two `theme-color` values and `favicon.svg`.
- Pydantic is already a dependency (`pydantic>=2.7.0`; `api/schemas.py`, `web/generate_page.py`). `setup.py` says `python_requires='>=3.9'`, so each new dependency needs a 3.9-compatible release.

Library facts (PyPI, 2026-10-09 [S]): pydantic 2.14.0 and pydantic-settings 2.15.0 need Python 3.10, and the last releases for 3.9 are pydantic 2.13.5 (the lock pin) and pydantic-settings 2.11.0. The same split holds for Flask-Limiter (3.11.0 for 3.9; limits 4.2), `webauthn` (2.7.1) and Flask-WTF (1.2.2); `argon2-cffi` 25.1.0 and `pyotp` 2.10.0 support 3.8 and later.

pydantic-settings adds python-dotenv and its environment-name mapping does not fit the `PLANETGEN_*` names: with `env_nested_delimiter="_"`, `PLANETGEN_TILE_CACHE_DIR` was applied but `PLANETGEN_TILE_CACHE_MAX_MB` was silently ignored [C]. Declare each option's variable in its field metadata instead. On Python 3.9 write `Optional[str]`, not `str | None`, which raises `TypeError` inside a model unless `eval_type_backport` is installed [C].

### A2. The model, and which file holds what

Nested `BaseModel` sections with `extra="forbid"` when validating a form or import and `extra="ignore"` plus one logged warning when loading (an older server must not crash on a newer file): `config_version`, `site`, `web`, `seo`, `og`, `icons`, `security`, `smtp`, `limits`, `logging`, `ratelimit`, `cache`, `jobs`, `wiki`, `redis`, `proxy_fix`, and `mysql` and `control_database` (never editable).

Field metadata goes in `json_schema_extra` under `x-` keys, declared through one helper; the keys survive `model_json_schema()` and `Field(ge=, le=, min_length=, pattern=)` become `minimum`, `maximum` and so on [C].

| Key | Meaning |
|---|---|
| `description` (Pydantic) | help text for the form, the docs table and the schema |
| `x-unit` | suffix on the input and in docs |
| `ge`, `le`, `min_length`, `max_length`, `pattern` (Pydantic) | validation and the HTML `min`, `max`, `maxlength` |
| `x-restart` | change applies only after a restart |
| `x-secret` | never sent to the browser or logged; redacted in export and audit |
| `x-category` | form section |
| `x-editable` | `false`: file-only, shown read-only |
| `x-min_role` | `admin` or `owner`: who may save it |
| `x-env` | environment variable that overrides it; field shown disabled |
| `x-widget` | optional override (`textarea`, `color`, `file`) |

**The permissions problem.** `docs/config.md` sets `config.json` to `root:<apache group>`, mode 640: the web process reads it but cannot write it, and an atomic replace also needs write access to the directory. ADM.43 cannot save as written. Two files, one model:

| File | Owner and mode | Holds |
|---|---|---|
| `config.json` | root, group read, 640 | database, control database, `secret_key`, the host-level tier |
| `settings.json` | web user, 0600, in `/var/lib/planetGen` or the platform equivalent (next to `jobs.dir`) | everything the page may change |

Load order, lowest to highest: model default, `config.json`, `settings.json`, `PLANETGEN_*` environment, explicit argument. A key in the wrong file (a database key in `settings.json`) is an error. Making `config.json` web-writable instead would give any web-process compromise the database password. The installers and `set-permissions.sh` must create the settings folder for the web user, as `examples/apache/create-cache-dir.sh` does for the cache and jobs folders.

**Three tiers.** Some options are code execution or request-forgery paths if an admin account is stolen: `jobs.python` is the interpreter the web app launches, so editing it from the browser is remote code execution as the web user.

| Tier | Options | Edited in the browser by |
|---|---|---|
| admin | site text, SEO, Open Graph, icons, `limits`, `ratelimit` numbers, tile cache, `jobs.keep`, page cache | any admin |
| Owner | SMTP, `base_url`, wiki URLs and credentials | the Owner, after re-entering the password |
| read-only | `mysql.*`, `control_database`, `jobs.python`, `jobs.dir`, log paths, `secret_key`, `proxy_fix.*`, `admin_cookie_insecure`, `redis.url`, `api_base_url` | nobody; the page shows the value and "edit config.json on the server" |

This enforces "everything except the database" by file permission, not only by hiding fields.

**ADM.18 stays separate.** The galaxy-creation JSON must reproduce a galaxy regardless of deployment. Give it its own model (`GalaxyCreationSettings`, own `schema_version`) with the same metadata conventions, never inside `config.json` or `settings.json` (see [reproducible-galaxies.md](reproducible-galaxies.md)).

### A3. The Admin form generated from the model

A Jinja macro walks `Settings.model_fields` and picks a Shoelace widget (UX.49) from the annotation: `sl-switch` for `bool`; `sl-input type=number` with `min`, `max`, `step` and a unit suffix for numbers; `sl-input` for `str`; `sl-select` for `Literal` and `Enum`; `sl-textarea` one per line for `List[str]` (login allowlist) and for long text (robots.txt); a password input where empty means keep, plus a "clear" box, for secrets; `sl-color-picker`; and an upload for icons and the OG image (the setting stores the file name).

- **Field names are dotted paths** (`tile_cache.max_mb`). `ValidationError.errors()` returns `loc` tuples that join to exactly those names, so errors map to fields with no extra code: `site_name | string_too_short`, `debug | bool_parsing`, `mysql.port | less_than_equal`, `unknown | extra_forbidden` [C].
- **Re-render the whole form** with submitted values and per-field messages in each control's `help-text` slot, with `aria-invalid`. That works without JavaScript; UX.49's browser test should submit one bad value to prove it.
- **Missing versus false.** An unchecked checkbox sends nothing, so post a hidden list of the rendered field names and apply only those. Disabled (environment-overridden) and read-only fields are never applied.
- **Validation layers:** per field (Pydantic); cross field (`model_validator`: `smtp.security=none` with a non-loopback host is refused); environmental checks that only warn (folder writable, Redis answers, SMTP test, PNG decodes), run from "Test" buttons, not on every save.
- **Concurrent edits.** The form carries a hash of the file as read; a stale hash is refused with a diff.
- **Recent authentication** for the Owner tier: the password again (valid 5 minutes) before saving SMTP, `base_url` or wiki settings. Each field also shows where its value comes from (default, file or environment).

### A4. Applying changes: live reload versus restart

The code already splits this way by accident (from the `load_config` call sites, read not run [R]; confirm each when ADM.42 is built):

| Option | Read | Treat as |
|---|---|---|
| `site_name`, `tile_cache.*`, `jobs.*`, `mysql.statement_timeout_seconds`, `wiki.*`, `redis.url` (queue) | on each use | live |
| `debug`, `log_file`, `log_dir`, `log_rotation`, `page_cache.*` | at logger or app start | restart (could become live) |
| `ratelimit.*`, `login_allowlist`, `admin_cookie_insecure`, `secret_key`, `proxy_fix.*`, `api_base_url`, `mysql.*`, `control_database` | at import into `Config` | restart |

1. `get_settings()` returns an immutable validated object cached on `(mtime_ns, size)` of both files plus a hash of the relevant environment variables. A `stat` costs 1.3 microseconds and validating the full model 8.5 microseconds [C], so per-request checking is free and every live field becomes live with no extra code.
2. At start store `STARTUP = get_settings()`. The page lists `x-restart` fields whose current value differs from `STARTUP` ("saved but not yet running: `ratelimit.default` running `200 per day`, saved `500 per day`"). No extra state is persisted.
3. A bad hand-edited file at runtime logs once and keeps the last good settings; at startup it stops the app with per-field messages (the existing `proxy_fix` behaviour, generalized). Add `python -m planetgen.cli.config check`.
4. Restart. Under Apache and mod_wsgi daemon mode (the documented setup), touching the WSGI script file restarts the daemon on the next request [S: modwsgi.readthedocs.io]; Unix daemon mode only, and the web user must be able to write the script. gunicorn takes SIGHUP [R]; waitress has no reload. Show "Restart now" only when the deployment is detected and permitted, otherwise a copy-paste command (`update.sh` already prints `sudo systemctl reload apache2`), and say that a restart drops in-flight requests for a few seconds.
5. Never apply a changed `secret_key` live: it invalidates every CSRF token, flash cookie and signed token.

### A5. Secrets

Secrets in scope: SMTP password, wiki token and password, `secret_key`, MySQL password. SMTP and wiki secrets go in `settings.json`; `secret_key` and the MySQL password stay in `config.json`. An environment variable is a supported override (`x-env`), not storage (mod_wsgi does not pass Apache `SetEnv` into `os.environ`). Encrypting in the control database is not worth it: the key would sit on the same host where the web user reads it, so it only protects a leaked dump, and a file outside the database is already absent from dumps. Secrets are never returned to the browser (the form shows "set" and takes a replacement), never written to the audit or debug log (record `changed`), and redacted in export and in `config check`.

### A6. Saving, backup, rollback, versioning, audit, import and export

Save pipeline [C: 20 saves in 0.6 s with fsyncs, mode 0640 kept, backups pruned]:

1. Parse the form to a nested dict, merge it onto the current web-layer values, validate the whole model; if invalid, re-render and change nothing.
2. Show a diff (secrets as `changed`); require confirmation for restart-flagged and Owner-tier fields.
3. Copy the current file to `settings.json.bak-<UTC stamp>` and keep the newest 10.
4. Write a temp file in the same directory, copy mode and owner, `fsync`, `os.replace` (atomic on POSIX, and on Windows on one volume), `fsync` the directory.
5. Write audit rows and bump the cache.

Rollback lists the backups with diffs and a "restore" button that runs the same pipeline, so a restore is itself backed up and audited; one that fails validation under the current model is refused with the reason. Save only `exclude_defaults` values so the file stays short and new defaults reach existing installs.

Versioning: top-level `"config_version": 1` (today's file has none; missing means 1). Migrations are ordered `dict -> dict` functions applied in memory; the file is rewritten at the next save, since the loader may not own it. A test loads every committed fixture of every past version.

Audit: `admin_audit_log` (`action="config.update"`, `target="config:tile_cache.max_mb"`, `detail` JSON `{old,new}`, secrets as `"(changed)"`) and the activity log's `DB` category, as other admin writes do; also `config.restore`, `config.import`, `smtp.test`. The 90-day pruning of failed-login rows does not apply to these.

Import and export: export downloads the effective web-layer settings with secrets removed. Import runs the migration chain, validates, shows the diff and follows the save pipeline; keys the importer may not edit are rejected with a message. Also write `config.schema.json` from `model_json_schema()` (3.7 ms [C]).

### A7. Generated docs and the drift tests

`python -m planetgen.cli.config docs` rewrites the table between `<!-- BEGIN GENERATED -->` markers in `docs/config.md` (option, type, default, unit, live or restart, environment variable, description) and writes `config.json.example` from the defaults; the permissions, precedence and activity-log prose stays hand-written. Tests: (1) regenerating matches the committed files; (2) no module outside the settings module calls `load_config`, reads `DEFAULT_CONFIG` or indexes `config["..."]` (the ADM.42 requirement; typed access such as `settings.tile_cache.max_mb` makes a missing option an `AttributeError`, stronger than a grep); (3) every `x-env` name appears in the docs, no non-editable field is in the generated form, no secret field is in export.

### A8. ADM.44: web, Open Graph and SEO settings

All feed templates and routes at request time, so none needs a restart. Existing: `site_name`.

| Setting | Type, default | Notes |
|---|---|---|
| `site.tagline`, `site.language` | str <=120 `""`; `en` | home page title; `<html lang>`, `og:locale` |
| `web.base_url` | absolute http(s) URL, trailing slash, `http://localhost/` | Owner tier; canonical, OG, sitemap, email links, WebAuthn origin; never derived from the `Host` header |
| `web.title_template`, `web.theme_color_light`, `_dark` | `{title} - {site_name}`; `#f6f7fb`, `#14151e` | hard-coded in templates today |
| `web.assets_dir` | platform default | host-level, file-only |
| `icons.favicon_svg`, `favicon_png`, `apple_touch_icon` (180x180), `manifest_enabled` | uploads; false | built-in `favicon.svg`; the manifest serves `/site.webmanifest` |
| `seo.indexing` | `off` / `top_pages` / `all`, `top_pages` | `off` serves `Disallow: /` and `noindex` everywhere |
| `seo.detail_pages` | `noindex` / `index`, `noindex` | system, sector and phenomenon pages |
| `seo.description`, `seo.keywords` | str <=300; list, empty | description for pages that give none (about 155 characters show in snippets [R]); Google has ignored keywords since 2009 [S], and the help text says so |
| `seo.json_ld` | bool, true | home page `WebSite` block |
| `seo.robots_extra`, `seo.robots_override` | text <=64 KiB, <=500 KiB | appended to, or replacing, the generated robots.txt |
| `seo.block_ai_crawlers` | bool, false | `Disallow: /` for GPTBot, ClaudeBot, CCBot, Google-Extended, Applebot-Extended [S, vendor guides] |
| `seo.sitemap_enabled`, `sitemap_max_urls`, `sitemap_cache_hours` | true; 1..50000 (50000); 1..168 (24) | |
| `seo.security_txt`, `humans_txt`, `privacy_note` | text, empty | `/.well-known/security.txt` (RFC 9116 wants `Contact` and `Expires` [R]) and `/humans.txt` when non-empty; the note shows on the account page |
| `og.enabled`, `og.site_name`, `og.image_alt` | true; falls back to `site_name`; str <=200 | |
| `og.default_image` | PNG or JPEG, built-in 1200x630 | width and height read from the file |
| `og.twitter_card`, `og.twitter_site` | `summary_large_image`; `@handle` or empty | |
| `og.generated_cards`, `og.card_cache_mb`, `og.card_cache_dir` | false; 100; `""` | per-object PNG cards; reuse the tile-cache pruner |

`seo.robots_override` or `robots_extra` can break crawling site-wide, so validate: size cap, parseable `User-agent`, `Allow`, `Disallow` and `Sitemap` lines, and a warning when the result blocks `/` while `seo.indexing` is not `off`.

**Open Graph and Twitter/X tags** come from one block in `base.html`: `og:title`, `og:type` (`website`), `og:url` (the canonical URL), `og:description`, `og:site_name`, `og:locale`, `og:image` (absolute https URL) with `:width`, `:height`, `:alt` and `:type`, plus `twitter:card` and `twitter:image:alt`.

- X falls back to the `og:` tags when the Twitter ones are missing [S, guide sites]. One source says Discord shows only a small thumbnail without `twitter:card=summary_large_image` [S, single source; test with a real Discord post].
- Image: 1200x630; X crops to about 2:1, minimum 300x157, at most 5 MB; PNG or JPEG, because WebP support is uneven in Slack and Discord [S].
- Discordbot fetches when a link is posted and may keep the preview long; sources disagree on whether it obeys robots.txt [S]. Do not `Disallow` page URLs that should unfurl, and add `?v=<content hash>` to `og:image` so a changed card refetches.
- Admin, account and login pages emit no Open Graph tags.

**Canonical URLs and pagination.** `canonical` is `web.base_url` plus the path with tracking and sort parameters stripped (`needs_canonical_redirect` in `web/searchpage.py` already redirects some old forms). Do not emit `rel=prev/next`; Google dropped it years ago [S]. List pages (`/sectors`, `/systems`, `/species`, `/polities`): page 1 is indexable with a self canonical; `?page=N`, filters and sorts get `noindex,follow` and a self canonical, not a canonical to page 1, which would hide the linked items. The lists reach billions of rows (about 10.5 billion sector rows, see [galaxy-disk-density.md](galaxy-disk-density.md)), so add `Disallow: /*?*page=` to the generated robots.txt (`*` is in RFC 9309).

**Indexing matrix** (defaults):

| URL class | Robots meta or header | robots.txt | Sitemap |
|---|---|---|---|
| `/`, `/classes*`, `/species*`, `/polities*`, `/galaxy` shell | `index,follow` | allow | yes |
| `/sectors`, `/systems` (page 1, no query) | `index,follow` | allow | yes |
| list pages with `page`, `sort` or filter queries | `noindex,follow` | `Disallow` for deep `page=` | no |
| `/system/<id>`, `/sector/<id>`, `/phenomenon/...` | `noindex,follow` (`index` if `seo.detail_pages=index`) | allow, so the tag is seen | no (yes if `index`, capped) |
| `/search`, `/nav`, `/galaxy/tiles`, `/galaxy/stage`, `/galaxy/territories`, `/galaxy/locate`, `/table/*`, scene and shape JSON | `X-Robots-Tag: noindex` | `Disallow` | no |
| `/admin*`, `/login*`, `/logout`, `/account*`, `/api/*` | `noindex,nofollow` | `Disallow` | no |

- A robots.txt `Disallow` stops crawling, not indexing, and hides any `noindex` on that URL, so choose one per URL [S: Search Engine Journal, Semetrical].
- Google reads robots.txt up to 500 KiB, caches it about 24 hours, treats a persistent 5xx as "stop crawling" [S, translated snippets] and ignores `Crawl-delay` [S].
- Not indexing generated detail pages follows third-party summaries of Google's scaled-content guidance [S]; `seo.detail_pages` switches it for a small galaxy.
- The page limiter allows 300 requests a minute per IP, which a crawler can hit; keep robots.txt and the sitemap cached.

**Sitemap.** `/sitemap.xml` is an index listing `/sitemap-pages.xml` and, only if detail pages are indexed, `/sitemap-systems-<n>.xml`. Limits: 50,000 URLs and 50 MB uncompressed per file [S, secondary quoting sitemaps.org and Bing 2016]; at about 100 bytes a URL, 50,000 URLs is about 5 MB [C], so the count binds. Emit `loc` and `lastmod` only (Google ignores `changefreq` and `priority`, and trusts `lastmod` only when consistently accurate [S]), from a real change time or omitted. Build files by keyset pagination (`WHERE id > ? ORDER BY id LIMIT 50000`) so boundaries do not shift, cache them `seo.sitemap_cache_hours`, and list the sitemap in robots.txt (pings are deprecated [S]).

**JSON-LD.** Google retired the sitelinks search box in November 2024 [S, trade press], so emit only a `WebSite` block (`name`, `url`), and no `Place` or `Product` markup for generated objects. An inline `application/ld+json` block is not executed and should not violate the CSP's `script-src` [R: confirm in the Playwright tests].

**Link-preview images (ties to VIEW.3).** Default: one static admin-uploaded PNG. Per-object cards drawn with Pillow (name, star types, planet count, orbit sketch) cost 125 ms and 39 KiB per 1200x630 card on the slow VM [C], need no new dependency (Pillow arrives through scikit-image) and are cached on disk; a WebGL screenshot would need headless Chromium, so no. Route `/og/<kind>/<id>.png`: `image/png`, `Cache-Control: public, max-age=86400`, an `ETag` from the object's change stamp (`apiclient.get_galaxy_changes`), its own limiter bucket, 404 for unknown ids, text only from stored names (already filtered). `render_card(kind, object_id) -> bytes` sits behind a registry and falls back to the static image on error, so VIEW.3's sky render can register later; ADM.44 does not wait for it. `og.generated_cards` is off by default to avoid extra write load.

**Icon upload.** Accept PNG and ICO; verify with Pillow, check the pixel size, store under a content-hash name in `web.assets_dir`, serve with a long cache and `nosniff`. Accept SVG only with `Content-Security-Policy: default-src 'none'; sandbox` on the asset response, or not at all (a script in an uploaded SVG opened directly runs on the site's origin). External icon URLs would be blocked by the CSP, so uploads are the only form.

### A9. How it fits the code

New package `planetgen/settings/` (model, loader, writer, migrations, doc generator); `util/appconfig.py` is deleted once callers move. `api/config.py` keeps `Config` but reads `STARTUP`. The Admin page sits beside `admin_pages.py` behind `_admin_or_403` with `Cache-Control: no-store`. Sitemap, robots, security.txt, humans.txt, manifest and OG routes go in a new blueprint `web/site_documents.py`. `base.html` takes `theme-color`, favicon and title from the settings.

## Part B: user accounts and security practice

### B1. State of the code (what USR builds on)

- Tables (`db/control_schema.sql`): `admin_users` (username, password_hash, must_change_credentials, last_login_at), `admin_sessions` (token SHA-256, `expires_at`), `admin_api_keys`, `admin_audit_log` (username kept as text, foreign key `SET NULL`), `admin_devices`, `admin_totp`, `admin_recovery_codes`. No email column, no role. The control schema is at v10 (`CONTROL_SCHEMA_VERSION`, `db/store.py`) and migrations are Alembic revisions (DB.11, done), so the account columns are the next revision.
- `admin/auth.py`: `PASSWORD_HASH_METHOD = "pbkdf2:sha256:600000"`, `MIN_PASSWORD_LENGTH = 12`, a bundled blocklist, no maximum length, `SESSION_TTL_HOURS = 12` fixed, tokens from `secrets.token_urlsafe(32)` stored as SHA-256, `needs_rehash` at login, and a dummy hash so unknown usernames cost the same time.
- Session cookie: HttpOnly, Secure (unless the development flag), `SameSite=Strict`. `create_session` always mints a new token, so a pre-set cookie cannot be fixated. `web/csrf.py` is a nonce cookie plus an HMAC bound to the session, checked on every unsafe page request. Headers: `nosniff`, `X-Frame-Options: DENY`, `Referrer-Policy: no-referrer` (a reset link cannot leak its token by Referer), a strict CSP, HSTS over HTTPS.
- `api/loginguard.py` applies the per-IP lockout and per-username backoff of `admin/throttle.py`; the counters are the project's own `RedisStore`, on the same Redis server Flask-Limiter counts requests on. `record_login_failure` writes `login.failed` and `login.locked`. `cli/lockouts.py` already has `--reset-two-factor USER` and `--forget-devices USER`.

### B2. Password storage

| Scheme | OWASP position [S: Password Storage Cheat Sheet, snippets] | Measured here [C] | Fit |
|---|---|---|---|
| Argon2id | first choice; 19 MiB, t=2, p=1 minimum (46 MiB, t=1 equivalent) | 263 ms (19 MiB, t=2); 221 ms (46 MiB, t=1); 1240 ms (argon2-cffi defaults 64 MiB, t=3, p=4) | recommended |
| scrypt | N=2^17, r=8, p=1 | 2464 ms, 128 MiB (Werkzeug's N=2^15: 446-510 ms) | too heavy at 5 threads |
| PBKDF2-HMAC-SHA256 | 600,000 iterations, listed last, for FIPS-140 | 1345 ms | in use; no memory hardness |

The brute-force design measured PBKDF2-600k at about 175 ms on the build machine; the research VM is about 8 times slower, so only ratios matter. Argon2id at the OWASP minimum is faster than PBKDF2-600k, uses 19 MiB per hash (95 MiB across 5 threads) and resists GPU cracking, which matters once user accounts make a stolen control database the main offline risk.

Migration: `PasswordHasher(time_cost=2, memory_cost=19456, parallelism=1)`; `verify_password` dispatches on the stored prefix (`$argon2id$` versus `pbkdf2:` or `scrypt:`) because Werkzeug cannot verify Argon2; `needs_rehash` is also true for old schemes and weaker Argon2 parameters (`check_needs_rehash`), so every user upgrades at next login with the code that exists; the `VARCHAR(255)` column fits (about 97 characters); recreate the dummy hash with the active scheme. Add a 1,024-character maximum to `validate_password_policy` (none today; NIST asks verifiers to accept at least 64 [S, secondary]). Put the cost parameters in the settings model within minimum limits.

NIST SP 800-63B-4 was finalized 2025-07-31 [S: csrc.nist.gov]; secondary summaries say 15 characters minimum when the password is the only factor and 8 when it is one of two, no composition rules, a blocklist, no forced rotation [S; normative text not read]. The 12-character policy fits the second-factor case, so the default in the brute-force design (section 4) stands. Suggested: 15 without TOTP, 12 with it; Boss's choice.

### B3. Reset, invite and transfer tokens

OWASP's Forgot Password cheat sheet [S, snippets] asks for the same message whether or not the account exists, uniform response time, secure random tokens, single use, expiry and rate limits. Recommended:

- Token `secrets.token_urlsafe(32)`, stored as SHA-256 only (a fast hash suits a 256-bit random value), like sessions and API keys.
- One table, `account_tokens(id, account_id, purpose ENUM('set_password','reset','email_change','owner_confirm','owner_accept'), token_hash CHAR(64) UNIQUE, created_at, expires_at, used_at, meta JSON)`. A new token for the same account and purpose invalidates older ones.
- Consume with one statement and require `rowcount == 1`, which closes the replay race (design pattern): `UPDATE account_tokens SET used_at=NOW() WHERE token_hash=? AND used_at IS NULL AND expires_at>NOW()`.
- The link opens a GET form; the password is set on POST with the CSRF token; `Cache-Control: no-store`. The token is in the path and reaches access logs, so lifetimes stay short.
- Lifetimes: reset 30 minutes (third-party guidance says 15-30 [S]); first-password link 24 hours; owner-transfer links 24 hours; email change 1 hour.
- A reset ends the account's other sessions and trusted devices (as `change_credentials` does), keeps API keys, and emails a "your password was changed" notice.
- No enumeration: "forgot password" always answers the same text, does the same work (hashing a dummy token) and sends from an RQ job so timing is flat.
- itsdangerous is not used: it cannot know whether a token was used [S], and the flows need use counts, revocation and listings anyway.

### B4. Sessions, cookies and CSRF

- Keep sessions in the control database: the properties of "sessions in Redis" (opaque id, hash at rest, immediate revoke, end-others on credential change) already exist, Redis has no native Windows build, and a flush would sign everyone out.
- Add a fixation test (a made-up cookie that then logs in gets a different cookie) and renew the token after any role or password change. Use the `__Host-` prefix when the cookie is Secure [R].
- Admins keep 12 hours fixed. For ordinary users propose `security.session_hours_user` default 720 (30 days), with an absolute cap and a sliding refresh limited to 7 days, as bounded settings.
- `SameSite=Strict` makes a link from a mail client a cross-site navigation with no session cookie. Reset and invite links need no session; USR.6's confirm link and email-change confirmation do, so they land on a neutral page that asks to sign in (carrying `next`) or shows a "continue" button doing a same-site POST. Test in a browser.
- Keep `web/csrf.py`; Flask-WTF 1.3.0 needs Python 3.10 and would duplicate it. Every new form gets the check because it is app-wide.

### B5. Lockout and throttling for the new flows (USR.1, SEC.1)

Extend the existing design (delays that expire, per-IP and per-username limits, trusted-device cookie, neutral messages; NIST prefers delays [S, secondary]) with `loginguard` namespaces `reset`, `invite` and `register`. Reset: 3 per hour per hashed address and 10 per hour per IP, counted whether or not the address has an account. Invite page: 20 per hour per IP (guessing 256 bits is hopeless, so this is about spam and load). Registration: 5 per hour per IP. Redis keys use a hash of the lowercased address so Redis never stores addresses. Extend `revoke_devices` to resets and role changes. A failed second factor and a failed owner-transfer password prompt count as failed logins (SEC.26 already counts TOTP).

### B6. TOTP and WebAuthn

TOTP is done on pyotp (plus or minus one step, no reuse of a step at or before `last_step`); keep it, require it for the Owner and refuse Owner transfer to an account without it. WebAuthn (py_webauthn; Python 3.10 for 3.0.1, 2.7.1 for 3.9; pulls cryptography, pyOpenSSL, cbor2, pyasn1) is in no USR item: defer it to an optional factor after USR.7, with the relying-party origin from `web.base_url` and a `webauthn_credentials` table (credential id, public key, sign count, transports, label, created, last used). The library's source was not reviewed [S: PyPI release dates only].

### B7. Role model (USR.2)

Answers to USR.2's open questions, with defaults:

- **One table.** Add to `admin_users`: `email VARCHAR(254) NULL` unique, `email_verified_at`, `role ENUM('user','admin','owner') NOT NULL DEFAULT 'admin'` (existing rows are admins), `disabled_at`, `password_changed_at`, `created_by`. Renaming to `accounts` (and `admin_user_id` to `account_id` across six tables) is cleaner but touches every query; add columns now, rename in a follow-up unless Boss wants it in the same PR.
- **Exactly one Owner, enforced by the database:** a generated column `owner_flag TINYINT AS (CASE WHEN role='owner' THEN 1 END) VIRTUAL` with a `UNIQUE` index (NULLs may repeat, one `1` only). Not testable here; verify on MySQL 8.4 and MariaDB 11.4 [R].
- **Which existing account is Owner:** the lowest id (the one `bootstrap_control_schema` seeded). The migration prints its choice and `python -m planetgen.cli.accounts --set-owner USER` (like `planetgen.cli.lockouts`) corrects it.
- **Demotion:** an admin may demote any other admin except the Owner, never themselves (so the last admin cannot vanish); only Owner transfer changes the top role. Every change is audited.
- **Disable and delete:** admins disable and re-enable ordinary users (ends sessions and keys); disabling an admin is Owner-only. Hard delete is by the account itself or the Owner and rewrites `admin_audit_log.admin_username` to `deleted-user-<id>`.
- **API keys:** admins only at first (API.6). Effective rights are the intersection of the key's scope and its owner's current role, evaluated per request (never cached in the key row), so a demoted admin's keys fall to the new role and a user-level key owned by an admin never exceeds read scope. The rule from PR #293 (any API key gets 403 on key management, credentials, 2FA and logout) stays. Scopes come from API.9.
- **Names:** usernames 3-64 characters from `[A-Za-z0-9_.-]`, unique case-insensitively (explicit `utf8mb4_bin` on a lowercased shadow column; collation names differ between MySQL 8.4 and MariaDB 11.4 [R]). Email stored lowercased and unique; `email_verified_at` is set by the first emailed link. Validate with `email.utils.parseaddr` and a length check; the confirmation email is the real validation.

### B8. Invites (USR.4)

Tables `invites(id, token_hash UNIQUE, created_by, created_at, expires_at NULL=never, max_uses NULL=infinite, uses INT, revoked_at, note, invitee_email NULL)` and `invite_uses(invite_id, account_id, used_at)`. Redeem with one atomic statement, `rowcount == 1`, then create the account in the same transaction:

```sql
UPDATE invites SET uses = uses + 1
WHERE token_hash = ? AND revoked_at IS NULL
  AND (expires_at IS NULL OR expires_at > NOW())
  AND (max_uses IS NULL OR uses < max_uses)
```

Boss's expiry list is a fixed enum in the model. Defaults for the open questions: the admin may type the invitee's email, and with SMTP configured the link is sent (and also shown to copy); without SMTP only the link is shown, which answers USR.3's no-SMTP question. Open invites are capped at 50 (a setting), each use of a multi-use link is recorded in `invite_uses`, and new accounts get role `user`.

### B9. SMTP (USR.3)

- Standard-library `smtplib` with `EmailMessage` (Flask-Mail 0.10.0, last release May 2024, wraps about 40 lines), plain text with `Date` and `Message-ID`.
- Modes `starttls` (587, default), `ssl` (465, `SMTP_SSL`), `none` (loopback only); `ssl.create_default_context()` in both TLS modes; class chosen from the mode, not the port; timeout 10 s; log host, port, error class and SMTP code, never the password.
- Send from an RQ job with a small retry, so a slow server holds no thread and timing does not reveal whether an account exists.
- Settings: `smtp.host`, `port` (587), `security`, `username`, `password` (secret), `from_address`, `from_name`, `timeout` (10). "Send test email" goes to the signed-in admin's address and shows the SMTP error text, never the password.
- **Owner only, with recent re-authentication, for SMTP, `base_url` and the From address.** The reset link is built from `base_url` and sent through this server; whoever sets either can read the Owner's reset email and take over the Owner account. This is the one place where admin power must not equal Owner power.
- The password lives in `settings.json` (A5). Deliverability (SPF, DKIM, DMARC) is the SMTP provider's job and the docs say so.

### B10. Password setting, reset and Owner flows (USR.5, USR.6)

Lifetimes are in B3. Changing the email address needs the current password and a confirmation link to the new address, plus a notice to the old one, which stays until the new one confirms. An Owner reset works by email, but with TOTP enrolled the reset page also asks for the code, and each one emails the Owner and writes an `AUTH` event. Break-glass for an Owner without password and email is the console: `python -m planetgen.cli.accounts --reset-password OWNER` and `--set-owner` (API.15's "console user is god").

**Owner transfer** is a state machine in `owner_transfers(id, from_id, to_id, state, confirm_token_hash, accept_token_hash, created_at, expires_at, ...)` with states `pending_owner`, `pending_target`, `done`, `cancelled`, `expired`, and at most one active transfer:

1. The Owner picks an existing admin, retypes the password twice, enters the TOTP code if enrolled, ticks the explicit confirmation; an email goes to the Owner.
2. The Owner clicks it; an email goes to the target.
3. The target signs in and accepts (the click alone is not enough: signing in proves password and TOTP).
4. One transaction locks both rows (`SELECT ... FOR UPDATE`), sets old to `admin` and new to `owner` (the unique index forbids two Owners), ends both accounts' sessions and trusted devices and writes the audit rows. The old Owner gets a final notice.

Defaults: the new Owner must already be an admin with TOTP; the transfer expires in 48 hours; the Owner may cancel before acceptance; unconfigured SMTP disables it.

### B11. Audit, data minimization and USR.7

Add to `admin_audit_log` and the activity log: `account.create|disable|enable|delete`, `role.change`, `invite.create|revoke|use`, `password.reset.request|complete`, `password.set`, `email.change.request|confirm`, `owner.transfer.start|confirm|accept|cancel|expire|done`, `config.update`, `smtp.test`. Log a token-hash prefix, never tokens or email bodies; role and owner events are not pruned with failed logins.

The account record is username, email, password hash, role and timestamps; no IP addresses, no tracking pixels, no third-party calls. The account page offers "download my data" JSON and "delete my account" (bookmarks, keys and sessions cascade), and `seo.privacy_note` gives the operator a place for a privacy note. GDPR specifics were not researched [R].

USR.7: bookmarks are keyed by account and a stable object reference (NAV.7's references); decision 4 of [galaxy-drilldown-navigation.md](galaxy-drilldown-navigation.md) is a prerequisite per the TODO. A broken bookmark after a regenerate shows the stored name and "no longer in this galaxy". If reading ever requires sign-in, `seo.indexing=off` must be forced; the default stays public read.

### B12. USR.8: the one-off system limit

Today `/admin/generate/system` and `/download` in `web/system_page.py` are admin-only (`_admin_or_redirect`, `_admin_or_403`). Flask-Limiter semantics (limits 4.2 under Flask-Limiter 3.11.0) [C]:

- The bucket is whatever `key_func` returns: `user:<account id>` when signed in, `ip:<address>` otherwise. Use `"30 per hour"` with `strategy="moving-window"`, which counts the last 3,600 seconds exactly (one entry per hit; trivial at 30). `fixed-window` can pass nearly double the limit across a boundary. `sliding-window-counter` approximates and, after 30 hits at t=0, lets 1 through at 60m01s and 29 at 119m.
- The 31st request is a 429; `Retry-After` and `X-RateLimit-*` appear because `api/config.py` sets `RATELIMIT_HEADERS_ENABLED`, and the page shows the wait in words. `exempt_when` takes a callable; admins and the Owner pass it (40 admin requests all 200).
- In-process API calls skip the default limits but not explicit `@limiter.limit(...)` ones (`api/limiter.py`), so decorate the web view directly, as `page_limit` does.

Count in Redis, not the control database as USR.8's default says: the data is ephemeral, SEC.30 already put the limiter on Redis, and a database row per generation is pure write load. A Redis outage falls back to per-process counting (`RATELIMIT_IN_MEMORY_FALLBACK_ENABLED`), fine for a fairness limit. Add a concurrency cap, since each request starts a generator process for about a second on a 5-thread server and a per-hour limit does not stop 30 requests in one second: for example 3 at once, a Redis counter with a TTL, "busy, try again". Settings: `limits.one_off_per_hour` (30), `limits.one_off_concurrent` (3).

## Corrections to existing statements

- [login-brute-force-protection.md](login-brute-force-protection.md) section 5 said "Planned: libraries" for work that is done (SEC.29, SEC.30); corrected in place.
- The research draft said `loginguard.py` uses Flask-Limiter's storage and named control schema v9; the lockout counters are the project's own `RedisStore` on the same Redis server, and the control schema is at v10.
- `docs/config.md` makes `config.json` read-only to the web process, so ADM.43's web-saved configuration is impossible as written (A2); its "`base_url` not yet read" changes with ADM.44.
- NIST SP 800-63B-4 asks for 15 characters when the password is the only factor; the 12-character policy holds only while a second factor is expected (B2).

## Evidence notes

[C] results used Python 3.9.25, pydantic 2.13.5, pydantic-settings 2.11.0, Flask-Limiter 3.11.0, limits 4.2, argon2-cffi 25.1.0, Werkzeug 3.1.9 and Pillow 11.3.0, on a slow, noisy VM; the scripts were scratch files and are not in the repository. [R] items to verify when paper and vendor access is allowed:

- Google: a long-lived `noindex` becoming `nofollow`; robots.txt cached 24 hours and read to 500 KiB (translated snippets only); the 155-character snippet length.
- Discord: use of `theme-color`, whether Discordbot obeys robots.txt, whether `twitter:card` is needed for a large image; X's card ratio (sources said 1.91:1 and 2:1).
- Inline `application/ld+json` under the CSP; RFC 9116 requiring `Contact` and `Expires`.
- The unique generated-column Owner index and username collations on MySQL 8.4 and MariaDB 11.4; `__Host-` cookie rules; gunicorn reload by SIGHUP; that table A4 matches every call site.
- NIST SP 800-63B-4's normative text, the OWASP wording on hashing reset tokens, the 15-30 minute reset lifetime; py_webauthn maintenance; the GDPR points in B11.

## Sources

Package pages (PyPI JSON, fetched 2026-10-09): pydantic, pydantic-settings, pyotp, argon2-cffi, webauthn, flask-limiter, limits, flask-wtf, flask-mail, itsdangerous, werkzeug, pillow (`https://pypi.org/pypi/<package>/json`).

Search results (secondary unless noted):

- OWASP, snippets only: https://cheatsheetseries.owasp.org/cheatsheets/Password_Storage_Cheat_Sheet.html, https://cheatsheetseries.owasp.org/cheatsheets/Forgot_Password_Cheat_Sheet.html
- NIST SP 800-63B-4: https://csrc.nist.gov/pubs/sp/800/63/B/4/final (listing), https://www.sakimura.org/en/2025/10/7710/, https://www.enzoic.com/blog/nist-sp-800-63b-rev4/
- Sitemaps: https://library.linkbot.com/what-are-the-url-and-file-size-limits-for-sitemaps-and-how-can-large-sites-adapt/, https://blogs.bing.com/webmaster/november-2016/increasing-the-size-limit-of-sitemaps-file-to-addr, https://developers.google.com/search/blog/2023/06/sitemaps-lastmod-ping (title only), https://dynomapper.com/blog/search-engine-optimization/xml-sitemaps-seo-best-practices/
- Search engines: https://www.searchenginejournal.com/google-on-robots-txt-when-to-use-noindex-vs-disallow, https://www.semetrical.com/insights/noindex-vs-disallow-why-robotstxt-wont-stop-low-quality-pages-from-getting-indexed/, https://www.searchenginejournal.com/google-stopped-supporting-relprev-next-in-search-indexing-years-ago/299689/, https://developers.google.com/crawling/docs/robots-txt/robots-txt-spec (translated pages only), https://webmasters.googleblog.com/2009/09/google-does-not-use-keywords-meta-tag.html, https://ppc.land/google-to-remove-sitelinks-search-box/, https://www.mediapost.com/publications/article/400424/
- mod_wsgi: https://modwsgi.readthedocs.io/en/latest/user-guides/reloading-source-code.html
- Link previews: https://support.discord.com/hc/en-us/articles/42500550752919-About-Discord-Link-Previews-and-the-Discordbot (title only), https://htmlcsstoimage.com/features/og-images-for-discord, https://env.dev/guides/opengraph-image-not-showing, https://screenhance.com/blog/og-image-size-guide, https://quickseo.ai/tools/twitter-card-validator, https://superblog.ai/blog/ai-crawlers-guide, https://www.texta.ai/blog/ai-crawler-user-agents
- itsdangerous: https://itsdangerous.palletsprojects.com/url_safe/

Repository files read: `docs/TODO.md` (ADM.42-44, ADM.18, USR.1-8, API.6, API.15, UX.49, VIEW.3), `docs/config.md`, `config.json.example`, the code under `src/planetgen/` named above, `db/control_schema.sql`, `setup.py` and `requirements.lock`.

See also [api-design-standards.md](api-design-standards.md) for the API conventions (keys, errors, paging, logging) that API.6 and the account items rely on.
