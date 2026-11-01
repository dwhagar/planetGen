# Deployment Configuration

This document describes `config.json`, the single per-deployment
configuration file for the whole planetGen project -- the command-line
tools (`planetgen`, `planetgen.db.query`, `planetgen.cli.migrate`, `planetgen.cli.orbits`,
`planetgen.cli.reset`, ...) and the Flask app (the API in
`../src/planetgen/api/` and the pages in `../src/planetgen/web/`) all read it -- and
[`../config.json.example`](../../config.json.example), the committed
template it's copied from.

## What it's for

Before this file existed, every one of those entry points read its own
handful of `PLANETGEN_*` environment variables (`planetgen.db.store.MySQLConfig`,
`../src/planetgen/api/config.py`, `../src/planetgen/web/lib/apiclient.py`),
each with its own hardcoded default. That works, but it means a
deployment that just wants "one MySQL server, one account, one API base
URL" still has to either edit code or set half a dozen `SetEnv`/
`EnvironmentFile` lines to get there. `config.json` is a single JSON file,
edited once per deployment, holding every one of those settings in one
place.

It is loaded by [`../src/planetgen/util/settings.py`](../../src/planetgen/util/settings.py):
one Pydantic model, `Settings`, holds every option with its type, default,
help text, unit, whether it is a secret, whether a change needs a restart
and which environment variable overrides it. A `config.json` that only sets
the fields it needs to change is enough; anything it omits keeps its
default value, a value of the wrong type stops the program with a message
naming the field, and `python3 -m planetgen.cli.config check` validates the
file without starting anything (an option name it doesn't know is an
error there, so a typo such as `ratelimt` is caught). The debug switch and
the two log locations are read by the standard-library-only
[`../src/planetgen/util/logpaths.py`](../../src/planetgen/util/logpaths.py) instead, because the
installers run it as root before any virtual environment exists.

## Location: repo root

`config.json` lives at the repo root -- next to `src/` -- rather than
inside `../src/html/` itself (its *example*
template, [`../config.json.example`](../../config.json.example), lives at the
repo root too, since a template holds only placeholder values, not
secrets). This mirrors how Apache's `DocumentRoot` for this application is
`../src/html/` alone (see
[`../examples/apache/planetgen.conf.example`](../../examples/apache/planetgen.conf.example)),
and the other web servers only serve `src/html/static/` from disk (see
[`deployment/`](../deployment/README.md)), so `config.json`, outside those
trees, can never be requested over HTTP no matter how the site's rules
are written. This matters because
`config.json` holds real MySQL credentials, unlike the placeholder fields
its now-removed predecessor (`webconfig.json`) only ever reserved for
them.

It also has to be unreadable by other local users. `install.sh` and
`update.sh` (through `examples/apache/set-permissions.sh`) set it to
`<apache user>:<apache group>`, mode `640`, whenever it exists (the whole
checkout belongs to that account): the web app (`www-data` for Apache and
the gunicorn unit alike) can read it, root can too, nobody else can. The macOS
guide sets the same rule by hand. A login user who runs the
generator or the maintenance scripts from a shell without `sudo` needs to
be in Apache's group (`sudo usermod -aG www-data <user>`, then log in
again) to read it. After creating or copying a new `config.json`, re-run
`sudo examples/apache/set-permissions.sh` (or `sudo ./update.sh`), or set
it by hand: `sudo chown www-data:www-data config.json && sudo chmod 640 config.json`.

## Precedence

Every setting below has up to four sources, checked in this order:

1. An explicit function/CLI argument (e.g. `planetgen sector --mysql-password ...`,
   or a test building its own `MySQLConfig(...)` directly) -- always wins.
2. The matching `PLANETGEN_*` environment variable, if set and not empty
   (a MySQL option set to an empty value counts).
3. `config.json`, the one file that holds all program configuration. It is your own
   copy of `config.json.example`: the example is only a template that documents every
   option and its default, and nothing ever writes to it or overwrites your copy
   (`update.sh` leaves `config.json` alone, and an option added in a later release has
   its default in the code, so an older copy keeps working). The Admin settings page
   (ADM.43) and `python3 -m planetgen.cli.config set` edit `config.json` itself; see
   [Saving settings](#saving-settings). If the file does not exist they say so and
   create nothing.
4. The built-in default (matching `config.json.example`'s values).

Environment variables are still checked *before* `config.json` rather than
being replaced by it, so a single shared `config.json` can still be
overridden per-process where that's needed -- most notably
`planetgen-orbits@.service`'s per-instance `PLANETGEN_MYSQL_DATABASE=%i`
(see [`../examples/maintenance/`](../../examples/maintenance/)), which relies on
one templated systemd unit overriding just the database name per timer
instance while everything else comes from the shared credentials file (or
now, `config.json`).

## Saving settings

The writer (`planetgen.util.settings_writer`) is the only thing that writes
`config.json`; the command line and the Admin settings page (`/admin/settings`,
below) both use it. The file must belong to the web server's user (`update.sh`
and `set-permissions.sh` see to that) so the page can replace it; when it doesn't,
the page lists the options and says so instead of offering a Save button.

- **Who may change what.** The web (an admin) changes the options the table marks
  "any admin"; the Owner (from 8.2) the "the Owner" ones; only the command line,
  run on the server, changes "command line only" options: the database, the
  secrets that open the host (`secret_key`), `jobs.python`, the log paths. A secret
  can be replaced but never read back: pages, diffs, audit rows and logs show
  `(set)`, `(empty)` or `(changed)`.
- **Checked first.** An option that doesn't exist, that the caller's role may not
  change, or whose value doesn't fit its type or range is refused, each at its own
  option name, and nothing at all is saved. What is stored is the value the model
  parsed (`"30"` is stored as `30`). An option left equal to its default is not
  stored.
- **Everything else stays.** The file is edited as parsed JSON: keys the model
  doesn't know and the order of the keys are kept; only the changed options move.
- **Written atomically.** A temporary file in the same folder, `fsync`, the old
  file's mode, a rename, then `fsync` of the folder: a reader sees the old file or
  the new one, never half of one. A lock file next to it keeps two saves (two
  admins, or the page and the command line) from interleaving.
- **Backed up.** Just before a save, restore or import replaces the file, a
  settings-only backup archive is taken ([backups](#backups) below): the old file,
  who is making the change and when. The newest `settings.backups_kept`
  (default 20, 1 to 200) archives are kept. A restore from the web puts back only
  the options an admin may change; the command line puts the whole file back.
- **Audited.** One `admin_audit_log` row per changed option (`config.update`,
  `config.restore`; target `config:<option>`; detail `{"old":..,"new":..}`) and
  one `DB` line in the activity log. A secret's values are written `(changed)`.
- **Restart.** The result lists the changed options marked "after a restart".

## The settings page

`/admin/settings` (Admin, Settings tab; ADM.84) lists every option of the model
from the model itself, so an option added with `settings.opt(...)` appears there
with its help, default, limits and unit and needs no page code. Each shows where
the value in use comes from: the default, `config.json`, or an environment
variable (which still wins over the file; the page and the review say so). Change
what you want and press Review the changes: a bad value is sent back with its
message at its option and nothing is saved; a good form shows the old and new
value of each changed option, marks the ones that need a restart, and Save these
changes writes it as described in [Saving settings](#saving-settings). Use the
default removes a key from the file; a secret's box is write-only (empty keeps
it, typing replaces it, Clear it removes it). While a saved option still needs a
restart of the web server, every admin page shows a banner naming it; the banner
goes when the server restarts. The page is described in [html-interface.md](html-interface.md).

## Backups

`planetgen.backup` (ADM.87) takes every backup as one `backups/<UTC time>.tar.bz2` in
the state folder (`/var/lib/planetGen/state`, `PLANETGEN_STATE_DIR` moves it), mode
0600, owned by the web user (`set-permissions.sh` keeps the folder 700). It holds
`manifest.json` (format 1, time, who, why, the version and, per part, its folder,
sizes and SHA-256) and one folder per part: `settings/` (a copy of `config.json`,
secrets included, which is why the archive is private and its download is for the
Owner) now, and `control/` (the control database) and `parity/` (the DB repair
parity data) when those parts arrive in 9.0. A save on the settings page takes an
archive of the settings part alone; `python3 -m planetgen.cli.backup take` takes
every part.

Retention is the option `settings.backups_kept`: the newest archives are kept, the
rest deleted, but never the newest archive that holds every part. A restore verifies
the archive first, backs up what it is about to replace in a safety archive, then
restores the whole archive or only the parts named; a part this build has no
provider for is skipped and reported, and the archive keeps it. The galaxy database is
not in an archive: it is regenerated from its seed or repaired from parity.

### The backups screen, export and import

`/admin/settings/backups` (ADM.85; a link on the settings page) is the one screen for
archives. It lists them newest first (when, who, why, which parts and how big), and
**Back up now** makes one of every part. **Compare and restore** opens an archive: it
is checked against its hashes, the settings part is compared with today's
`config.json` (every option a restore would change, a secret as `(changed)`, a backup
that no longer fits the model refused with the reason), and the parts ticked are
restored: the whole archive, or one part. The restore makes a safety archive of what
it replaces first, writes a `config.restore` audit row per changed option and a
`backup.restore` row, and from the web puts back only the options an admin may change.
**Download** is the Owner's, from 8.2 (an archive holds every secret); until then
archives are read on the server.

**Export** downloads the options set in `config.json` that an admin may change, as
JSON in the file's own layout with `config_version`: never a secret and never a
command-line-only option. **Import** takes such a file and checks it exactly as an
edit is checked (each option's role, type and range; unknown keys and keys an admin
may not set are errors; a file from a newer release is refused), shows the review
diff, and saves it as `config.import` with a backup first. It sets the options the
file names and leaves every other option as it is; one bad value changes nothing.

### Web, link previews and search engines

The `site`, `web`, `icons`, `seo` and `og` sections (ADM.44) set what the pages tell
browsers, search engines and chat programs. All of them apply to the next page view
with no restart; the settings page shows them under Site, Web pages, Icons, Search
engines and Link previews. They are all read at request time by `planetgen.web.seo`.

- **The address.** Canonical links, `og:url`, `og:image`, the sitemap and the robots.txt
  `Sitemap:` line are built from `base_url` (the Owner changes it), never from the
  request's Host header. While it is still `http://localhost/` none of those absolute
  addresses is written and `/sitemap.xml` is a 404.
- **The title and tab.** `web.title_template` (`{title} - {site_name}`), `site.language`
  (the `lang` of every page and the Open Graph locale), `site.tagline` (under the home
  page heading, and the description of pages that give none), the theme colours.
- **What may be indexed.** `seo.indexing` is `off` (robots.txt disallows `/`, every page
  says `noindex`), `top_pages` (the home page, the class, species, polity and galaxy
  pages and the lists) or `all` (also the system, sector and phenomenon pages);
  `seo.detail_pages = index` does the same for the generated pages alone. Millions of
  generated pages look like scaled content to a search engine, so they stay out by
  default. The data and tool routes (search, navigation, the Galaxy Map's nodes, the
  `/scene` and `/table` data) answer `X-Robots-Tag: noindex` and are disallowed in
  robots.txt; the admin, login and account pages say `noindex,nofollow`. A page that has
  a query (`?page=2`, a filter) is not indexed but its links are followed, and its
  canonical link drops tracking and sort parameters. A page the code has not classified
  is treated as a tool route, and a test fails until it is classified on purpose
  (`planetgen.web.seo`).
- **robots.txt.** Generated from the routes: every tool and private route is disallowed
  with `*` wildcards, `block_ai_crawlers` adds the AI crawlers, `robots_extra` is
  appended and `robots_override` replaces the whole file. A line that is not a
  directive is refused where it is typed; an override that disallows `/` while
  `indexing` is not `off` is allowed but warned about, since it hides the site.
- **Sitemaps.** `/sitemap.xml` is an index of `/sitemap-pages.xml` and, when the
  generated pages may be indexed, `/sitemap-systems-<n>.xml` and
  `/sitemap-sectors-<n>.xml` of at most `seo.sitemap_max_urls` addresses each, built
  from the database 500 at a time and kept for `seo.sitemap_cache_hours`. Phenomena are
  not listed. The number of files comes from a count that may be an estimate, so the
  last may be empty.
- **Other files.** `/.well-known/security.txt` and `/humans.txt` serve the text in
  `seo.security_txt` and `seo.humans_txt` (404 while empty); `/site.webmanifest` is served
  while `icons.manifest_enabled`; `seo.privacy_note` is shown on the account page.
- **Link previews.** `og.enabled` adds Open Graph and Twitter tags to the public pages
  (never the admin, login or account ones). The picture is `og.default_image`, else the
  built-in 1200 by 630 one. `og.generated_cards` is for a per-system and per-sector card:
  with it on, `/og/<system|sector>/<id>.png` draws the card with a renderer installed
  by `planetgen.web.seo_pages.register_card_renderer` (the drawer is VIEW.3; with none,
  every page uses the default picture), caches it in `og.card_cache_dir` (the state
  folder's `og-cards` by default) up to `og.card_cache_mb`, and falls back to the default
  picture when it can't draw.
- **Icons and pictures.** The icon options hold the name of an uploaded file, not a
  path: upload on the Images page (`/admin/settings/images`, linked from the settings
  page). A picture is checked on upload (PNG, JPEG, ICO or a sandboxed SVG; size and
  squareness where they matter; at most 1 MiB), stored as `<16 hex digits of its
  SHA-256>.<ext>` in `web.assets_dir` (the state folder's `assets` by default) and served
  at `/site-assets/<name>` for a year (its name is its content). Uploading and removing
  save the setting through the same writer, so they are backed up and audited. These site pictures take no GIF; an account's own picture (docs/reference/api.md, `POST /api/auth/avatar`) is the one place a GIF is accepted.

## Options

The table is generated from the settings model (`src/planetgen/util/settings.py`, ADM.42); change an option there and run `python3 -m planetgen.cli.config docs`. Next to each option: its type, default, whether a change applies at once or after a restart, who may change it from the web ("command line only" means on the server only) and the environment variable that overrides it.

<!-- BEGIN GENERATED OPTIONS (python -m planetgen.cli.config docs) -->

| Option | Setting | Type | Default | Applies | Edited | Environment | What it does |
|---|---|---|---|---|---|---|---|
| `config_version` | Settings file version | whole number (format version) | 1 | at once | command line only | none | The settings file format version; a file without one is version 1. |
| `site_name` | Site name | text | planetGen | at once | any admin | none | Display name for this deployment, shown in the header, the home page heading, every page title and the footer. |
| `base_url` | Site address | text | http://localhost/ | at once | the Owner | none | The address this site is reached at, as an absolute http or https address with a trailing slash. Canonical links, link previews, the sitemap and emailed links are built from it, never from the request's Host header. |
| `api_base_url` | API address | text | empty | after a restart | command line only | `PLANETGEN_API_BASE_URL` | Base URL of the Flask API's /api mount point (for example https://example.org/api). The update asks it whether the site answers, and web/lib/apiclient.py calls it outside the Flask app. There is no default: nothing guesses a loopback address, so it must be set for either to work. |
| `debug` | Debug log | true or false | false | after a restart | command line only | `PLANETGEN_DEBUG` | Debug mode: every decision the generator makes, every random roll, every SQL statement, queue event and request, written to the console of each planetgen command (stdout) and to log_file. The output is large. |
| `debug_areas` | Debug areas | text | empty | after a restart | command line only | `PLANETGEN_DEBUG_AREAS` | Limit debug mode to some areas of output: decision, draw, sql, queue, cache, web, lib (comma separated; a leading - leaves one out). Empty writes everything but branch, which logs every conditional jump of the program's own code and is written only when named. |
| `log_file` | Debug log file | text | /var/log/planetgen.log | after a restart | command line only | `PLANETGEN_LOG_FILE` | Where the debug log goes. |
| `record_file` | Record file | text | empty | after a restart | command line only | `PLANETGEN_RECORD_FILE` | A file every log line of a planetgen command is also written to as a record: one per line, the time in UTC. Empty writes none; the --record-file option does the same for one command. |
| `record_format` | Record file form | "logfmt" or "json" | logfmt | after a restart | command line only | `PLANETGEN_RECORD_FORMAT` | The form of the record file's lines: logfmt (key=value) or json (one object per line). |
| `log_dir` | Activity log folder | text | empty | after a restart | command line only | `PLANETGEN_LOG_DIR` | The folder of the always-on activity log; empty means /var/log/planetgen on Linux, /Library/Logs/planetgen on macOS. |
| `log_rotation` | Log rotation | "auto" or "system" or "app" | auto | after a restart | command line only | none | How the activity log is rotated: system (logrotate or newsyslog), app (the program, 100 MB per file, 30 old copies) or auto (system where an installed rotation config exists, else app). |
| `settings.backups_kept` | Settings backups kept | whole number (backups) | 20 | at once | any admin | none | How many backup archives are kept. One is made just before every settings save, restore or import replaces the file, and by Back up now; the oldest is deleted once there are more than this, but never the newest archive that holds every part. |
| `site.tagline` | Tagline | text | empty | at once | any admin | none | A short line shown under the home page's heading and used as the description of pages that give none of their own. |
| `site.language` | Language | text | en | at once | any admin | none | The pages' language as a code (en, de, en-GB): the `lang` of every page and the Open Graph locale. |
| `web.title_template` | Page title format | text | {title} - {site_name} | at once | any admin | none | How a page's title is written in the browser tab and in search results; {title} is the page's own title and {site_name} the site name. |
| `web.theme_color_light` | Browser colour, light theme | text | #f6f7fb | at once | any admin | none | The browser's toolbar colour on phones for the light theme, as #rrggbb. |
| `web.theme_color_dark` | Browser colour, dark theme | text | #14151e | at once | any admin | none | The same for the dark theme. |
| `web.assets_dir` | Uploaded images folder | text | empty | at once | command line only | `PLANETGEN_ASSETS_DIR` | Where uploaded icons and images are kept; empty means the assets folder in the state folder. |
| `icons.favicon_svg` | Site icon (SVG) | text | empty | at once | any admin | none | The site icon as an uploaded SVG; empty uses the built-in one. Set it by uploading on the Images page. |
| `icons.favicon_png` | Site icon (PNG or ICO) | text | empty | at once | any admin | none | A PNG or ICO icon for browsers that do not take an SVG; empty offers none. Set it by uploading on the Images page. |
| `icons.apple_touch_icon` | Home-screen icon | text | empty | at once | any admin | none | The 180 by 180 PNG iPhones and iPads use for a home-screen bookmark. Set it by uploading on the Images page. |
| `icons.manifest_enabled` | Installable site | true or false | false | at once | any admin | none | Serve /site.webmanifest, which names the site, its icon and its theme colour for browsers that can install it. |
| `seo.indexing` | Search engine indexing | "off" or "top_pages" or "all" | top_pages | at once | any admin | none | Whether search engines may index the site: off (robots.txt disallows everything and every page says noindex), top_pages (the home, class, species, polity, galaxy and list pages), or all (also the system, sector and phenomenon pages). |
| `seo.detail_pages` | Index generated pages | "noindex" or "index" | noindex | at once | any admin | none | Whether the generated system, sector and phenomenon pages may be indexed when indexing is top_pages. Millions of generated pages look like scaled content to search engines, so they stay out by default; a small galaxy can switch this to index. |
| `seo.description` | Search description | text | empty | at once | any admin | none | The description search engines show for pages that give none of their own (about 155 characters are shown); empty uses the tagline. |
| `seo.keywords` | Keywords | list of text | empty list | at once | any admin | none | Keywords for the keywords meta tag. Google has ignored it since 2009; some smaller engines still read it. |
| `seo.json_ld` | Structured data | true or false | true | at once | any admin | none | Add a WebSite block (schema.org) to the home page. |
| `seo.robots_extra` | Extra robots.txt lines | text | empty | at once | any admin | none | Lines added to the end of the generated robots.txt. |
| `seo.robots_override` | Replace robots.txt | text | empty | at once | any admin | none | Replaces the generated robots.txt entirely when not empty. A mistake here can hide the whole site from search engines. |
| `seo.block_ai_crawlers` | Block AI crawlers | true or false | false | at once | any admin | none | Disallow the crawlers that collect text to train AI models (GPTBot, ClaudeBot, CCBot, Google-Extended, Applebot-Extended) in robots.txt. |
| `seo.sitemap_enabled` | Sitemap | true or false | true | at once | any admin | none | Serve /sitemap.xml and list it in robots.txt. |
| `seo.sitemap_max_urls` | Sitemap size | whole number (addresses) | 50000 | at once | any admin | none | Most addresses in one sitemap file (the limit is 50,000). |
| `seo.sitemap_cache_hours` | Sitemap refresh | whole number (hours) | 24 | at once | any admin | none | How long a built sitemap is reused. |
| `seo.security_txt` | security.txt | text | empty | at once | any admin | none | The text of /.well-known/security.txt (RFC 9116 needs a Contact and an Expires line); empty serves none. |
| `seo.humans_txt` | humans.txt | text | empty | at once | any admin | none | The text of /humans.txt; empty serves none. |
| `seo.privacy_note` | Privacy note | text | empty | at once | any admin | none | A short note shown on the account page about what the site keeps. |
| `og.enabled` | Link previews | true or false | true | at once | any admin | none | Add Open Graph and Twitter tags to the public pages, so a link shows a title, a description and a picture when shared. |
| `og.site_name` | Name in previews | text | empty | at once | any admin | none | The site name in link previews; empty uses the site name. |
| `og.image_alt` | Preview picture description | text | empty | at once | any admin | none | A description of the default preview picture for people who cannot see it. |
| `og.default_image` | Preview picture | text | empty | at once | any admin | none | The picture link previews show (a PNG or JPEG about 1200 by 630); empty uses the built-in one. Set it by uploading on the Images page. |
| `og.twitter_card` | X (Twitter) card style | "summary" or "summary_large_image" | summary_large_image | at once | any admin | none | How X (Twitter) draws a preview: a large picture or a small one. |
| `og.twitter_site` | X (Twitter) account | text | empty | at once | any admin | none | The site's X (Twitter) account, as @name; empty names none. |
| `og.generated_cards` | Preview picture per page | true or false | false | at once | any admin | none | Draw a preview picture per system and sector instead of one for the whole site. Off until a card drawer is installed; every page then uses the default picture. |
| `og.card_cache_mb` | Preview picture disk use | number (MB) | 100 | at once | any admin | none | How much disk the drawn preview pictures may use. |
| `og.card_cache_dir` | Preview picture folder | text | empty | at once | command line only | none | Where drawn preview pictures are cached; empty means the state folder. |
| `mysql.host` | Database host | text | 127.0.0.1 | after a restart | command line only | `PLANETGEN_MYSQL_HOST` | The MySQL server's host name or address. |
| `mysql.port` | Database port | whole number (port number) | 3306 | after a restart | command line only | `PLANETGEN_MYSQL_PORT` | The server's port. |
| `mysql.user` | Database user | text | planetgen | after a restart | command line only | `PLANETGEN_MYSQL_USER` | The account every entry point connects with; give it whatever grants the most demanding caller needs. |
| `mysql.password` | Database password | text | empty | after a restart | command line only | `PLANETGEN_MYSQL_PASSWORD` | That account's password. **Secret.** |
| `mysql.database` | Galaxy database | text | planetgen | after a restart | command line only | `PLANETGEN_MYSQL_DATABASE` | The galaxy database the pages show and the generators fill. |
| `mysql.database_prefix` | Galaxy database prefix | text | planetgen | after a restart | command line only | `PLANETGEN_MYSQL_DATABASE_PREFIX` | The schema-name prefix the database picker lists, for a deployment with several galaxies on one server. |
| `mysql.statement_timeout_seconds` | Longest page query | number (seconds) | 10 | at once | any admin | none | The longest any statement on the web interface's and API's read-only connections may run (0 turns the limit off). A query it stops gives a 504 page, or a QUERY_TIMEOUT error from the API. |
| `control_database` | Control database | text | planetgen_control | after a restart | command line only | `PLANETGEN_CONTROL_DATABASE` | The MySQL schema holding admin identities, sessions, API keys and the audit log. Never listed or selectable as a galaxy. |
| `redis.url` | Redis address | text | redis://127.0.0.1:6379/0 | after a restart | command line only | `PLANETGEN_REDIS_URL` | The Redis server the work queue, the rate limits and the login lockouts use. |
| `queue.interactive_workers` | Workers for short jobs | whole number (workers) | 1 | at once | any admin | none | How many workers are kept for short jobs: edits, estimates and slow page answers. They are never counted in the bulk worker count, so a generation run on every core leaves them be, and they run at a higher priority than the bulk workers. |
| `queue.maintenance_workers` | Workers for upkeep during a run | whole number (workers) | 1 | at once | any admin | none | How many workers run the cache warm-up, the node prefetch and the database check while a generation run is live. They are taken out of the generation worker count (a generation worker finishes its unit and leaves when they are needed), so the host's total does not grow; with no generation run live they may use the whole count. |
| `queue.max_host_workers` | Most workers on this host | whole number (workers) | 0 | at once | any admin | none | The most worker processes planetgen may run on this host at once, every queue and command together: a start that would pass it is refused and said in the log. 0 sizes it from the host: 60 percent of its memory in workers of 300 MB each, no more than four per core, and never fewer than the cores plus four. It stops a burst of requests, or several commands starting workers together, from filling the memory and getting MySQL killed. |
| `queue.long_workers` | Workers for long jobs | whole number (workers) | 1 | at once | any admin | none | How many workers run the jobs expected to take longer than the long job limit (a neighbourhood or a large settle, which can run for hours). They have a lane of their own so a long job can never hold the workers kept for short ones. |
| `queue.long_job_seconds` | Long job limit | number (seconds) | 60 | at once | any admin | none | A job expected to take longer than this goes to the long lane. A job with no estimate counts as short. |
| `queue.unit_ceiling_seconds` | Unit length ceiling | number (seconds) | 300 | at once | any admin | none | The longest a unit of a generation stage is expected to run. A stage is cut into units no longer than this at the default galaxy size; a unit that runs longer is written to the run's log, naming the stage. A stage that is one unit by nature is not held to it. A unit this long is also the most a worker waits to hand its place to a lane that needs it. |
| `queue.wait_minutes` | Longest wait for a run | number (minutes) | 60 | at once | any admin | none | How long a generation run, or an API job that generates sectors, waits in line when another generation run holds the machine, before it fails naming that run (0 does not wait). |
| `queue.max_waiting` | Waiting list length | whole number (entries) | 8 | at once | any admin | none | How many entries the one waiting list holds: generation runs and API jobs that wait for the machine, and Generate page jobs that wait for the running job. One more than this is refused, naming the list. |
| `ratelimit.default` | API request limit | text | 200 per day;50 per hour | after a restart | any admin | `PLANETGEN_RATELIMIT_DEFAULT` | Flask-Limiter's default limit for the API. |
| `ratelimit.storage_uri` | Where limits are counted | text | empty | after a restart | any admin | `PLANETGEN_RATELIMIT_STORAGE_URI` | Where limits count; empty means the Redis server in redis.url, memory:// counts per process (right with exactly one worker). |
| `ratelimit.pages.search` | Search limit | text | 30 per minute | after a restart | any admin | none | Limit per client address on /search and /galaxy/locate. Empty turns it off. |
| `ratelimit.pages.galaxy` | Galaxy Map page limit | text | 60 per minute | after a restart | any admin | none | Limit on /galaxy, the Galaxy Map page. Empty turns it off. |
| `ratelimit.pages.galaxy_tiles` | Map data limit | text | 600 per minute | after a restart | any admin | none | Limit on /galaxy/nodes, /galaxy/stage, /galaxy/window, /galaxy/near and /galaxy/territories, fetched as the map's camera moves. Empty turns it off. |
| `ratelimit.pages.health` | Health check limit | text | 60 per minute | after a restart | any admin | none | Limit on /api/health. Empty turns it off. |
| `ratelimit.pages.other` | Other pages limit | text | 300 per minute | after a restart | any admin | none | Limit on every other page, counted together. Empty turns it off. |
| `api_log.internal` | Log the pages' own API calls | "all" or "off" | all | at once | any admin | none | Whether the calls the web pages make to the API in-process (about six for each page view) are written to the activity log as API lines: all, or off. |
| `login_allowlist` | Never-locked-out addresses | list of text | empty list | after a restart | the Owner | `PLANETGEN_LOGIN_ALLOWLIST` | Client addresses or networks (CIDR) never locked out after failed logins; loopback never is either. Environment: comma- or space-separated. |
| `admin_cookie_insecure` | Admin login over plain HTTP | true or false | false | after a restart | command line only | `PLANETGEN_ADMIN_COOKIE_INSECURE` | Send the admin session cookie over plain HTTP. Only for local development without TLS; a production deployment must never set it. |
| `secret_key` | Form signing key | text | empty | after a restart | command line only | `PLANETGEN_SECRET_KEY` | Signs the CSRF tokens on the pages' forms. Empty makes one and keeps it in the state folder (mode 0600), so every worker and restart shares it. **Secret.** |
| `password_hash.memory_kib` | Password hashing memory | whole number (KiB) | 19456 | after a restart | command line only | `PLANETGEN_PASSWORD_HASH_MEMORY_KIB` | Memory one admin password check uses, in KiB (Argon2id). 19456 (19 MiB) with time_cost 2 and parallelism 1 is OWASP's minimum; raise it where the host has the memory to spare. A stored hash made with other numbers is re-hashed at that admin's next sign-in. |
| `password_hash.time_cost` | Password hashing passes | whole number (passes) | 2 | after a restart | command line only | `PLANETGEN_PASSWORD_HASH_TIME_COST` | Passes over that memory for one admin password check (Argon2id). |
| `password_hash.parallelism` | Password hashing threads | whole number (threads) | 1 | after a restart | command line only | `PLANETGEN_PASSWORD_HASH_PARALLELISM` | Threads one admin password check uses (Argon2id). |
| `proxy_fix.x_for` | Proxies for the visitor's address | whole number (proxies) | 0 | after a restart | command line only | `PLANETGEN_PROXY_FIX_X_FOR` | Reverse proxies to trust for X-Forwarded-For, the client address (0 = off). |
| `proxy_fix.x_proto` | Proxies for the scheme | whole number (proxies) | 0 | after a restart | command line only | `PLANETGEN_PROXY_FIX_X_PROTO` | Reverse proxies to trust for X-Forwarded-Proto, the scheme (0 = off). |
| `proxy_fix.x_host` | Proxies for the host name | whole number (proxies) | 0 | after a restart | command line only | `PLANETGEN_PROXY_FIX_X_HOST` | Reverse proxies to trust for X-Forwarded-Host (0 = off). |
| `cache.max_mb` | Disk cache size | number (MB) | 1024 | at once | any admin | `PLANETGEN_CACHE_MAX_MB` | How much disk the site's stored caches (today the Galaxy Map's nodes and stages) may use, in all. When it is full the warm-up stops, the nodes it built ahead stay, and the nodes a visit asked for are dropped, the least recently visited first (0 turns the disk caches off). |
| `cache.warm_on_change` | Rebuild caches after changes | true or false | true | at once | any admin | `PLANETGEN_WARM_ON_CHANGE` | Start the warm-up by itself after anything that changes what the site shows (a generation run, the orbital update, an edit), so the caches are rebuilt without waiting for the next update. It rebuilds only what the change touched. |
| `cache.warm_settle_seconds` | Wait before rebuilding | number (seconds) | 60 | at once | any admin | none | How long the automatic warm-up waits after the first change before it starts, so a burst of edits is one warm-up and not one for each. |
| `tile_cache.dir` | Map cache folder | text | empty | at once | command line only | `PLANETGEN_TILE_CACHE_DIR` | Where the Galaxy Map's nodes are cached on disk; empty means /var/cache/planetgen/tiles. |
| `tile_cache.warm_zoom_floor` | Closest map block built ahead | whole number (shells) | 9 | at once | any admin | none | Finest block level of the zoom run `planetgen warm` builds from the opening view in to the galaxy's centre, on top of the populated blocks it builds in full (9 is blocks of nine shells, about 36 pc across; the run covers blocks of that size and bigger; 243 or more builds no run). |
| `tile_cache.prefetch_per_visit` | Nodes prefetched per visit | whole number (nodes) | 128 | at once | any admin | none | Most nodes one request for nodes queues to be built in the background: the children of the nodes it asked for, then their parents and the parents' other children (0 turns the background build off). |
| `tile_cache.prefetch_per_minute` | Nodes prefetched per minute | whole number (nodes) | 600 | at once | any admin | none | Most nodes the background build queues in a minute, all requests together. |
| `tile_cache.prefetch_max_queued` | Prefetch queue limit | whole number (jobs) | 8 | at once | any admin | none | The background build queues nothing while this many of its jobs (16 nodes each) or other upkeep work already wait. |
| `page_cache.enabled` | Page answer cache | true or false | true | after a restart | any admin | `PLANETGEN_PAGE_CACHE` | Keep the API's public answers in memory, per WSGI process. |
| `page_cache.max_entries` | Cached answers | whole number (answers) | 2000 | after a restart | any admin | none | Most answers kept. |
| `page_cache.max_mb` | Answer cache memory | number (MB) | 64 | after a restart | any admin | none | Most memory the answers may take. |
| `page_cache.stamp_seconds` | Change check interval | number (seconds) | 15 | after a restart | any admin | none | How often the galaxy's content stamp is checked. |
| `page_cache.max_age_seconds` | Answer lifetime in memory | number (seconds) | 300 | after a restart | any admin | none | Nothing is kept longer than this. |
| `page_cache.shared_seconds` | Answer lifetime in Redis | number (seconds) | 21600 | after a restart | any admin | none | How long an answer is kept in Redis, where every web process and the warm-up share it. Each answer is served only while the data it read has not changed, so this only bounds the memory. |
| `page_cache.warm_max_objects` | Pages built ahead | whole number (objects) | 5000 | at once | any admin | none | Most sector pages, system pages and phenomenon pages the warm-up builds (each, nearest the galaxy centre first); the rest are built on their first visit, and the warm-up says how many it left out. |
| `page_cache.stored_counts` | Use stored table counts | true or false | true | at once | any admin | `PLANETGEN_COUNT_CACHE` | Table pages show stored totals and filter-menu counts, refreshed in the background, instead of counting whole tables in the request (off: exact counts at any price). |
| `jobs.dir` | Job files folder | text | empty | at once | command line only | `PLANETGEN_JOBS_DIR` | Where each background job's command lines, status and output are kept; empty means /var/lib/planetGen/jobs. |
| `jobs.keep` | Finished jobs kept | whole number (jobs) | 20 | at once | any admin | none | How many finished jobs are kept. |
| `jobs.python` | Python for jobs | text | empty | at once | command line only | `PLANETGEN_PYTHON` | The interpreter that runs planetgen and planetgen.cli.reset for the jobs; empty means the web app's own Python. Editing it from the web would be code execution, so it is file-only. |
| `parity.enabled` | Sector repair file | true or false | false | after a restart | any admin | `PLANETGEN_PARITY` | Keep a parity file beside the database that can rebuild a damaged sector: every saved or edited sector updates the parity of its group, which costs a few milliseconds a save and about 7 percent of the sectors' size in disk space. Off, saves write no parity and `planetgen refresh-parity` builds it later. |
| `parity.dir` | Repair file folder | text | empty | at once | command line only | `PLANETGEN_PARITY_DIR` | Where the parity file is kept; empty means /var/lib/planetGen/parity. Keep it off the disk the database is on if you can: a parity file on the same failing disk protects nothing. |
| `parity.group_size` | Sectors per repair group | whole number (sectors) | 32 | after a restart | any admin | none | Sectors in one parity group (G). Any `parity_blocks` damaged sectors of a group are rebuilt from the rest, which are read to do it: a bigger group costs less disk and more reading. A change applies to the groups made afterwards. |
| `parity.parity_blocks` | Repair blocks per group | whole number (blocks) | 2 | after a restart | any admin | none | Parity blocks per group (m): how many damaged sectors of one group can be rebuilt. A change applies to the groups made afterwards. |

<!-- END GENERATED OPTIONS -->

## Details of each option

| Field | Purpose |
|---|---|
| `site_name` | Display name for this deployment, e.g. `"planetGen"`. The pages show it in the header, the home page heading, every `<title>`, and the footer (when it isn't the default). |
| `base_url` | The base URL this deployment is served from, e.g. `"http://localhost/"` or `"https://planetgen.example.com/"`. Not yet read by any page -- reserved for future absolute-URL generation (e.g. constructing shareable links) that can't be derived from a request alone. |
| `api_base_url` | Base URL of the Flask API's `/api` mount point, the address `update.sh` asks whether the site answers (`api_base_url` plus `/health`) and the one `../src/planetgen/web/lib/apiclient.py` calls outside the Flask app (the pages call the API in-process). There is no default and nothing falls back to a loopback address (OPS.86): set it to the address Apache or gunicorn serves the API at, for example `https://example.org/api`; the update stops before it changes the database when it is empty. Equivalent to `PLANETGEN_API_BASE_URL`. |
| `debug` | Debug mode. When true, every part of planetGen writes everything it decides: each `planetgen` command (and each other command: `python3 -m planetgen.cli.*`) writes every decision, random draw (with its source line and the probabilities it refers to), SQL statement, run, stage, unit, worker and cache event to **stdout**, with an ISO timestamp to the millisecond, and the same lines go to `log_file` along with the `html/` pages', the API's and the workers' (every web request, API call and admin access check, and every error with its traceback). The lines of the queue workers that do a run's work are relayed to the stdout of the command that started the run, as they are written. A job the Generate page starts shows them in its job log. The output is large (a single system writes about 265 KB) and the run takes a quarter to two thirds longer, so `install.sh`/`update.sh` install a logrotate config for the file (daily, or past 100 MB; 7 compressed copies kept); `debug_areas` limits what is written. A command whose result goes to stdout (`--output -`) writes its debug lines to stderr; `--quiet` still silences the console. The web interface never shows tracebacks on its pages, debug or not; with debug on, a 500 page says the traceback is in the debug log. See `../src/planetgen/util/log.py`. Equivalent to `PLANETGEN_DEBUG` (`1` on, `0` off). |
| `debug_areas` | Limits debug mode to some areas: `decision`, `draw`, `sql`, `queue`, `cache`, `web`, `lib` (comma separated; a leading `-` leaves one out, as in `-draw`). Empty, the default, writes everything but `branch`, which logs every conditional jump in planetGen's own code (hundreds of thousands of lines for one system, about 15 times the run time) and is written only when named (`branch`, or `all,branch`). Equivalent to `PLANETGEN_DEBUG_AREAS` and to `--debug-areas`. |
| `log_file` | Where the debug log goes; default `/var/log/planetgen.log`. With `debug` on, `install.sh`/`update.sh` create it (owned by Apache's user and group, mode 0660) and point the logrotate config at it. A login user who runs the generator from a shell must be in Apache's group (`sudo usermod -aG www-data <user>`, then log in again) to append to it; root can always write to it. A process that can't open it carries on without a debug log and prints one warning to stderr. Equivalent to `PLANETGEN_LOG_FILE`. |
| `log_dir` | The folder of the always-on **activity log** (`planetgen.log`; see [The activity log](#the-activity-log) below). Empty (the default) means the platform's standard place: `/var/log/planetgen` on Linux, `/Library/Logs/planetgen` on macOS. `install.sh`/`update.sh` create it (root and Apache's group, mode 2770; the file Apache's user and group, mode 0660). Equivalent to `PLANETGEN_LOG_DIR`. |
| `log_rotation` | How the activity log is rotated: `"auto"` (the default), `"system"` or `"app"`. `"system"` leaves it to logrotate or newsyslog (the program reopens the file after it is moved); `"app"` makes the program rotate it itself (100 MB per file, 30 old copies). `"auto"` means `"system"` where `install.sh`/`update.sh` installed `/etc/logrotate.d/planetgen-log` or `/etc/newsyslog.d/planetgen-log.conf`, else `"app"`. No environment variable. |
| `mysql.host`/`mysql.port`/`mysql.user`/`mysql.password`/`mysql.database` | The one MySQL connection every entry point uses -- the generation CLIs, `install.sh`/`planetgen.cli.migrate`, and the Flask API's reads and writes alike (see `planetgen.db.store.MySQLConfig`). There's no separate write-capable override: give this account whatever grants the most demanding caller needs. Equivalent to `PLANETGEN_MYSQL_HOST`/`_PORT`/`_USER`/`_PASSWORD`/`_DATABASE`. |
| `mysql.statement_timeout_seconds` | The longest any one statement on the web interface's and API's read-only connections may run, in seconds (default 10; 0 turns the limit off). Uses MySQL's `MAX_EXECUTION_TIME`. A query stopped by it gives a 504 "Took too long" page, or a `QUERY_TIMEOUT` error from the API, instead of holding a web worker until Apache's request timeout. Generation, migrations and admin writes are never limited. Separately, the web interface's and API's control-database connections wait at most 5 seconds for a row another transaction has locked (`WEB_LOCK_WAIT_SECONDS` in `tuning.py`, not an option); a login or session refresh then goes ahead without noting its last-login time, and other requests answer a 503 "Busy" page or API error that retries by itself. |
| `mysql.database_prefix` | The schema-name prefix `list_databases`/`resolve_database` filter by, for a deployment with more than one game database on one server (e.g. `planetgen`, `planetgen_alpha`, ...). Equivalent to `PLANETGEN_MYSQL_DATABASE_PREFIX`. |
| `control_database` | Name of the separate MySQL schema holding admin identities/sessions/API keys/audit log (see `../src/planetgen/db/control_schema.sql`), reached via the account above. Never listed by `/api/databases` or selectable with `?db=`, even when it shares `mysql.database_prefix` (as the default `planetgen_control` does). Equivalent to `PLANETGEN_CONTROL_DATABASE`. |
| `redis.url` | The Redis server the work queue runs on (PERF.24): a `planetgen` run with more than one worker queues its tasks there (without Redis it runs one sector at a time), and the admin Generate page's jobs run there (without Redis none starts). Default `"redis://127.0.0.1:6379/0"`; `PLANETGEN_REDIS_URL` overrides it. `install.sh`/`update.sh` check that it answers and, on Linux with a local address, install and start `redis-server` when it doesn't; on macOS they say to run `brew install redis && brew services start redis`. |
| `ratelimit.default`/`ratelimit.storage_uri` | Flask-Limiter's default rate limit and storage backend for the API (see `../src/planetgen/api/config.py`). `storage_uri` is empty by default, which means the Redis server in `redis.url`; the login lockouts (SEC.1, SEC.21) count there too (SEC.30), so every worker shares them and a restart keeps them. `memory://` counts in each process instead (right with exactly one worker). A Redis outage falls back to counting in memory. Equivalent to `PLANETGEN_RATELIMIT_DEFAULT`/`PLANETGEN_RATELIMIT_STORAGE_URI`. |
| `ratelimit.pages` | Per-client-IP limits on the HTML pages and `/api/health` (see `../src/planetgen/api/limiter.py` and [`html-interface.md`](html-interface.md#rate-limits)): `search` (`/search` and `/galaxy/locate`, default `"30 per minute"`), `galaxy` (`/galaxy`, `"60 per minute"`), `galaxy_tiles` (the Galaxy Map's data and scene routes: `/galaxy/tiles`, `/galaxy/stage`, `/galaxy/window`, `/galaxy/near`, `/galaxy/territories`, the uncharted and nebula routes and `/sector/<id>/scene`, `"600 per minute"`), `health` (`/api/health`, `"60 per minute"`) and `other` (every other page, counted together, `"300 per minute"`). A name left out keeps its default; an empty string turns that limit off. No environment variable. |
| `proxy_fix.x_for`/`proxy_fix.x_proto`/`proxy_fix.x_host` | How many reverse proxies to trust for the client address (`X-Forwarded-For`), the scheme (`X-Forwarded-Proto`) and the host (`X-Forwarded-Host`); werkzeug's `ProxyFix`, applied in `create_app` (`../src/planetgen/web/app.py`). All `0` (the default) means off: the app uses the connection's own address and scheme, which is right under Apache + mod_wsgi. Behind nginx, Caddy or Apache's `mod_proxy` in front of gunicorn, set `x_for` and `x_proto` to `1`: otherwise every visitor shares the proxy's address (one rate-limit budget for everyone) and the app never sends `Strict-Transport-Security`. Only turn it on when the app server is reachable through that proxy alone (a Unix socket or `127.0.0.1`), since whoever connects can set these headers. A value that isn't a whole number 0 or above stops the app at startup. See [`deployment/README.md`](../deployment/README.md#behind-a-reverse-proxy-proxy_fix). Equivalent to `PLANETGEN_PROXY_FIX_X_FOR`/`PLANETGEN_PROXY_FIX_X_PROTO`/`PLANETGEN_PROXY_FIX_X_HOST`. |
| `login_allowlist` | Client addresses or networks (`"198.51.100.7"`, `"10.0.0.0/8"`, IPv6 too) that are never locked out after failed logins (SEC.1); loopback (`127.0.0.0/8`, `::1`) never is either. Put your own fixed address here if you'd rather never be locked out, and a reverse proxy's address if `proxy_fix` can't be used. A bad entry is skipped with a warning. Equivalent to `PLANETGEN_LOGIN_ALLOWLIST` (comma- or space-separated). |
| `admin_cookie_insecure` | When true, the admin session cookie is sent over plain HTTP. Only for local development without TLS in front of the app (e.g. `python src/html/wsgi.py`) -- a production deployment must never set this. With it the cookie is named `pg_admin_session`; otherwise it is `__Host-pg_admin_session` (SEC.32). Equivalent to `PLANETGEN_ADMIN_COOKIE_INSECURE=1`. |
| `secret_key` | Signs the CSRF tokens on the Flask-served pages' forms (`src/planetgen/web/csrf.py`). Set it to a long random string (e.g. `python3 -c "import secrets; print(secrets.token_hex(32))"`) and keep it private. If empty, the app makes a random one at startup and logs a warning: forms still work, but a form left open across an app restart fails once and has to be resubmitted. Equivalent to `PLANETGEN_SECRET_KEY`. |
| `password_hash.memory_kib`/`password_hash.time_cost`/`password_hash.parallelism` | The cost of the Argon2id hash an admin password is stored as (SEC.32): memory in KiB, passes over it, and threads. The defaults, 19456 (19 MiB), 2 and 1, are OWASP's minimum for Argon2id, about a tenth of a second on a small host; raise the memory where the host has room (each sign-in attempt holds that much while it checks, and the app server runs several at once). A hash stored with other numbers (or by an older release, as PBKDF2 or scrypt) is replaced at that admin's next successful sign-in. Equivalent to `PLANETGEN_PASSWORD_HASH_MEMORY_KIB`/`PLANETGEN_PASSWORD_HASH_TIME_COST`/`PLANETGEN_PASSWORD_HASH_PARALLELISM`. |
| `tile_cache.dir` | Where the web pages (the Flask-served `/galaxy` and `/galaxy/nodes`, in the WSGI daemon) keep their on-disk cache of 3D Galaxy Map nodes (see `../src/planetgen/web/lib/tilecache.py`), (its size is `cache.max_mb`, below). An empty `dir` (the default) means `/var/cache/planetgen/tiles`, which `install.sh`/`update.sh` create for Apache's user (via `examples/apache/create-cache-dir.sh`, which also creates a `dir` you set here); if it's missing and Apache can't create it, the cache falls back to a private (mode 0700) `planetgen-tiles` folder in the system temp directory; an existing one that isn't a real directory owned by Apache's user, or that other users can write to, is refused (the cache is then off) rather than reused. Equivalent to `PLANETGEN_TILE_CACHE_DIR`. |
| `cache.max_mb` | How much disk the stored caches may use in all (default 1024; `0` turns the disk caches off; `PLANETGEN_CACHE_MAX_MB`; it replaces `tile_cache.max_mb`, which is no longer read). Today only the Galaxy Map's nodes are stored; the object index and the stored lists will share it and are kept first. When it is full `planetgen warm` stops (it says how many nodes it left out), the nodes it built ahead stay, and the nodes a visit asked for are dropped, the least recently visited first (PERF.133, PERF.139). |
| `page_cache.enabled`/`max_entries`/`max_mb`/`stamp_seconds`/`max_age_seconds` | The defaults are `true`, 2000, 64, 15 and 300, as in `config.json.example`. The web pages keep the API's public answers in memory, per WSGI process (see `../src/planetgen/web/lib/pagecache.py`), so a repeat visit doesn't query the database. Any successful write through the API clears it; a new galaxy stamp (sectors or systems added, edited or deleted, by `planetgen` jobs too) is noticed within `stamp_seconds`; nothing is kept longer than `max_age_seconds`. `PLANETGEN_PAGE_CACHE=off` (or `enabled: false`) turns it off. `stored_counts` (default `true`; `PLANETGEN_COUNT_CACHE` overrides it) is separate: table pages show stored totals and filter-menu counts, refreshed in the background, instead of counting rows on every page view (`../src/planetgen/db/countcache.py`). |
| `jobs.dir`/`jobs.keep`/`jobs.python` | The admin Generate page's background jobs (see `../src/planetgen/web/jobs.py`). `dir` is where each job's command lines, status and output are kept; empty (the default) means `/var/lib/planetGen/jobs`, which `install.sh`/`update.sh` create for Apache's user (via `examples/apache/create-cache-dir.sh`), falling back to a private (mode 0700) `planetgen-jobs` folder in the system temp directory, refused under the same conditions as the tile cache's fallback. `keep` is how many finished jobs are kept (default 20). `python` is the interpreter that runs `planetgen` and `planetgen.cli.reset`; empty means the web app's own Python (under mod_wsgi, `<sys.prefix>/bin/python3`). Equivalent to `PLANETGEN_JOBS_DIR`/`PLANETGEN_PYTHON`; `keep` has no environment variable. |

## The activity log

Separate from the debug log, planetGen always writes a short record of
who did what to `planetgen.log` in `log_dir` (SEC.28,
`../src/planetgen/admin/activity_log.py`): sign-ins, refused requests,
every change to the database and every API call. It is on whatever `debug` says; with
`debug` on, each line is copied into the debug log too. Passwords,
session tokens and API keys are never written.

Every line has the same shape, so a tool such as fail2ban can match it
([`deployment/fail2ban.md`](../deployment/fail2ban.md) has a ready filter and jail):

```text
2026-10-01T08:00:00Z planetgen[1234]: AUTH login.failed ip=203.0.113.5 user="admin"
2026-10-01T08:00:09Z planetgen[1234]: DB sector.update ip=203.0.113.5 user="boss" target=sector:42 detail="{'name': 'Home'}"
```

That is: the time in UTC, `planetgen[<process id>]:`, the category, the
event, `ip=` (the client address, or `-` with none, as on the command
line), `user=` (always quoted; the admin, the username as typed, or the
login name that ran a command; `-` when nobody), then any further
`key=value` fields. The address always comes before anything a visitor
typed and is written only when it is a real IPv4 or IPv6 address (else
`-`). Every value that isn't a plain number or token is quoted, with `"`
and `\` escaped and control characters written as `\xNN`, so a crafted
username can't start a new line or fake a field; quoted values are cut
at 200 characters.

| Category | Events |
|---|---|
| `AUTH` | `login.ok`; `login.failed` (wrong username or password); `login.locked` (refused unchecked while the address or username is locked, with `scope=ip` or `scope=user` and `retry_after=`); `lockout.start` (a failure that started a lock, with `scope=`, `subject=` and `seconds=`); `login.ratelimited` (past the per-address limit); `logout`; `credentials.changed` (with `old_user=` after a rename); `password.failed` (a wrong current password on `/account`); `login.password_ok` (right password, two-factor code still to come); `totp.failed` (a wrong two-factor code); `totp.recovery_used` (with `left=`); `totp.enabled`, `totp.disabled`; `apikey.create`, `apikey.revoke` |
| `AUTHZ` | `login.required` (an admin route with no session or key); `session.invalid` (an expired, ended or unknown session cookie); `apikey.invalid` (an unknown or revoked API key); `credentials.unchanged` (an admin route refused until the first credentials are changed); `admin.required` (an admin-only page action refused); `csrf.failed` (a form token that doesn't match); each with `path=` |
| `DB` | every write through the web interface or API (`sector.create`, `system.update`, `star.rename`, `facility.create`, `lockout.lift`, ...: the same action names as `admin_audit_log`, with `target=` and `detail=`; `planetgen.cli.lockouts` writes `lockout.lift`, `devices.revoke` and `totp.reset` too); `migrate` (each schema migration step, `db=`, `from_version=`, `to_version=`) |
| `GEN` | `generate.start` and `generate.finish` for each `planetgen` run that writes the database (`command=`, `db=`; the finish line adds `status=ok/failed/interrupted`, `seconds=`, `sectors=`, `systems=`, `phenomena=`); `job.start` and `job.cancel` for the admin Generate page |
| `API` | `call` for every request to `/api/...` (API.15), written when the answer is ready: `user=` (the owner of the API key, the signed-in user, `god` for a command-line tool, or `-`), `via=` (`key`, `session`, `anonymous` or `console`), `method=`, `route=` (the matched route, `-` for none), `path=`, `status=` (the HTTP response code), `ms=`, `bytes_in=`, `bytes_out=`, `rid=` (also sent back as `X-Request-Id`), `ua=`, `client_api=` (the API version the client sent), and when they apply `key_id=` and `key_prefix=` (never the key), `internal=1` (a call a web page made to the API in-process, about six for each page view) and `ratelimited=1`. A command-line tool that does what an API route does writes the same line with `user="god"`, `via=console` and `status=200` or `500` |

The API lines use HTTP response codes, for parity with the web server's own logs; a request the web server refuses before it reaches the program (a 413 from `LimitRequestBody`, a failed TLS handshake) appears only in the web server's log. `api_log.internal` set to `off` leaves out the calls the web pages make to the API in-process; the default logs all of them. Request bodies, tokens and `Authorization` headers are never written. The address is the client's full address, kept as long as the file is (the rotation below: 30 copies of up to 100 MB).

Failed and locked sign-ins and wrong current passwords also go into the
control database's `admin_audit_log` (kept 90 days), which the admin
Stats page lists under "Failed sign-ins".

Rotation: `install.sh`/`update.sh` (as root, through
`examples/apache/setup-debug-log.sh`) install the system's own rotation,
so nobody has to edit these by hand. For a setup without the installer,
the samples are:

Linux, `/etc/logrotate.d/planetgen-log` (the user and group are the web
server's: `www-data` on Debian and Ubuntu, `apache` on RHEL and Fedora):

```conf
/var/log/planetgen/*.log {
    daily
    maxsize 100M
    rotate 30
    missingok
    notifempty
    compress
    delaycompress
    dateext
    dateformat -%Y%m%d-%s
    su root www-data
    create 0660 www-data www-data
}
```

macOS, `/etc/newsyslog.d/planetgen-log.conf` (newsyslog runs every hour
by itself; `$D0` rotates at midnight, the size column in KB rotates
sooner past 100 MB, `J` compresses with bzip2, `N` means no process to
signal):

```conf
# logfilename                              [owner:group]  mode count size(KB) when flags
/Library/Logs/planetgen/planetgen.log      _www:_www      660  30    102400   $D0  JN
```

The program reopens the file after either tool moves it, so no restart
or `copytruncate` is needed. A process that can't open the file
carries on without it and prints one warning to stderr (only when the
folder exists, so a development checkout stays quiet). A login user who
runs `planetgen` from a shell needs to be in the web server's group to
append to it, as for the debug log.

## Why the real file is gitignored but the example isn't

`config.json` (the real, per-deployment file) is listed in
[`../.gitignore`](../../.gitignore), next to the existing `*.db`/`db.sqlite3`
entries -- it's deployment-specific configuration that holds real MySQL
credentials, and committing it would either leak those values or force
every deployment to share one repo-tracked file. `../config.json.example`
is the opposite: a template with placeholder values, meant to be
committed so a fresh checkout always has a shape to copy from, exactly
the same split already used for
`../examples/apache/planetgen.conf.example` (committed template, edited
into a site-specific vhost file that itself doesn't live in the repo).

## Setup

Copy the template at the repo root and edit it for this deployment:

```bash
cp config.json.example config.json
```

Then fill in `mysql.host`/`mysql.user`/`mysql.password`/`mysql.database`
(and `control_database`, if this deployment names its control schema
something other than the default -- see
[`deployment/README.md`](../deployment/README.md#mysql-accounts)) for this
deployment's actual database server. `site_name`/`base_url` are cosmetic
and safe to leave as-is.

If `config.json` doesn't exist at all, `planetgen.util.settings.get_settings()`
falls back to the built-in defaults shown in `config.json.example`, so
every entry point keeps working (against `127.0.0.1:3306` as user
`planetgen`) without this step -- exactly the environment-variable-only
behavior this project had before `config.json` existed.

## Relationship to `PLANETGEN_*` environment variables

`config.json` doesn't replace the `PLANETGEN_*` environment variables --
see "Precedence" above. In short: `config.json` is where a deployment
sets its defaults once; environment variables (the web server's or
service's own environment -- not an Apache `SetEnv`, which never reaches
`os.environ` under mod_wsgi -- a systemd `EnvironmentFile`, a launchd
`EnvironmentVariables`, or a shell export) remain the way to override
one of those defaults for a single process, most importantly the
per-instance database name `planetgen-orbits@.service` needs.

### `PLANETGEN_WORKERS`: not a `config.json` option

The number of worker processes a generation run uses is set by `--workers N` on the `planetgen` command, else by the
`PLANETGEN_WORKERS` environment variable (a number, or `auto`); with neither, a run uses 80% of the cores (one fewer when
the database is on the same machine, so it keeps a core). A positive number is used as given. More than one worker needs
the Redis server in `redis.url`; without it a run generates one sector at a time. The admin Generate page's jobs run `planetgen`, so
they use it too. The test suite sets it to 1 ([`testing.md`](testing.md)).
