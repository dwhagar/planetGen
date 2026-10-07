# Deployment Configuration

This document describes `config.json`, the single per-deployment
configuration file for the whole planetGen project -- the command-line
tools (`generate.py`, `planetgen.db.query`, `planetgen.cli.migrate`, `planetgen.cli.orbits`,
`planetgen.cli.reset`, ...) and the Flask app (the API in
`../src/planetgen/api/` and the pages in `../src/html/web/`) all read it -- and
[`../config.json.example`](../config.json.example), the committed
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

It is loaded by [`../src/stellarObjects/appconfig.py`](../src/stellarObjects/appconfig.py),
a small, dependency-free (standard library `json`/`os`/`copy` only) loader
module, and merged onto its own built-in defaults -- a `config.json` that
only sets the fields it needs to change is enough; anything it omits
keeps its default value.

## Location: repo root

`config.json` lives at the repo root -- next to `src/` -- rather than
inside `../src/html/` itself (its *example*
template, [`../config.json.example`](../config.json.example), lives at the
repo root too, since a template holds only placeholder values, not
secrets). This mirrors how Apache's `DocumentRoot` for this application is
`../src/html/` alone (see
[`../examples/apache/planetgen.conf.example`](../examples/apache/planetgen.conf.example)),
and the other web servers only serve `src/html/static/` from disk (see
[`deployment/`](deployment/README.md)), so `config.json`, outside those
trees, can never be requested over HTTP no matter how the site's rules
are written. This matters because
`config.json` holds real MySQL credentials, unlike the placeholder fields
its now-removed predecessor (`webconfig.json`) only ever reserved for
them.

It also has to be unreadable by other local users. `install.sh` and
`update.sh` (through `examples/apache/set-permissions.sh`) set it to
`root:<apache group>`, mode `640`, whenever it exists: root can edit it,
the web app (running as that group's user, `www-data` for Apache and the
gunicorn unit alike) can read it, nobody else can. The macOS and Windows
guides set the same rule by hand. A login user who runs the
generator or the maintenance scripts from a shell without `sudo` needs to
be in Apache's group (`sudo usermod -aG www-data <user>`, then log in
again) to read it. After creating or copying a new `config.json`, re-run
`sudo examples/apache/set-permissions.sh` (or `sudo ./update.sh`), or set
it by hand: `sudo chown root:www-data config.json && sudo chmod 640 config.json`.

## Precedence

Every setting below has up to four sources, checked in this order:

1. An explicit function/CLI argument (e.g. `generate.py sector --mysql-password ...`,
   or a test building its own `MySQLConfig(...)` directly) -- always wins.
2. The matching `PLANETGEN_*` environment variable, if set.
3. `config.json`.
4. The built-in default (matching `config.json.example`'s values).

Environment variables are still checked *before* `config.json` rather than
being replaced by it, so a single shared `config.json` can still be
overridden per-process where that's needed -- most notably
`planetgen-orbits@.service`'s per-instance `PLANETGEN_MYSQL_DATABASE=%i`
(see [`../examples/maintenance/`](../examples/maintenance/)), which relies on
one templated systemd unit overriding just the database name per timer
instance while everything else comes from the shared credentials file (or
now, `config.json`).

## Fields

| Field | Purpose |
|---|---|
| `site_name` | Display name for this deployment, e.g. `"planetGen"`. The pages show it in the header, the home page heading, every `<title>`, and the footer (when it isn't the default). |
| `base_url` | The base URL this deployment is served from, e.g. `"http://localhost/"` or `"https://planetgen.example.com/"`. Not yet read by any page -- reserved for future absolute-URL generation (e.g. constructing shareable links) that can't be derived from a request alone. |
| `api_base_url` | Base URL of the Flask API's `/api` mount point, used by `../src/planetgen/web/lib/apiclient.py` only outside the Flask app (the pages call the API in-process). Defaults to `http://127.0.0.1/api`; override when the API is deployed at a different host/port. Equivalent to `PLANETGEN_API_BASE_URL`. |
| `debug` | When true, every part of planetGen -- the generator CLI, the maintenance scripts, the `html/` pages and the API -- writes a verbose debug log to `log_file`: every decision the generator makes and why, every random roll (with the source line that asked for it and the probabilities or thresholds that line refers to), every SQL statement, every web request, API call and admin access check, and every error with its traceback, each line timestamped to the millisecond and tagged with its process. Off when missing. The log grows fast (a single sector writes megabytes), so `install.sh`/`update.sh` install a logrotate config for it (daily, or past 100 MB; 7 compressed copies kept). The web interface never shows tracebacks on its pages, debug or not; with debug on, a 500 page says the traceback is in the debug log. See `../src/stellarObjects/log.py`. Equivalent to `PLANETGEN_DEBUG` (`1` on, `0` off). |
| `log_file` | Where the debug log goes; default `/var/log/planetgen.log` (set it explicitly on Windows). With `debug` on, `install.sh`/`update.sh` create it (owned by Apache's user and group, mode 0660) and point the logrotate config at it. A login user who runs the generator from a shell must be in Apache's group (`sudo usermod -aG www-data <user>`, then log in again) to append to it; root can always write to it. A process that can't open it carries on without a debug log and prints one warning to stderr. Equivalent to `PLANETGEN_LOG_FILE`. |
| `log_dir` | The folder of the always-on **activity log** (`planetgen.log`; see [The activity log](#the-activity-log) below). Empty (the default) means the platform's standard place: `/var/log/planetgen` on Linux, `/Library/Logs/planetgen` on macOS, `logs` under the checkout on Windows. `install.sh`/`update.sh` create it (root and Apache's group, mode 2770; the file Apache's user and group, mode 0660); `install.ps1`/`update.ps1` create it and give the app's account write access. Equivalent to `PLANETGEN_LOG_DIR`. |
| `log_rotation` | How the activity log is rotated: `"auto"` (the default), `"system"` or `"app"`. `"system"` leaves it to logrotate or newsyslog (the program reopens the file after it is moved); `"app"` makes the program rotate it itself (100 MB per file, 30 old copies). `"auto"` means `"system"` where `install.sh`/`update.sh` installed `/etc/logrotate.d/planetgen-log` or `/etc/newsyslog.d/planetgen-log.conf`, else `"app"` (always on Windows). No environment variable. |
| `mysql.host`/`mysql.port`/`mysql.user`/`mysql.password`/`mysql.database` | The one MySQL connection every entry point uses -- the generation CLIs, `install.sh`/`planetgen.cli.migrate`, and the Flask API's reads and writes alike (see `planetgen.db.store.MySQLConfig`). There's no separate write-capable override: give this account whatever grants the most demanding caller needs. Equivalent to `PLANETGEN_MYSQL_HOST`/`_PORT`/`_USER`/`_PASSWORD`/`_DATABASE`. |
| `mysql.statement_timeout_seconds` | The longest any one statement on the web interface's and API's read-only connections may run, in seconds (default 10; 0 turns the limit off). Uses `max_statement_time` on MariaDB and `MAX_EXECUTION_TIME` on MySQL. A query stopped by it gives a 504 "Took too long" page, or a `QUERY_TIMEOUT` error from the API, instead of holding a web worker until Apache's request timeout. Generation, migrations and admin writes are never limited. |
| `mysql.database_prefix` | The schema-name prefix `list_databases`/`resolve_database` filter by, for a deployment with more than one game database on one server (e.g. `planetgen`, `planetgen_alpha`, ...). Equivalent to `PLANETGEN_MYSQL_DATABASE_PREFIX`. |
| `control_database` | Name of the separate MySQL schema holding admin identities/sessions/API keys/audit log (see `../src/planetgen/db/control_schema.sql`), reached via the account above. Never listed by `/api/databases` or selectable with `?db=`, even when it shares `mysql.database_prefix` (as the default `planetgen_control` does). Equivalent to `PLANETGEN_CONTROL_DATABASE`. |
| `redis.url` | The Redis server the work queue and the rate limits will use (OPS.21; nothing reads it yet), default `"redis://127.0.0.1:6379/0"`. `install.sh`/`update.sh` check that it answers and, on Linux with a local address, install and start `redis-server` when it doesn't; on macOS they say to run `brew install redis && brew services start redis`. `install.ps1`/`update.ps1` only check: Redis has no native Windows build, so use Memurai or Redis in WSL ([`deployment/windows.md`](deployment/windows.md#redis)). No environment variable. |
| `ratelimit.default`/`ratelimit.storage_uri` | Flask-Limiter's default rate limit and storage backend for the API (see `../src/planetgen/api/config.py`). `storage_uri` needs a shared backend (e.g. `redis://...`) once a deployment runs more than one worker process. Equivalent to `PLANETGEN_RATELIMIT_DEFAULT`/`PLANETGEN_RATELIMIT_STORAGE_URI`. |
| `ratelimit.pages` | Per-client-IP limits on the HTML pages and `/api/health` (see `../src/planetgen/api/limiter.py` and [`html-interface.md`](html-interface.md#rate-limits)): `search` (`/search` and `/galaxy/locate`, default `"30 per minute"`), `galaxy` (`/galaxy`, `"60 per minute"`), `galaxy_tiles` (`/galaxy/tiles`, `/galaxy/stage` and `/galaxy/territories`, `"600 per minute"`), `health` (`/api/health`, `"60 per minute"`) and `other` (every other page, counted together, `"300 per minute"`). A name left out keeps its default; an empty string turns that limit off. No environment variable. |
| `proxy_fix.x_for`/`proxy_fix.x_proto`/`proxy_fix.x_host` | How many reverse proxies to trust for the client address (`X-Forwarded-For`), the scheme (`X-Forwarded-Proto`) and the host (`X-Forwarded-Host`); werkzeug's `ProxyFix`, applied in `create_app` (`../src/planetgen/api/app.py`). All `0` (the default) means off: the app uses the connection's own address and scheme, which is right under Apache + mod_wsgi. Behind nginx, Caddy, IIS or Apache's `mod_proxy` in front of gunicorn or waitress, set `x_for` and `x_proto` to `1`: otherwise every visitor shares the proxy's address (one rate-limit budget for everyone) and the app never sends `Strict-Transport-Security`. Only turn it on when the app server is reachable through that proxy alone (a Unix socket or `127.0.0.1`), since whoever connects can set these headers. A value that isn't a whole number 0 or above stops the app at startup. See [`deployment/README.md`](deployment/README.md#behind-a-reverse-proxy-proxy_fix). Equivalent to `PLANETGEN_PROXY_FIX_X_FOR`/`PLANETGEN_PROXY_FIX_X_PROTO`/`PLANETGEN_PROXY_FIX_X_HOST`. |
| `login_allowlist` | Client addresses or networks (`"198.51.100.7"`, `"10.0.0.0/8"`, IPv6 too) that are never locked out after failed logins (SEC.1); loopback (`127.0.0.0/8`, `::1`) never is either. Put your own fixed address here if you'd rather never be locked out, and a reverse proxy's address if `proxy_fix` can't be used. A bad entry is skipped with a warning. Equivalent to `PLANETGEN_LOGIN_ALLOWLIST` (comma- or space-separated). |
| `admin_cookie_insecure` | When true, the admin session cookie is sent over plain HTTP. Only for local development without TLS in front of the app (e.g. `python src/html/wsgi.py`) -- a production deployment must never set this. Equivalent to `PLANETGEN_ADMIN_COOKIE_INSECURE=1`. |
| `secret_key` | Signs the CSRF tokens on the Flask-served pages' forms (`src/html/web/csrf.py`). Set it to a long random string (e.g. `python3 -c "import secrets; print(secrets.token_hex(32))"`) and keep it private. If empty, the app makes a random one at startup and logs a warning: forms still work, but a form left open across an app restart fails once and has to be resubmitted. Equivalent to `PLANETGEN_SECRET_KEY`. |
| `tile_cache.dir`/`tile_cache.max_mb` | Where the web pages (the Flask-served `/galaxy` and `/galaxy/tiles`, in the WSGI daemon) keep their on-disk cache of 3D Galaxy Map tiles (see `../src/planetgen/web/lib/tilecache.py`), and roughly how big that cache may grow before its oldest files are pruned. An empty `dir` (the default) means `/var/cache/planetgen/tiles` (set it explicitly on Windows), which `install.sh`/`update.sh` create for Apache's user (via `examples/apache/create-cache-dir.sh`, which also creates a `dir` you set here); if it's missing and Apache can't create it, the cache falls back to a private (mode 0700) `planetgen-tiles` folder in the system temp directory; an existing one that isn't a real directory owned by Apache's user, or that other users can write to, is refused (the cache is then off) rather than reused. `max_mb` of `0` turns the disk cache off. Equivalent to `PLANETGEN_TILE_CACHE_DIR`/`PLANETGEN_TILE_CACHE_MAX_MB`. |
| `page_cache.enabled`/`max_entries`/`max_mb`/`stamp_seconds`/`max_age_seconds` | Optional (not in `config.json.example`; the defaults are `true`, 2000, 64, 15 and 300). The web pages keep the API's public answers in memory, per WSGI process (see `../src/planetgen/web/lib/pagecache.py`), so a repeat visit doesn't query the database. Any successful write through the API clears it; a new galaxy stamp (sectors or systems added, edited or deleted, by `generate.py` jobs too) is noticed within `stamp_seconds`; nothing is kept longer than `max_age_seconds`. `PLANETGEN_PAGE_CACHE=off` (or `enabled: false`) turns it off. |
| `jobs.dir`/`jobs.keep`/`jobs.python` | The admin Generate page's background jobs (see `../src/html/web/jobs.py`). `dir` is where each job's command lines, status and output are kept; empty (the default) means `/var/lib/planetgen/jobs` (set it explicitly on Windows), which `install.sh`/`update.sh` create for Apache's user (via `examples/apache/create-cache-dir.sh`), falling back to a private (mode 0700) `planetgen-jobs` folder in the system temp directory, refused under the same conditions as the tile cache's fallback. `keep` is how many finished jobs are kept (default 20). `python` is the interpreter that runs `generate.py` and `planetgen.cli.reset`; empty means the web app's own Python (under mod_wsgi, `<sys.prefix>/bin/python3`). Equivalent to `PLANETGEN_JOBS_DIR`/`PLANETGEN_PYTHON`; `keep` has no environment variable. |
| `wiki.wikijs.base_url`/`.api_token` | The target Wiki.js instance's root URL and a Personal API Token (Admin -> API Access) -- see `../src/wikiClient/wikijs.py`. Leaving `base_url` empty (the default) means Wiki.js isn't offered as an "Upload to Wiki" target at all. Equivalent to `PLANETGEN_WIKIJS_BASE_URL`/`PLANETGEN_WIKIJS_API_TOKEN`. |
| `wiki.mediawiki.base_url`/`.username`/`.password` | The target MediaWiki instance's API entry point directory (everything up to, not including, `api.php`) and a [Bot Password](https://www.mediawiki.org/wiki/Special:BotPasswords) (`username` in `"User@BotName"` form) -- see `../src/wikiClient/mediawiki.py`. Leaving `base_url` empty (the default) means MediaWiki isn't offered as an upload target. Equivalent to `PLANETGEN_MEDIAWIKI_BASE_URL`/`PLANETGEN_MEDIAWIKI_USERNAME`/`PLANETGEN_MEDIAWIKI_PASSWORD`. |

## The activity log

Separate from the debug log, planetGen always writes a short record of
who did what to `planetgen.log` in `log_dir` (SEC.28,
`../src/stellarObjects/activitylog.py`): sign-ins, refused requests and
every change to the database. It is on whatever `debug` says; with
`debug` on, each line is copied into the debug log too. Passwords,
session tokens and API keys are never written.

Every line has the same shape, so a tool such as fail2ban can match it
([`deployment/fail2ban.md`](deployment/fail2ban.md) has a ready filter and jail):

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
| `GEN` | `generate.start` and `generate.finish` for each `generate.py` run that writes the database (`command=`, `db=`; the finish line adds `status=ok/failed/interrupted`, `seconds=`, `sectors=`, `systems=`, `phenomena=`); `job.start` and `job.cancel` for the admin Generate page |

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
or `copytruncate` is needed. On Windows there is no system rotation: the
program rotates `logs\planetgen.log` itself (`log_rotation` `"app"`),
keeping up to 30 numbered copies. A process that can't open the file
carries on without it and prints one warning to stderr (only when the
folder exists, so a development checkout stays quiet). A login user who
runs `generate.py` from a shell needs to be in the web server's group to
append to it, as for the debug log.

## Why the real file is gitignored but the example isn't

`config.json` (the real, per-deployment file) is listed in
[`../.gitignore`](../.gitignore), next to the existing `*.db`/`db.sqlite3`
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
[`deployment/README.md`](deployment/README.md#mysql-accounts)) for this
deployment's actual database server. `site_name`/`base_url` are cosmetic
and safe to leave as-is.

If `config.json` doesn't exist at all, `planetgen.util.appconfig.load_config()`
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
`EnvironmentVariables`, a WinSW `<env>`, or a shell export) remain the way to override
one of those defaults for a single process, most importantly the
per-instance database name `planetgen-orbits@.service` needs.
