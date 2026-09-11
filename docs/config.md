# Deployment Configuration

This document describes `config.json`, the single per-deployment
configuration file for the whole planetGen project -- the generation
CLIs (`sectorGen.py`/`systemGen.py`/`galaxyGen.py`/`queryDb.py`/
`migrateDb.py`/`updateOrbits.py`), the Flask API (`../src/html/api/`), and the
`html/` CGI browser all read it -- and
[`../config.json.example`](../config.json.example), the committed
template it's copied from.

## What it's for

Before this file existed, every one of those entry points read its own
handful of `PLANETGEN_*` environment variables (`stellarObjects._db.MySQLConfig`,
`../src/html/api/config.py`, `../src/html/lib/apiclient.py`, `../src/html/lib/page.py`),
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

`config.json` lives at the repo root -- a sibling of `../src/html/`, `db/`,
and `src/` -- rather than inside `../src/html/` itself (its *example*
template, [`../config.json.example`](../config.json.example), lives at the
repo root too, since a template holds only placeholder values, not
secrets). This mirrors how Apache's `DocumentRoot` for this application is
`../src/html/` alone (see
[`../examples/apache/planetgen.conf.example`](../examples/apache/planetgen.conf.example)),
so `config.json`, one level above `../src/html/`, can never be requested
over HTTP no matter how the vhost or `.htaccess` rules are written --
there's no path traversal or misconfiguration that reaches it, because
it's outside the tree Apache serves at all. This matters because
`config.json` holds real MySQL credentials, unlike the placeholder fields
its now-removed predecessor (`webconfig.json`) only ever reserved for
them.

## Precedence

Every setting below has up to three sources, checked in this order:

1. An explicit function/CLI argument (e.g. `sectorGen.py --mysql-password ...`,
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
| `site_name` | Display name for this deployment, e.g. `"planetGen"`. Not yet read by any page -- reserved for a future UI element (page title, header) once the web interface wants to show it. |
| `base_url` | The base URL this deployment is served from, e.g. `"http://localhost/"` or `"https://planetgen.example.com/"`. Not yet read by any page -- reserved for future absolute-URL generation (e.g. constructing shareable links) that can't be derived from a CGI request alone. |
| `api_base_url` | Base URL of the Flask API's `/api` mount point that the `html/` CGI browser talks to (see `../src/html/lib/apiclient.py`). Defaults to `http://127.0.0.1/api`; override when the API is deployed at a different host/port. Equivalent to `PLANETGEN_API_BASE_URL`. |
| `debug` | When true, an unhandled exception in the `html/` CGI browser renders its traceback on the page instead of a generic 500 (see `../src/html/lib/page.py`). Never enable this in production -- it leaks file paths and query text. Equivalent to `PLANETGEN_DEBUG`. |
| `mysql.host`/`mysql.port`/`mysql.user`/`mysql.password`/`mysql.database` | The default (read-only-capable) MySQL connection -- used by every entry point that doesn't need a different account (see `stellarObjects._db.MySQLConfig`). Equivalent to `PLANETGEN_MYSQL_HOST`/`_PORT`/`_USER`/`_PASSWORD`/`_DATABASE`. |
| `mysql.database_prefix` | The schema-name prefix `list_databases`/`resolve_database` filter by, for a deployment with more than one game database on one server (e.g. `planetgen`, `planetgen_alpha`, ...). Equivalent to `PLANETGEN_MYSQL_DATABASE_PREFIX`. |
| `mysql_write.host`/`.port`/`.user`/`.password`/`.database` | Overrides for the write-capable account the Flask API's write endpoints (and `/api/auth/`) use -- an empty string (the default) for any field means "inherit the matching `mysql` value", so a single-account local/dev setup needs no extra configuration. Equivalent to `PLANETGEN_MYSQL_WRITE_HOST`/etc. See [`apache-deployment.md`](apache-deployment.md#mysql-accounts). |
| `control_database` | Name of the separate MySQL schema holding admin identities/sessions/API keys/audit log (see `../src/stellarObjects/control_schema.sql`), reached via the write-capable account above. Equivalent to `PLANETGEN_CONTROL_DATABASE`. |
| `ratelimit.default`/`ratelimit.storage_uri` | Flask-Limiter's default rate limit and storage backend for the API (see `../src/html/api/config.py`). `storage_uri` needs a shared backend (e.g. `redis://...`) once a deployment runs more than one worker process. Equivalent to `PLANETGEN_RATELIMIT_DEFAULT`/`PLANETGEN_RATELIMIT_STORAGE_URI`. |
| `admin_cookie_insecure` | When true, the admin session cookie is sent over plain HTTP. Only for local development without TLS in front of the app (e.g. `python src/html/wsgi.py`) -- a production deployment must never set this. Equivalent to `PLANETGEN_ADMIN_COOKIE_INSECURE=1`. |

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
(and `mysql_write.*`/`control_database` if this deployment uses a
separate write-capable account -- see
[`apache-deployment.md`](apache-deployment.md#mysql-accounts)) for this
deployment's actual database server. `site_name`/`base_url` are cosmetic
and safe to leave as-is.

If `config.json` doesn't exist at all, `stellarObjects.appconfig.load_config()`
falls back to the built-in defaults shown in `config.json.example`, so
every entry point keeps working (against `127.0.0.1:3306` as user
`planetgen`) without this step -- exactly the environment-variable-only
behavior this project had before `config.json` existed.

## Relationship to `PLANETGEN_*` environment variables

`config.json` doesn't replace the `PLANETGEN_*` environment variables --
see "Precedence" above. In short: `config.json` is where a deployment
sets its defaults once; environment variables (Apache `SetEnv`, a
systemd `EnvironmentFile`, or a shell export) remain the way to override
one of those defaults for a single process, most importantly the
per-instance database name `planetgen-orbits@.service` needs.
