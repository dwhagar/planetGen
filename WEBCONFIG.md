# WEBCONFIG.md

This document describes `webconfig.json`, the site-level configuration
file for the planetGen web interface ([`html/`](html/README.md)), and
`webconfig.json.example`, the committed template it's copied from.

## What it's for

Everything the web interface currently reads is either hard-coded or
supplied by Apache as an environment variable (`PLANETGEN_DB_DIR`,
`PLANETGEN_DEBUG` -- see [`html/README.md`](html/README.md)). There is no
place to record settings that describe *this deployment* -- a site name to
show in the UI, the base URL it's served from -- without either hard-coding
them or adding yet another Apache `SetEnv` line per setting. `webconfig.json`
is that place: a single JSON file, edited once per deployment, holding
settings that belong to the site itself rather than to how Apache happens
to invoke the CGI scripts.

It is loaded by [`stellarObjects/webconfig.py`](stellarObjects/webconfig.py),
a small, dependency-free (standard library `json`/`os` only) loader
module. **As of this writing, nothing calls it yet** -- it exists so the
config file's shape and loading behavior are settled before any feature
depends on them, not because any page currently reads `site_name` or
`base_url`.

## Location: repo root, not `html/`

`webconfig.json` lives at the repo root -- a sibling of `html/`, `db/`, and
`stellarObjects/` -- rather than inside `html/` itself. This deliberately
mirrors how `db/` is already kept out of the served webroot: Apache's
`DocumentRoot` for this application is `html/` alone (see
[`apache/planetgen.conf.example`](apache/planetgen.conf.example)), so
anything placed at the repo root, one level above `html/`, can never be
requested over HTTP no matter how the vhost or `.htaccess` rules are
written -- there's no path traversal or misconfiguration that reaches it,
because it's outside the tree Apache serves at all. The same reasoning is
why `db/`'s `.db` files (which can contain a full generated galaxy) live
next to `html/` rather than under it, and why `stellarObjects/webconfig.py`
resolves the repo root the same way
[`html/lib/dbutil.py`](html/lib/dbutil.py)'s `_PROJECT_ROOT`/
`DEFAULT_DB_DIR` already do (walking up from the module's own file via
`os.path.dirname`), rather than introducing a second convention for
locating "the project root."

This matters more for `webconfig.json` than it might for a settings file
with no sensitive contents, because the field list below includes
placeholders for future database credentials (see "Fields" below) -- a file
that might one day hold a password should never be reachable from inside
the served webroot in the first place, regardless of what Apache's own
`<Directory>` rules say.

## Fields

| Field | Purpose |
|---|---|
| `site_name` | Display name for this deployment, e.g. `"planetGen"`. Not yet read by any page -- reserved for a future UI element (page title, header) once the web interface wants to show it. |
| `base_url` | The base URL this deployment is served from, e.g. `"http://localhost/"` or `"https://planetgen.example.com/"`. Not yet read by any page -- reserved for future absolute-URL generation (e.g. constructing shareable links) that can't be derived from a CGI request alone. |
| `db_username` | **Unused placeholder.** Reserved for a possible future database backend other than SQLite. SQLite needs no credentials today, and nothing in this repo reads this field. |
| `db_password` | **Unused placeholder.** Same reasoning as `db_username` -- present only as scaffolding for a possible future non-SQLite backend. |
| `db_name` | **Unused placeholder.** Same reasoning as `db_username`/`db_password`. |

The three `db_*` fields are intentionally inert. They exist so that *if*
this project ever moves from SQLite (file-based, no credentials, located
via `PLANETGEN_DB_DIR`) to a client/server database, the config file
shape to hold its connection details already exists and doesn't require a
breaking change to `webconfig.json`'s structure. Nothing currently
generates, validates, or consumes them -- do not wire a feature to them
without first deciding on (and documenting here) which database backend
they'd actually apply to.

## Why the real file is gitignored but the example isn't

`webconfig.json` (the real, per-deployment file) is listed in
[`.gitignore`](.gitignore), next to the existing `*.db`/`db.sqlite3`
entries -- it's deployment-specific configuration, potentially holding
credentials once the `db_*` fields are ever put to use, and committing it
would either leak those values or force every deployment to share one
repo-tracked file. `webconfig.json.example` is the opposite: a template
with placeholder values, meant to be committed so a fresh checkout always
has a shape to copy from, exactly the same split already used for
`apache/planetgen.conf.example` (committed template, edited into a
site-specific vhost file that itself doesn't live in the repo).

## Setup

Copy the template and edit it for this deployment:

```bash
cp webconfig.json.example webconfig.json
```

Then edit `webconfig.json`'s `site_name`/`base_url` to match this
deployment. Leave the `db_*` fields as empty strings -- they're not read by
anything yet (see "Fields" above).

If `webconfig.json` doesn't exist at all, `stellarObjects.webconfig.load_webconfig()`
falls back to built-in defaults matching `webconfig.json.example`'s shape,
so nothing currently breaks by skipping this step -- it only matters once a
future feature actually reads `site_name`/`base_url` and a deployment wants
something other than the defaults.

## Relationship to `PLANETGEN_DB_DIR`/`PLANETGEN_DEBUG`

`webconfig.json` and the `PLANETGEN_DB_DIR`/`PLANETGEN_DEBUG` environment
variables (documented in [`html/README.md`](html/README.md)) solve
different problems and aren't interchangeable:

- The environment variables are **Apache-level overrides**, set via
  `SetEnv` in the vhost config -- they exist because CGI scripts only ever
  see what the web server's process environment hands them, and because
  `PLANETGEN_DB_DIR` in particular needs to vary per-vhost (e.g. multiple
  sites on one server, each pointed at a different `db/` directory) without
  editing any file inside the repo checkout itself.
- `webconfig.json` is a **site-level settings file**, edited once per
  deployment and read directly by application code (once something reads
  it) rather than passed through Apache at all -- it's the right place for
  settings that describe the site itself (its name, its own base URL)
  rather than settings that describe how this particular Apache instance
  is wired up to run the CGI scripts.

In short: environment variables configure *how Apache runs the scripts*;
`webconfig.json` configures *the site those scripts present*.
