# planetGen Web Interface

A small, dependency-free web interface for browsing the SQLite databases
described in [`db/README.md`](../db/README.md) -- pick a `.db` file,
drill into its sectors and star systems, and view (or copy) the
rendered wikitext/Markdown page saved for each one.

This is a read-only browser, not the Phase 5 web application described in
[`TODO.md`](../TODO.md) -- there's no backend framework, no search/filter
UI, and no writing to the database. It exists so a generated galaxy can be
looked at from a browser today, on nothing more than Apache2 and a system
Python 3.

## How it works

Every page here is a plain [CGI](https://en.wikipedia.org/wiki/Common_Gateway_Interface)
script (`#!/usr/bin/env python3`, executed fresh by Apache on every
request) using only the standard library -- no Flask/FastAPI, no `pip
install` needed beyond the `planetGen` package itself. That keeps
deployment to "copy the files, set permissions, enable a vhost" with
nothing else to install or run.

| File | Purpose |
|---|---|
| `index.py` | Lists every `.db` file in the database directory. |
| `browse.py` | A chosen database's sectors and standalone systems. |
| `sector.py` | One sector's systems (name, quadrant, star type). |
| `system.py` | One system's stars/planets/moons/belts, and the full rendered wikitext or Markdown page (toggle via a link), shown in a `<textarea>` for easy copy-paste into the wiki. |
| `lib/dbutil.py` | Read-only database access and HTML-escaping helpers. Not web-accessible -- see the Apache config note below. |
| `lib/page.py` | Shared CGI response/HTML-shell helpers. Not web-accessible. |
| `static/style.css` | Shared stylesheet, served directly (not through CGI). |

All database access goes through `sqlite3`'s `file:...?mode=ro` URI mode,
so these scripts cannot write to a database even if a query were buggy.
Database and system names pulled from generated data are HTML-escaped
before being placed in a page; a requested `?db=` filename is validated
against the actual directory listing (exact basename match only), which
is what prevents it from being used for path traversal.

## Locating the database directory

By default, each script looks for `db/` as a sibling of `html/` --
matching this repo's own layout, and the recommended deployment layout
(`/var/lib/planetGen/html` and `/var/lib/planetGen/db` side by side; see
`../apache/planetgen.conf.example`). Set the `PLANETGEN_DB_DIR`
environment variable (e.g. via `SetEnv` in the Apache vhost) to point
somewhere else.

## Deploying

1. Copy the repo (or at least `html/`, `db/`, `stellarObjects/`,
   `install.sh`, `setup.py`, and `apache/`) to the server, e.g.
   `/var/lib/planetGen/`.
2. From that directory, run `sudo ./install.sh` -- installs the Python
   package, pre-fetches the NLTK `words` corpus into a shared
   world-readable location (so it works under Apache's `www-data`, not
   just whatever user happens to run the CLI tools), makes the CGI
   scripts executable, enables Apache's CGI module, and sets `html/`/`db/`
   ownership for Apache via `apache/set-permissions.sh`. See
   [`../apache/README.md`](../apache/README.md) for what each step does
   and how to re-run pieces of it individually.
3. `install.sh` prints one remaining manual step: copy
   `apache/planetgen.conf.example` to
   `/etc/apache2/sites-available/planetgen.conf`, edit it (at minimum,
   `ServerName`), then `sudo a2ensite planetgen && sudo systemctl reload
   apache2`. This is deliberately not automated -- the vhost's
   ServerName/TLS/logging are your call, not something to silently create
   or overwrite.

Re-run `sudo ./install.sh` any time after pulling updates -- every step
is idempotent (the corpus fetch skips itself if already present, `chmod
+x`/`a2enmod`/permission-setting are all safe to repeat).

## Local testing without Apache

Every script is a normal CGI program: it reads `QUERY_STRING` from the
environment and writes an HTTP response (status + headers + body) to
stdout. That makes them runnable directly for a quick smoke test without
standing up Apache at all:

```bash
PLANETGEN_DB_DIR=/path/to/db QUERY_STRING="db=planetgen.db" python3 browse.py
```
