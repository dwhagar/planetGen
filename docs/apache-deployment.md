# Apache2 Deployment Files

Deployment tooling for the web interface in [`../src/html/`](html-interface.md) --
not itself part of the served site. The files it describes live in
[`../examples/apache/`](../examples/apache/) (grouped there as example/
template deployment config, alongside the example system files under
`../examples/systems/`); nothing in that directory is copied to the web
server's document root.

Most deployments just need `sudo ../install.sh` (first-time setup) or
`sudo ../update.sh` (pulling a later update) from the repo root -- both
call `set-permissions.sh` and `create-cache-dir.sh` (below) as steps. Use
`set-permissions.sh` directly only to re-apply permissions on its own
(e.g. after manually editing files as root, which resets ownership).

| File | Purpose |
|---|---|
| [`planetgen.conf.example`](../examples/apache/planetgen.conf.example) | Example Apache2 virtual host: `DocumentRoot` at `/var/lib/planetGen/src/html`, CGI handling for `../src/html/*.py`, a `WSGIScriptAlias` mounting the Flask app (the API under `/api` plus the HTML pages that have moved off CGI, via `../src/html/wsgi.py`, see [`api.md`](api.md#deploying-behind-apache-mod_wsgi) and [`html-interface.md`](html-interface.md#flask-pages-the-pages-are-moving-off-cgi)) at `/`, with `Alias /static/` and a `ScriptAliasMatch` for the remaining CGI pages, in its own `WSGIDaemonProcess` (real isolation from the CGI browser's own requests, not Apache's embedded/shared-process mode), and access rules (denying direct requests to `../src/html/lib/` and `../src/html/api/` -- the latter is imported by `mod_wsgi`, never meant to be requested/executed directly as CGI). Copy to `/etc/apache2/sites-available/`, edit, and enable with `a2ensite` -- the one step `install.sh` deliberately leaves manual. |
| [`set-permissions.sh`](../examples/apache/set-permissions.sh) | Detects the user/group Apache2 is actually configured to run as (from `/etc/apache2/envvars`, a running `apache2` process, or falling back to the Debian/Ubuntu default `www-data:www-data` if neither is found -- e.g. because apache2 isn't started yet), `chown`s the deployed `../src/html/` directory to that user:group, and makes every `*.py` file anywhere under it (at any subdirectory depth, `../src/html/lib/` included) executable. Prints a count of how many `.py` files it fixed, so a wrong path is obvious rather than silently matching nothing. Its optional `db-dir` argument is a pre-MySQL-port leftover (harmless no-op if that directory doesn't exist -- the database is a MySQL server now, not local files; see [`database-schema.md`](database-schema.md)). Bash, not Python -- Linux-only deployment step, safe to re-run any time as root. |
| [`create-cache-dir.sh`](../examples/apache/create-cache-dir.sh) | Creates the 3D Galaxy Map's on-disk tile cache (see [`config.md`](config.md)'s `tile_cache`: `PLANETGEN_TILE_CACHE_DIR`, else `tile_cache.dir`, else `/var/cache/planetgen/tiles`) and `chown`s it to Apache's user, detected the same way `set-permissions.sh` does (both source [`apache-identity.sh`](../examples/apache/apache-identity.sh)). Does nothing when `tile_cache.max_mb` is `0`. Takes the directory as an argument when it's set only by a `SetEnv` in the vhost, which the script can't see. Safe to re-run any time as root. |

See [`../install.sh`](../install.sh) for the full install script,
[`../update.sh`](../update.sh) for pulling later updates, and
[`html-interface.md`](html-interface.md#deploying) for the deployment
walkthrough.

## Static files, compression and security headers

- **Modules.** The example vhost needs `a2enmod cgid wsgi headers
  deflate`. `install.sh` enables `cgid headers deflate` itself; `deflate`
  is new (it is on by default in a stock Debian/Ubuntu Apache, so this is
  usually a no-op). Run `sudo a2enmod headers deflate && sudo systemctl
  reload apache2` once on an existing server.
- **Caching `static/`.** Every page links its CSS/JS/favicon as
  `static/<file>?v=<release version>` (`src/html/lib/fmt.py`'s
  `static_url`), and the map modules pass that same `?v=` on to the
  modules they import. The vhost's `<Directory .../static>` block sends
  `Cache-Control: public, max-age=31536000, immutable` for a request that
  carries `?v=` and `Cache-Control: no-cache` (always revalidate) for one
  that doesn't, so a release is picked up on the next page view without
  anyone clearing their cache.
- **Compression.** `AddOutputFilterByType DEFLATE` covers HTML, CSS,
  JavaScript, JSON and SVG. `three.module.min.js` goes from about 740 KB
  to about 190 KB on the wire.
- **Security headers have one source per response.** The HTML pages get
  `Content-Security-Policy`, `X-Frame-Options`, `X-Content-Type-Options`
  and `Referrer-Policy` from `src/html/lib/page.py` (`SECURITY_HEADERS`),
  and the API gets its own from `src/html/api/app.py`, so both work
  without this vhost. The page CSP is `default-src 'self'; base-uri
  'self'; form-action 'self'; frame-ancestors 'none'; object-src 'none'`.
  Earlier versions of the example vhost also set three of these
  site-wide with `Header always set ...`; **delete those three lines from
  an existing `/etc/apache2/sites-available/planetgen.conf`** when you
  update it. They only duplicated the application's own values. The one
  header Apache still adds is `X-Content-Type-Options: nosniff` on
  `static/`, which no application code serves.

## Routing: Flask pages, CGI pages and static files

The pages are moving from CGI to the Flask app one at a time (see
[`html-interface.md`](html-interface.md#flask-pages-the-pages-are-moving-off-cgi)),
so the vhost routes each URL to one of three places:

| URL | Served by | Directive |
|---|---|---|
| `/static/...` | Apache, straight from `src/html/static/` | `Alias /static/ .../src/html/static/` |
| `/<name>.py` (a page still on CGI) | `mod_cgid` | `ScriptAliasMatch "^/((?!wsgi\.py$)[a-z_]+\.py)$" ".../src/html/$1"` |
| everything else: `/`, `/sectors`, `/systems`, `/api/...`, each page as it moves | the Flask app in the `planetgen-api` daemon | `WSGIScriptAlias / .../src/html/wsgi.py` |

mod_wsgi applies `WSGIScriptAlias` after mod_alias's `Alias`/
`ScriptAliasMatch`, so those two win for the paths they match even
though `WSGIScriptAlias /` matches every URL. The regex matches only a
top-level `name.py` made of lowercase letters and underscores, never
`wsgi.py` and never anything under `lib/` or `api/` (those reach Flask
and get its 404 page). `index.py` and `browse.py` are now 301 shims to
`/`, so old links and bookmarks keep working. The CGI pages still call
the API over HTTP at `PLANETGEN_API_BASE_URL` (default
`http://127.0.0.1/api`), which the `/` mount serves as before.

Tested with Apache 2.4.58 + mod_wsgi (Ubuntu 24.04): `apache2ctl
configtest`, then real requests for each row above, a CGI login followed
by a Flask page showing the admin menu, a POST to `browse.py` (301, then
GET `/`), and the HTTPS section's `:80` redirect with the `/` mount (only
`/api` reaches the app).

### Updating an existing server

1. `sudo ./update.sh` (pulls, re-runs `install.sh`).
2. Add a `secret_key` to `config.json` (it signs the Flask pages' form
   tokens; see [`config.md`](config.md)):
   `python3 -c "import secrets; print(secrets.token_hex(32))"`.
3. Check `config.json`'s `mysql.database` (or `PLANETGEN_MYSQL_DATABASE`
   in Apache's environment) names the database the site should show. The
   old `index.py` showed the alphabetically first database matching
   `database_prefix`; the Flask pages show this one.
4. Edit `/etc/apache2/sites-available/planetgen.conf` (and the `:443`
   block, if you have one) to match the new example:
   - add `Alias /static/ /var/lib/planetGen/src/html/static/` and the
     `ScriptAliasMatch` line above, before the `WSGIScriptAlias`;
   - change `WSGIScriptAlias /api /var/lib/planetGen/src/html/wsgi.py ...`
     to `WSGIScriptAlias / /var/lib/planetGen/src/html/wsgi.py ...`
     (same options; don't keep both);
   - remove `DirectoryIndex index.py` from `<Directory .../src/html>`.
5. `sudo apache2ctl configtest && sudo systemctl reload apache2`.
6. Check: `/` shows the new header, `/index.py` answers 301 to `/`,
   `/sector.py?...` still works, `/static/style.css` is served, and
   `curl http://127.0.0.1/api/health` still answers.

## MySQL accounts

One MySQL account (`PLANETGEN_MYSQL_*`, or `config.json`'s `mysql`
section -- see [`config.md`](config.md) for setting it once in a shared
`config.json` instead of repeating it across every vhost/service file) is
used everywhere: the generation CLIs, `install.sh`/`migrateDb.py` (which
need full `CREATE`/`ALTER`/DML grants -- that's the account schema DDL
runs against), and the deployed Flask API's reads *and* writes alike
(sector/system create/update/delete, and everything under `/api/auth/`,
need at least `INSERT`/`UPDATE`/`DELETE`/`SELECT`). There's no separate
write-capable account to configure -- `WRITE_MYSQL_CONFIG` simply reuses
`MYSQL_CONFIG` (see `html/api/config.py`). If a deployment still wants the
*deployed* Apache/`mod_wsgi` process on a different, less-privileged
account than the one `install.sh`/`migrateDb.py` run as, point its own
`PLANETGEN_MYSQL_*` at that account separately (see the vhost example's
comments) -- but remember it now needs write grants too if the API's
write endpoints are reachable, not just `SELECT`.

That same account also needs to create and seed the **control schema**
(`PLANETGEN_CONTROL_DATABASE`, default `planetgen_control`, see
[`database-schema.md`](database-schema.md#the-control-schema)) --
`migrateDb.py` does this automatically, alongside its usual
content-schema migration, every time it's run (i.e. every
`install.sh`/`update.sh`).

**The admin web UI (`/login`, `/admin`, `/admin/stats`, `/account`) requires
HTTPS** — its session cookie is `Secure` by default and simply won't be
sent by the browser over plain HTTP. Terminate TLS in front of this vhost
(e.g. `certbot --apache`) before relying on it; see
[`api.md`](api.md#deploying-behind-apache-mod_wsgi)'s note on
`PLANETGEN_ADMIN_COOKIE_INSECURE` (or `config.json`'s
`admin_cookie_insecure`, see [`config.md`](config.md)) for the
local-development-only escape hatch.
