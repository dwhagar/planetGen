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
| [`planetgen.conf.example`](../examples/apache/planetgen.conf.example) | Example Apache2 virtual host: `DocumentRoot` at `/var/lib/planetGen/src/html`, a `WSGIScriptAlias` mounting the Flask app (the API under `/api` plus every page, via `../src/html/wsgi.py`, see [`api.md`](api.md#deploying-behind-apache-mod_wsgi) and [`html-interface.md`](html-interface.md#flask-pages)) at `/`, with `Alias /static/` for the static files, in its own `WSGIDaemonProcess` (not Apache's embedded/shared-process mode), and access rules (denying direct requests to `../src/html/lib/` and `../src/html/api/`). Copy to `/etc/apache2/sites-available/`, edit, and enable with `a2ensite` -- the one step `install.sh` deliberately leaves manual. |
| [`set-permissions.sh`](../examples/apache/set-permissions.sh) | Detects the user/group Apache2 is actually configured to run as (from `/etc/apache2/envvars`, a running `apache2` process, or falling back to the Debian/Ubuntu default `www-data:www-data` if neither is found -- e.g. because apache2 isn't started yet), `chown`s the deployed `../src/html/` directory to that user:group, and makes every `*.py` file anywhere under it (at any subdirectory depth, `../src/html/lib/` included) owner/group readable (and executable, a harmless leftover from the CGI pages). Prints a count of how many `.py` files it fixed, so a wrong path is obvious rather than silently matching nothing. Its optional `db-dir` argument is a pre-MySQL-port leftover (harmless no-op if that directory doesn't exist -- the database is a MySQL server now, not local files; see [`database-schema.md`](database-schema.md)). Bash, not Python -- Linux-only deployment step, safe to re-run any time as root. |
| [`create-cache-dir.sh`](../examples/apache/create-cache-dir.sh) | Creates the 3D Galaxy Map's on-disk tile cache (see [`config.md`](config.md)'s `tile_cache`: `PLANETGEN_TILE_CACHE_DIR`, else `tile_cache.dir`, else `/var/cache/planetgen/tiles`) and `chown`s it to Apache's user, detected the same way `set-permissions.sh` does (both source [`apache-identity.sh`](../examples/apache/apache-identity.sh)). Does nothing when `tile_cache.max_mb` is `0`. Takes the directory as an argument when it's set only by a `SetEnv` in the vhost, which the script can't see. Safe to re-run any time as root. |

## Managed Python

Newer Debian and Ubuntu releases (Debian 12+, Ubuntu 23.04+, including
Ubuntu 26.04 LTS) mark their system Python as *externally managed*
(PEP 668: an `EXTERNALLY-MANAGED` file in the stdlib directory), and pip
refuses to install into it. `install.sh`'s first step
([`../scripts/install-python-deps.sh`](../scripts/install-python-deps.sh))
checks for that file and picks a path, printing which one it took:

- **Not managed:** `pip install` of the package and its `api` extra, as
  before.
- **Managed:** the libraries come from apt (`python3-flask`,
  `python3-nltk`, `python3-pymysql`, `python3-dbutils`,
  `python3-werkzeug`, `python3-rich`, `python3-flask-limiter`). Any the
  distribution doesn't package, or packages below the version `setup.py`
  asks for, are pip-installed into a venv at `/opt/planetgen/venv`
  (`PLANETGEN_VENV_DIR`), and a `planetgen-venv.pth` file in the system
  site-packages puts that venv first on the system Python's `sys.path` so
  mod_wsgi sees it without any vhost change. The
  `.pth` is removed again once apt provides everything. planetGen itself
  runs straight from the checkout (every entry point adds `src/` to
  `sys.path`), with a `/usr/local/bin/planetgen` wrapper standing in for
  pip's console script. (The wrapper is written on the unmanaged path
  too, so the CLI always runs the checkout's code.)

`PLANETGEN_PYTHON_MODE=managed` or `unmanaged` overrides the detection.

`update.sh` reinstalls nothing. It runs the same script with `--check`,
which imports each library with the system Python (the venv `.pth`
included) and compares its version with `setup.py`'s floor, then
installs only what is missing, too old or fails to import, the same way
`install.sh` would on that host: apt and then the venv (kept, not
rebuilt) on a managed Python, pip on an ordinary one. It prints one line
per library (`present`, `installed`, `upgraded`, `repaired` or `failed`)
and stops the update if anything is still unusable. Its last step imports
the web app as Apache's user, so a library www-data can't read shows up
there rather than as a 500. `sudo ./install.sh` is still the full
reinstall.

### Which Python and libraries Apache uses

Nothing in the vhost points at the venv, and nothing needs to. mod_wsgi
embeds the system Python it was built against, and that Python reads
`planetgen-venv.pth` from its own site-packages at startup, so the venv
comes first on `sys.path` in the `planetgen-api` daemon exactly as it
does for `python3` in a shell. Don't add `python-home=/opt/planetgen/venv`
to `WSGIDaemonProcess`: the venv holds only the libraries apt lacks, and
the apt ones can drop out of the path.

The one thing that has to line up is the Python version: mod_wsgi
(`libapache2-mod-wsgi-py3`) must be built for the same Python that
`install.sh`/`update.sh` set up. Both scripts check this and warn if
they differ. By hand:

    ldd /usr/lib/apache2/modules/mod_wsgi.so | grep libpython   # e.g. libpython3.12
    python3 --version                                           # must match
    sudo -u www-data python3 -c "import flask; print(flask.__file__)"

If they differ, install the `libapache2-mod-wsgi-py3` that matches, or run
the scripts with that Python (`sudo PYTHON=/usr/bin/python3.12
./update.sh`).

To see what the running site really uses, open the admin Stats page
(`/admin/stats`): **Python** shows the daemon's version and prefix, and
**Libraries from** shows the directory it imports Flask from
(`/opt/planetgen/venv/lib/python3.X/site-packages` when the venv
fallback is in use, `/usr/lib/python3/dist-packages` when apt provides
everything). After an update, `sudo systemctl reload apache2` restarts
the daemon so it picks up anything newly installed.

See [`../install.sh`](../install.sh) for the full install script,
[`../update.sh`](../update.sh) for pulling later updates, and
[`html-interface.md`](html-interface.md#deploying) for the deployment
walkthrough.

## Static files, compression and security headers

- **Modules.** The example vhost needs `a2enmod wsgi headers deflate`
  (no CGI module any more). `install.sh` enables all three itself (`wsgi`
  when `libapache2-mod-wsgi-py3` is installed; it warns otherwise).
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
  and `Referrer-Policy` from `src/html/web/__init__.py` (`SECURITY_HEADERS`),
  and the API gets its own from `src/html/api/app.py`, so both work
  without this vhost. The page CSP is `default-src 'self'; base-uri
  'self'; form-action 'self'; frame-ancestors 'none'; object-src 'none'`.
  Earlier versions of the example vhost also set three of these
  site-wide with `Header always set ...`; **delete those three lines from
  an existing `/etc/apache2/sites-available/planetgen.conf`** when you
  update it. They only duplicated the application's own values. The one
  header Apache still adds is `X-Content-Type-Options: nosniff` on
  `static/`, which no application code serves.

## Routing: the Flask app and static files

| URL | Served by | Directive |
|---|---|---|
| `/static/...` | Apache, straight from `src/html/static/` | `Alias /static/ .../src/html/static/` |
| everything else: `/`, `/sectors`, `/system/<id>`, `/api/...`, and the old `/<name>.py` page URLs | the Flask app in the `planetgen-api` daemon | `WSGIScriptAlias / .../src/html/wsgi.py` |

mod_wsgi applies `WSGIScriptAlias` after mod_alias's `Alias`, so
`/static/` wins for the paths it matches even though `WSGIScriptAlias /`
matches every URL. The old CGI page URLs (`/index.py`, `/sector.py?id=5`,
...) reach the app, which answers each with a 301 to the page that
replaced it (see
[`html-interface.md`](html-interface.md#flask-pages), "Old URLs"); any
other `.py` name gets the app's 404 page.

### Updating an existing server

1. `sudo ./update.sh` (pulls, then checks and fills in only what's missing).
2. Make sure `config.json` has a `secret_key` (it signs the pages' form
   tokens; see [`config.md`](config.md)):
   `python3 -c "import secrets; print(secrets.token_hex(32))"`, and that
   its `mysql.database` (or `PLANETGEN_MYSQL_DATABASE` in Apache's
   environment) names the database the site should show.
3. Edit `/etc/apache2/sites-available/planetgen.conf` (and the `:443`
   block, if you have one) to match the example:
   - remove the `ScriptAliasMatch "^/((?!wsgi\.py$)[a-z_]+\.py)$" ...`
     line (until you do, Apache answers old `.py` URLs with its own 404
     instead of letting the app redirect them);
   - in `<Directory /var/lib/planetGen/src/html>`, change `Options
     +ExecCGI -Indexes` to `Options -Indexes` and remove `AddHandler
     cgi-script .py`;
   - then remove the `<Files "wsgi.py"> SetHandler wsgi-script </Files>`
     block (it only overrode that `AddHandler`; keep it if you keep the
     `AddHandler` line);
   - in `<Directory .../src/html/static>`, `Options -ExecCGI -Indexes`
     can become `Options -Indexes`;
   - remove any `SetEnv PLANETGEN_*` lines (they only ever reached the
     CGI pages; the app reads `config.json` or Apache's own
     environment);
   - keep `Alias /static/ ...` and `WSGIScriptAlias / ...`; an older
     `WSGIScriptAlias /api ...` should be `WSGIScriptAlias / ...`.
4. `sudo apache2ctl configtest && sudo systemctl reload apache2`. The CGI
   module is no longer needed: `sudo a2dismod cgid` (then reload again)
   if nothing else on the server uses it.
5. Check: `/` shows the site, `/index.py` and `/sector.py?id=1` answer
   301 to `/` and `/sector/1`, `/static/style.css` is served, and `curl
   http://127.0.0.1/api/health` still answers.

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
`PLANETGEN_MYSQL_*` at that account separately in Apache's own
environment (see the vhost example's comments) -- but remember it now needs write grants too if the API's
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
