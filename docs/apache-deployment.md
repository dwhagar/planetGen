# Apache2 Deployment Files

Deployment tooling for the web interface in [`../src/html/`](html-interface.md) --
not itself part of the served site. The files it describes live in
[`../examples/apache/`](../examples/apache/) (grouped there as example/
template deployment config, alongside the example system files under
`../examples/systems/`); nothing in that directory is copied to the web
server's document root.

Most deployments just need `sudo ../install.sh` (first-time setup) or
`sudo ../update.sh` (pulling a later update) from the repo root -- both
call `set-permissions.sh` (below) as one of their steps. Use
`set-permissions.sh` directly only to re-apply permissions on its own
(e.g. after manually editing files as root, which resets ownership).

| File | Purpose |
|---|---|
| [`planetgen.conf.example`](../examples/apache/planetgen.conf.example) | Example Apache2 virtual host: `DocumentRoot` at `/var/lib/planetGen/src/html`, CGI handling for `../src/html/*.py`, a `WSGIScriptAlias` mounting the Flask API (`../src/html/api/`, via `../src/html/wsgi.py`, see [`api.md`](api.md#deploying-behind-apache-mod_wsgi)) at `/api` in its own `WSGIDaemonProcess` (real isolation from the CGI browser's own requests, not Apache's embedded/shared-process mode), and access rules (denying direct requests to `../src/html/lib/` and `../src/html/api/` -- the latter is imported by `mod_wsgi`, never meant to be requested/executed directly as CGI). Copy to `/etc/apache2/sites-available/`, edit, and enable with `a2ensite` -- the one step `install.sh` deliberately leaves manual. |
| [`set-permissions.sh`](../examples/apache/set-permissions.sh) | Detects the user/group Apache2 is actually configured to run as (from `/etc/apache2/envvars`, a running `apache2` process, or falling back to the Debian/Ubuntu default `www-data:www-data` if neither is found -- e.g. because apache2 isn't started yet), `chown`s the deployed `../src/html/` directory to that user:group, and makes every `*.py` file anywhere under it (at any subdirectory depth, `../src/html/lib/` included) executable. Prints a count of how many `.py` files it fixed, so a wrong path is obvious rather than silently matching nothing. Its optional `db-dir` argument is a pre-MySQL-port leftover (harmless no-op if that directory doesn't exist -- the database is a MySQL server now, not local files; see [`database-schema.md`](database-schema.md)). Bash, not Python -- Linux-only deployment step, safe to re-run any time as root. |

See [`../install.sh`](../install.sh) for the full install script,
[`../update.sh`](../update.sh) for pulling later updates, and
[`html-interface.md`](html-interface.md#deploying) for the deployment
walkthrough.
