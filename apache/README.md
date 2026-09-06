# Apache2 Deployment Files

Deployment tooling for the web interface in [`../html/`](../html/README.md) --
not itself part of the served site. Nothing in this directory is copied
to the web server's document root.

Most deployments just need `sudo ../install.sh` from the repo root --
it calls `set-permissions.sh` (below) as one of its steps. Use
`set-permissions.sh` directly only to re-apply permissions on its own
(e.g. after manually editing files as root, which resets ownership).

| File | Purpose |
|---|---|
| `planetgen.conf.example` | Example Apache2 virtual host: `DocumentRoot` at `/var/lib/planetGen/html`, CGI handling for `html/*.py`, and access rules (denying direct requests to `html/lib/`). Copy to `/etc/apache2/sites-available/`, edit, and enable with `a2ensite` -- the one step `install.sh` deliberately leaves manual. |
| `set-permissions.sh` | Detects the user/group Apache2 is actually configured to run as (from `/etc/apache2/envvars`, a running `apache2` process, or falling back to the Debian/Ubuntu default `www-data:www-data` if neither is found -- e.g. because apache2 isn't started yet) and sets ownership/permissions on the deployed `html/` and `db/` directories accordingly. Bash, not Python -- Linux-only deployment step, safe to re-run any time as root. |

See [`../install.sh`](../install.sh) for the full install script and
[`../html/README.md`](../html/README.md#deploying) for the deployment
walkthrough.
