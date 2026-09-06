# Apache2 Deployment Files

Deployment tooling for the web interface in [`../html/`](../html/README.md) --
not itself part of the served site. Nothing in this directory is copied
to the web server's document root.

| File | Purpose |
|---|---|
| `planetgen.conf.example` | Example Apache2 virtual host: `DocumentRoot` at `/var/lib/planetGen/html`, CGI handling for `html/*.py`, and access rules (denying direct requests to `html/lib/`). Copy to `/etc/apache2/sites-available/`, edit, and enable with `a2ensite`. |
| `set-permissions.sh` | Detects the user/group Apache2 is actually configured to run as (from `/etc/apache2/envvars`, or a running `apache2` process as a fallback) and sets ownership/permissions on the deployed `html/` and `db/` directories accordingly. Bash, not Python -- Linux-only deployment step, run once (and again after any redeploy) as root. |

See [`../html/README.md`](../html/README.md#deploying) for the full
deployment walkthrough.
