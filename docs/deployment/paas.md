# VPS and platform-as-a-service hosting

## A VPS (recommended)

A small virtual server from any provider, running Debian or Ubuntu, is
the simplest place to host planetGen. Follow one of the Linux guides
([Apache](apache.md), [nginx](nginx.md) or [Caddy](caddy.md)) exactly as
on your own hardware. MySQL or MariaDB can run on the same server, or you
can use the provider's managed MySQL (MySQL 8.0.16+ or MariaDB 10.4+).

The app runs as one process with five threads. Generating sectors from
the admin Generate page is the heaviest work it does; the Galaxy Map's
tile cache on disk is capped by `tile_cache.max_mb` (default 200 MB).

Open ports 80 and 443 (and SSH) only. MySQL should listen on localhost, or on a
private network if it runs elsewhere.

## Platform-as-a-service

planetGen needs three things a typical platform-as-a-service makes
awkward:

- MySQL or MariaDB;
- a persistent disk for the Generate jobs directory (and, less
  critically, the tile cache, which can be rebuilt);
- a job runner that keeps running after the web request that started it.

As checked in September 2026:

| Platform | Persistent disk | MySQL | Notes |
|---|---|---|---|
| **Render** | Paid services only; one instance per disk; deploys stop the old instance first (brief downtime) | No managed MySQL (managed Postgres only) | You end up running MySQL yourself anyway |
| **Fly.io** | Volumes: one per Machine, on one host in one region, not replicated, daily snapshots | No managed MySQL | Works with a single Machine; backups are yours |
| **Railway** | One volume per service, no replicas with a volume, brief downtime on redeploy, size limited by plan | A MySQL template | The closest to one-click, within the plan's disk size |

planetGen runs as one process by design, so the scaling these platforms
sell doesn't help, and each needs its own build and start configuration
that this project doesn't ship. A VPS with a Linux guide is simpler and
usually cheaper.

If you use one anyway:

- Start the app with gunicorn: `gunicorn --pythonpath src/html
  --workers 1 --worker-class gthread --threads 5 --bind 0.0.0.0:$PORT
  wsgi:application`.
- Configure through `PLANETGEN_*` environment variables
  ([`config.md`](../config.md)) rather than a committed `config.json`.
  Set `PLANETGEN_SECRET_KEY`, the `PLANETGEN_MYSQL_*` variables, and
  `PLANETGEN_JOBS_DIR` and `PLANETGEN_TILE_CACHE_DIR` on the persistent
  disk.
- The platform's router is a reverse proxy. Set
  `PLANETGEN_PROXY_FIX_X_FOR` and `PLANETGEN_PROXY_FIX_X_PROTO` to the
  number of proxies the platform documents in front of your app (usually
  1), or every visitor shares one rate-limit budget
  ([`README.md`](README.md#behind-a-reverse-proxy-proxy_fix)).
- Run `python3 -m `planetgen.cli.migrate` once per deploy (a release or pre-deploy
  command) and read the first admin password from its output.
- Serve over HTTPS; the admin pages need it.
