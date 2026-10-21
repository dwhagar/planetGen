# Deploying planetGen

planetGen's web interface is one Flask app: every page and the JSON API
(`/api`). It runs from a git checkout, reads `config.json` at the repo
root, and talks to one MySQL server and one Redis server
(`redis.url` in [`config.json`](../config.md); the Generate page's jobs,
the work queue, the rate limits and the login lockouts run on it). This page compares the
supported ways to host it and covers what they share. Each guide has the
full steps for its platform.

## Choose a platform

| Platform | Guide | Web server | App server | Install scripts | Admin Generate page | HTTPS | Good for |
|---|---|---|---|---|---|---|---|
| Linux (Debian, Ubuntu) | [Apache](apache.md) | Apache2 | mod_wsgi (daemon mode) | `install.sh`, `update.sh` do everything but the vhost | Works | certbot `--apache` | The reference setup; least manual work |
| Linux | [nginx](nginx.md) | nginx | gunicorn under systemd | `install.sh`, `update.sh`, plus a systemd unit and a site file | Works | certbot `--nginx` | Servers that already run nginx |
| Linux | [Caddy](caddy.md) | Caddy | gunicorn under systemd | As for nginx | Works | Automatic | The shortest web server config |
| macOS | [macOS](macos.md) | Homebrew nginx | gunicorn under launchd | `install.sh` | Works | certbot or your own certificate | A Mac you already have. macOS Server is discontinued |
| A rented server | [VPS and PaaS](paas.md) | Any Linux option | Any Linux option | As on Linux | Works | As on Linux | Hosting you don't run at home |

Whatever the platform, on Linux add [fail2ban](fail2ban.md) to ban
addresses that keep guessing admin passwords.

## What every setup has in common

- **One process, five threads.** Apache's `WSGIDaemonProcess
  processes=1 threads=5`, gunicorn `--workers 1 --threads 5`. The rate limits and login lockouts count on
  Redis (`ratelimit.storage_uri` empty = the server in `redis.url`), so
  more processes share them; `memory://` is only right with one process
  ([`api.md`](../api.md#rate-limiting)).
- **Threads, because of the live job log.** The admin Generate page streams
  a running job's output over Server-Sent Events (`/admin/generate/jobs/<id>/stream`,
  ADM.22). Each open stream holds one of the five threads, so the app
  server must be threaded (the settings above are; a single-threaded
  worker would be held by one viewer). A response ends after about 40
  seconds, under mod_wsgi's `request-timeout=60`, and the browser
  reconnects by itself, resuming where it stopped. A proxy in front must
  not buffer the stream (the response sends `X-Accel-Buffering: no` for
  nginx).
- **The same entry point.** Every app server loads `src/html/wsgi.py`
  and its `application` object. gunicorn, from planetGen's own virtual environment
  (`/usr/local/planetgen/venv`, [apache.md](apache.md#planetgens-own-python)):
  `/usr/local/planetgen/venv/bin/gunicorn --pythonpath
  <checkout>/src/html wsgi:application`.
- **Static files from the web server.** `/static/` maps to
  `src/html/static/`, with `X-Content-Type-Options: nosniff`. A URL with
  `?v=` (every page adds the release version) is cached for a year
  (`public, max-age=31536000, immutable`); one without is `no-cache`.
  Nothing else in the checkout is served from disk. Every other URL goes
  to the app, including the old `/<name>.py` page URLs, which the app
  redirects.
- **Security headers from the app.** The app sets
  `Content-Security-Policy`, `X-Frame-Options`, `X-Content-Type-Options`,
  `Referrer-Policy`, and `Strict-Transport-Security` on HTTPS requests.
  Don't repeat them in the web server.
- **HTTPS.** The admin session cookie is `Secure`, so `/login` and the
  admin pages only work over HTTPS.
- **A 60-second limit** on a request, compression for text responses,
  and a 1 MB request body limit (the API only takes small JSON bodies).
- **Runtime directories** owned by the account the app runs as, mode
  750: the Galaxy Map tile cache (`tile_cache.dir`, default
  `/var/cache/planetgen/tiles`) and the Generate jobs (`jobs.dir`,
  default `/var/lib/planetGen/jobs`). On Linux, `install.sh` and
  `update.sh` create both (`examples/apache/create-cache-dir.sh`), and
  move jobs left in the old default, `/var/lib/planetgen/jobs`, into the
  new one (a running job stays until it finishes).
- **The checkout belongs to the web user.** Everything under the
  checkout (`/var/lib/planetGen`), including `config.json` and the
  `state/` directory the app writes its maintenance flag to, is owned by
  the account the app runs as (`www-data:www-data` on Debian and Ubuntu):
  directories 750, files 640, and scripts and `*.py` files 750. Nobody
  else can read it. A fresh `install.sh` and an `update.sh` end with the
  same ownership because both run `set-permissions.sh` on the whole
  checkout. The virtual environment (`/usr/local/planetgen/venv`) and the
  NLTK corpus stay root's, since `install.sh` and `update.sh` run their
  Python as root. The web server is kept out of the checkout by the
  directory rules in the example Apache site: nothing is served from disk
  but `static/`, and nothing runs but the WSGI app. The debug log, when
  `debug` is on, belongs to the app's account, mode 0660. On Linux
  `install.sh` and `update.sh` set all of this (`set-permissions.sh`,
  `setup-debug-log.sh`); the macOS guide gives the manual equivalent.
- **Redis settings.** Redis holds the work queue, rate limits and
  login lockouts, and out of the box keeps nothing on disk between
  saves, so a crash loses all of it. [`examples/redis/planetgen.conf`](../../examples/redis/planetgen.conf)
  turns on the append-only file (`appendonly yes`, `appendfsync
  everysec`), `maxmemory-policy noeviction`, a `maxmemory 512mb` guard and
  a loopback-only `bind`. Add those lines to your Redis config (on Debian
  and Ubuntu `/etc/redis/redis.conf`, or end it with `include
  /path/to/planetgen.conf`) and restart Redis. `install.sh` and
  `update.sh` check the running server and, for each setting that is off,
  print the exact line to add; they never edit your Redis config. One
  Redis serves everything. It needs no password while it binds loopback;
  the example shows the ACL user to add later.
- **The NLTK `words` corpus** somewhere the app's account can read:
  `/usr/local/share/nltk_data` on Linux and macOS, or a folder named by
  `NLTK_DATA`.
- **An optional population pass.** After the database step, the install
  and update scripts ask whether to run `planetgen population`
  (species, civilizations, territories): y/N, default N after 30
  seconds, and skipped when there is no terminal or console to ask on.
  `POPULATION=1` runs it without asking. It
  can be run by hand any time.

## Behind a reverse proxy: `proxy_fix`

With Apache and mod_wsgi, Apache owns the client's connection, so the app
sees the client's real address and whether the request came over HTTPS.
With nginx, Caddy or Apache's `mod_proxy` in front of gunicorn,
the app only sees the proxy. Then:

- every visitor has the proxy's address, so they all share one
  rate-limit budget, and a few visitors can get everyone a 429;
- every request looks like plain HTTP, so the app never sends
  `Strict-Transport-Security`;
- the debug log records the proxy's address for every request.

`config.json`'s `proxy_fix` fixes this. Each field is the number of
proxies to trust for one header:

```json
"proxy_fix": {"x_for": 1, "x_proto": 1, "x_host": 0}
```

`x_for` takes the client address from `X-Forwarded-For`, `x_proto` the
scheme from `X-Forwarded-Proto`, `x_host` the host from
`X-Forwarded-Host`. All three are 0 by default, which changes nothing.
The environment variables `PLANETGEN_PROXY_FIX_X_FOR`, `_X_PROTO` and
`_X_HOST` override the file. See [`config.md`](../config.md).

Rules:

- Use 1 for one proxy on the same machine, which is every guide here.
  Leave `x_host` at 0: every example passes the original `Host` header
  through.
- Only turn it on when the app server can be reached **only** through
  that proxy: a Unix socket, or `127.0.0.1`. Anyone who can reach the app
  server directly can set these headers and choose their own address.
- The proxy must set or append the headers. nginx, Caddy and Apache's
  `mod_proxy` do, as configured in the examples.
- Leave it at 0 under Apache + mod_wsgi.

Check it: request a page over HTTPS and look for
`Strict-Transport-Security` in the response. With `debug` on, the debug
log's `API request: ... from <address>` lines should show the client's
address, not `127.0.0.1`.

## MySQL accounts

planetGen needs MySQL 8.0.19 or later (MariaDB is not supported). One
account is used everywhere: the generation CLIs, `planetgen.cli.migrate` (run by
`install.sh` and `update.sh`), and the web app's reads and writes. There
is no separate read-only or write account to configure:
`WRITE_MYSQL_CONFIG` and `CONTROL_MYSQL_CONFIG` reuse `MYSQL_CONFIG`
(`src/planetgen/api/config.py`). Set it once in `config.json`'s `mysql`
section (or `PLANETGEN_MYSQL_*`; see [`config.md`](../config.md)).

That account creates and changes the schema, so it needs full rights on
the game database (`mysql.database`, default `planetgen`), on the
**control schema** (`control_database`, default `planetgen_control`,
which `planetgen.cli.migrate` creates: admin logins, sessions, API keys, audit
log; see [`database-schema.md`](../database-schema.md#the-control-schema)),
and on any other game databases that share the prefix. For example:

```sql
CREATE DATABASE planetgen CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci;
CREATE USER 'planetgen'@'localhost' IDENTIFIED BY 'a long random password';
GRANT ALL PRIVILEGES ON `planetgen%`.* TO 'planetgen'@'localhost';
```

Use `'planetgen'@'%'` (or the web host's address) when MySQL runs on
another machine. If you want the web app on a less privileged account
than the one `install.sh` runs as, give the web app's process its own
`PLANETGEN_MYSQL_USER`/`_PASSWORD`. That account still needs `SELECT`,
`INSERT`, `UPDATE` and `DELETE` on the game and control schemas, because
the admin pages and the API's write endpoints write through it.

## The first admin login

The first `planetgen.cli.migrate` run against an empty control schema creates the
admin `admin` with a random password and prints it once, in a boxed
block (on Linux, in `install.sh`'s step 2/8). Only its hash is stored.
Log in at `https://<your site>/login`; you are sent to `/account` to
choose your own username and password, and the other admin pages stay
locked until you do. If the password is lost, see
[`api.md`](../api.md#resetting-the-admin-login).

## Updating a running site

`sudo ./update.sh` takes the **update guard** (`planetgen.cli.update_guard`, `planetgen/admin/updateguard.py`)
before it changes anything, so the only thing that runs is the update:

1. **It lists the running and queued work** (a Generate job, a `planetgen` run from the command line, jobs on the
   shared queues, an orbital update, a Galaxy Map warm-up from an earlier update) and **stops, having changed and
   cancelled nothing, unless you agree to cancel it.** On a terminal it asks (y/N). With no terminal (the
   maintenance timer, cron) it refuses with exit status 1 and names the work; add `--cancel-running` to cancel it:
   `sudo ./update.sh --cancel-running`. Each piece is stopped through its own cancel path and the update waits (up
   to 10 minutes) until none is left; one that does not stop is reported and the update stops before any change.
2. **It closes the site**: maintenance mode on with the cause "update" and a lock, `update.lock`, in the state
   folder (`/var/lib/planetGen/state`). Visitors get the "being updated" page; admins can still sign in.
3. **Nothing new starts** until the update ends: the Generate page and the API refuse to queue work ("the site is
   being updated"), the queue's workers exit with status 75, and the orbital update exits with 75 without moving
   anything (its timer simply tries again at its next firing).
4. **A good update opens the site again** at its end (a closing an admin made earlier stays; the conversion to the
   spherical sector grid keeps the site closed until the new Galaxy Map is built). **A failed update leaves the
   site closed**, says which step failed and how to open it, because the database may be half migrated:
   `python3 -m planetgen.cli.maintenance off`. The lock is released, so `./update.sh` can be run again.

`python3 -m planetgen.cli.migrate` run by hand does the same when it has something to migrate (`--cancel-running`
works there too). `update.sh` never edits your Redis configuration or your Apache site.

## After deploying

Run the [server checklist](../server-checklist.md). It has the commands
for each platform.
