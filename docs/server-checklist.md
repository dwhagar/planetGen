# Server checklist: confirm the site loads

Run this after deploying or updating, on any platform. It started as the
check for the "site won't load" outages (request timeouts, oversized
responses, the web server being killed for memory after Galaxy Map
traffic), and covers the rest of a healthy install too.

Replace `HOST` with your site's hostname, `DB_NAME` with your database,
and `/var/lib/planetGen` with your checkout if it lives elsewhere
(`C:\srv\planetGen` in the Windows guide).

## Commands for your platform

The steps below name these by their row:

| | Apache + mod_wsgi | nginx or Caddy + gunicorn (Linux) | macOS (launchd) | Windows |
|---|---|---|---|---|
| **reload the app** | `sudo systemctl reload apache2` | `sudo systemctl reload planetgen-gunicorn` | `sudo launchctl kill SIGHUP system/org.planetgen.gunicorn` | IIS: `Restart-WebAppPool planetgen`; service: `Restart-Service planetgen` |
| **app error log** | `/var/log/apache2/planetgen_error.log` | `journalctl -u planetgen-gunicorn` | `/usr/local/planetgen/log/gunicorn.log` | `C:\ProgramData\planetgen\logs` |
| **app processes** | `ps -o rss,cmd -C apache2` | `ps -o rss,cmd -u www-data` | `ps -o rss,command -U _www` | Task Manager, `python.exe` / `waitress-serve.exe` |
| **web user** | `www-data` | `www-data` | `_www` | the app pool or service account |
| **request time limit** | `request-timeout=60` in `WSGIDaemonProcess` | nginx `proxy_read_timeout 60s`; Caddy `response_header_timeout 60s` | nginx `proxy_read_timeout 60s` | IIS `requestTimeout`; Caddy as Linux; Apache `ProxyPass ... timeout=60` |

## 1. Is the new code deployed?

    cd /var/lib/planetGen && git log -1 --oneline
    grep __version__ src/planetgen/_version.py

Pass: the version you meant to deploy (the **Version** badge at the top
of `README.md` on `main`).

## 2. Has the database been migrated?

Reloading or restarting the app does not apply schema migrations. Only
`planetgen.cli.migrate` (or `update.sh`/`install.sh`, and `update.ps1`/`install.ps1`
on Windows, which call it) does.

    curl -s https://HOST/api/health

Pass: `"schema_current": true`. If not: run `python3 -m planetgen.cli.migrate`
with the same database settings the site uses (`config.json`, or the
same `PLANETGEN_MYSQL_*` variables), then check again.
`python3 -m planetgen.cli.migrate --status` prints the current and target
versions without changing anything.

If you use more than one database, check each one. `planetgen.cli.migrate` only
migrates the one it is pointed at:

    curl -s "https://HOST/api/health?db=OTHER_DB_NAME"
    python3 -m planetgen.cli.migrate --mysql-database OTHER_DB_NAME

## 3. Are the spatial indexes there?

    mysql DB_NAME -e "SELECT table_name, index_name FROM information_schema.statistics
      WHERE table_schema = DATABASE() AND index_name LIKE 'idx_%_center' AND seq_in_index = 1"

Pass: nine rows: `sectors`, `black_holes`, `neutron_stars`, `nebulae`,
`supernova_remnants`, `quasars`, `rogue_planets`, `interstellar_comets`
and `asteroid_fields`.

## 4. Is the web server config current?

Check the **request time limit** row for your platform. On Apache:

    grep -n "WSGIDaemonProcess planetgen-api" /etc/apache2/sites-available/planetgen.conf

Pass: the 60-second limit is there (compare with your platform's example
under `examples/`). Without it, one runaway request can hold an app
thread for a long time.

Behind nginx, Caddy, IIS or Apache's `mod_proxy` (every setup except
Apache + mod_wsgi), also check that `config.json` has `proxy_fix` set
(see [`deployment/README.md`](deployment/README.md#behind-a-reverse-proxy-proxy_fix)):

    curl -sI https://HOST/ | grep -i strict-transport

Pass: a `Strict-Transport-Security` header. Without `proxy_fix`, the app
thinks every request is plain HTTP and every visitor shares one
rate-limit budget.

## 5. Is the tile cache writable?

The Galaxy Map caches tiles on disk, in `/var/cache/planetgen/tiles` by
default (see `tile_cache` in [`config.md`](config.md)).

    ls -ld /var/cache/planetgen/tiles /var/lib/planetGen/jobs

Pass: both folders exist and are owned by the **web user**. The tile
cache fills up after you open the Galaxy Map.

## 5a. Are the code and config.json locked down?

    ls -ld src/html src/html/wsgi.py config.json /var/log/planetgen.log

Pass: `src/html` and `wsgi.py` are owned by `root` with the web user's
group (`www-data`, or `_www` on macOS), not by the web user; `config.json`
is `-rw-r----- root www-data` (mode 640: it holds the database password
and `secret_key`); the debug log, if there is one, is `-rw-rw----` owned
by the web user and group (mode 660), never world-writable. If not, on
Linux: `sudo ./update.sh`, or `sudo examples/apache/set-permissions.sh`
and `sudo examples/apache/setup-debug-log.sh`. On macOS and Windows,
repeat the permissions step of your guide. Anyone who runs the generator
from a shell without `sudo` must be in the web user's group to read
`config.json` and append to the debug log.

## 5b. Has the first admin login been changed?

Log in at `https://HOST/login`. On a fresh install the username is
`admin` and the password is the random one `planetgen.cli.migrate` printed once
(there is no default password).

Pass: after logging in you land on `/admin`, not the forced "Change
Credentials" page. If you're sent to `/account`, choose a new username
and password now. If nobody has the printed password, see "Resetting the
admin login" in [`api.md`](api.md#resetting-the-admin-login).

## 6. Time the Galaxy Map's API calls

    time curl -s -o /dev/null "https://HOST/api/galaxy/stamp"
    time curl -s -o /dev/null "https://HOST/api/galaxy/tiles?tiles=0/0/0/0"

Pass: each finishes in a few seconds at most.

## 7. Reproduce real use while watching the logs

Follow the **app error log** for your platform, for example on Apache:

    sudo tail -f /var/log/apache2/planetgen_error.log | grep -E "Timeout|Truncated|MemoryError|Traceback"

In a second terminal, watch the app's memory with the **app processes**
command, for example:

    watch -n 2 "ps -o rss,cmd -C apache2 | sort -n | tail -3"

Open the Galaxy Map, click empty space well away from the core so the view
re-centers there, then zoom and rotate for a minute. While the map updates,
open another page (for example Browse) in a second tab.

Pass: nothing new in the error log, memory stays flat instead of climbing
by hundreds of MB, and the other page loads normally.
