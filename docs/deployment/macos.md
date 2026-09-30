# macOS

Apple discontinued macOS Server in 2022, and nothing Apple ships replaces
it for hosting a Python web app. A Mac hosts planetGen the way a Linux
server does, with Homebrew supplying MariaDB or MySQL, nginx and Python.
gunicorn runs the app under launchd, the macOS service manager. The job
runner behind the admin Generate page is POSIX code and works on macOS.

`install.sh`, `update.sh` and `install-maintenance-timer.sh` run on
macOS too (under its bash 3.2). See [The short way](#the-short-way);
the steps they automate are written out below for reference.

Paths match the Linux guides where they can:

| Path | What |
|---|---|
| `/var/lib/planetGen` | The checkout (on macOS `/var` is `/private/var`) |
| `/var/cache/planetgen/tiles`, `/var/lib/planetgen/jobs` | Tile cache and Generate jobs (the defaults) |
| `/var/log/planetgen.log` | Debug log (the default) |
| `/usr/local/planetgen/venv` | Python virtual environment |
| `/usr/local/planetgen/log` | gunicorn, nginx and orbit-update logs |
| `/usr/local/share/nltk_data` | NLTK corpus (same as Linux) |

The site runs as `_www`, macOS's built-in web server account.

## Files

| File | Goes to |
|---|---|
| [`examples/macos/org.planetgen.gunicorn.plist`](../../examples/macos/org.planetgen.gunicorn.plist) | `/Library/LaunchDaemons/` |
| [`examples/macos/planetgen-nginx.conf`](../../examples/macos/planetgen-nginx.conf) | `$(brew --prefix)/etc/nginx/servers/planetgen.conf` |
| [`examples/macos/org.planetgen.orbits.planetgen.plist`](../../examples/macos/org.planetgen.orbits.planetgen.plist) | `/Library/LaunchDaemons/` |

## The short way

Install the Homebrew packages (step 1), create the database, clone the
checkout to `/var/lib/planetGen`, write `config.json` (step 4), then:

    sudo /var/lib/planetGen/install.sh

It makes the venv at `/usr/local/planetgen/venv` from Homebrew's
`python3` (3.10 or later) with the libraries and gunicorn from
`requirements-server.lock` (checked by hash), fetches the NLTK corpus,
runs the migration, applies the permissions of step 6 for `_www`, sets
up `newsyslog` for the debug log, and installs and starts the gunicorn
daemon of step 7. nginx (step 8) stays yours. `sudo ./update.sh` later
pulls and checks everything the same way, and
`sudo examples/maintenance/install-maintenance-timer.sh [database ...]`
installs the monthly orbit update (and update.sh, unless
`--skip-update-timer`) as launchd daemons.

## Install

1. **Homebrew packages:**

       brew install python@3.12 nginx mariadb     # or mysql@8.4 instead of mariadb
       brew services start mariadb

   Homebrew's MariaDB lets your macOS user in as a database admin
   (`sudo mariadb` or `mariadb -u $(whoami)`). Create the database and
   account as in [`README.md`](README.md#mysql-accounts).

2. **Checkout and Python:**

       sudo git clone https://github.com/dwhagar/planetGen.git /var/lib/planetGen
       sudo mkdir -p /usr/local/planetgen/log
       sudo $(brew --prefix)/bin/python3.12 -m venv /usr/local/planetgen/venv
       sudo /usr/local/planetgen/venv/bin/pip install "/var/lib/planetGen[api]" gunicorn
       sudo /usr/local/planetgen/venv/bin/pip uninstall -y planetGen

   A virtual environment, because Homebrew's Python refuses system-wide
   `pip install` (PEP 668). The uninstall removes the installed copy of
   planetGen itself: the app runs from the checkout, as on Linux.

3. **NLTK corpus:**

       sudo /usr/local/planetgen/venv/bin/python -c "import nltk; nltk.download('words', download_dir='/usr/local/share/nltk_data')"
       sudo chmod -R a+rX /usr/local/share/nltk_data

4. **config.json:** `sudo cp /var/lib/planetGen/config.json.example
   /var/lib/planetGen/config.json`, then set `mysql.*`, `secret_key`
   (`python3 -c "import secrets; print(secrets.token_hex(32))"`) and the
   proxy setting, because nginx sits in front of gunicorn
   ([`README.md`](README.md#behind-a-reverse-proxy-proxy_fix)):

       "proxy_fix": {"x_for": 1, "x_proto": 1, "x_host": 0}

   The default paths work on macOS. `jobs.python` can stay empty:
   gunicorn runs from the venv, so the app finds its own Python.

5. **Schema and first admin login:**

       cd /var/lib/planetGen && sudo /usr/local/planetgen/venv/bin/python src/migrateDb.py

   Note the `admin` password it prints.

6. **Permissions.** The same rules `set-permissions.sh`,
   `create-cache-dir.sh` and `setup-debug-log.sh` apply on Linux: the code
   tree belongs to root and is only readable by the web group, and the
   runtime directories belong to the web user.

       sudo chown -R root:_www /var/lib/planetGen/src/html
       sudo find /var/lib/planetGen/src/html -type d -exec chmod 750 {} +
       sudo find /var/lib/planetGen/src/html -type f -exec chmod 640 {} +
       sudo find /var/lib/planetGen/src/html -type f -name '*.py' -exec chmod 750 {} +
       sudo chown root:_www /var/lib/planetGen/config.json && sudo chmod 640 /var/lib/planetGen/config.json
       sudo install -d -o _www -g _www -m 750 /var/cache/planetgen/tiles /var/lib/planetgen/jobs
       sudo touch /var/log/planetgen.log && sudo chown _www:_www /var/log/planetgen.log && sudo chmod 660 /var/log/planetgen.log
       sudo chown _www:_www /usr/local/planetgen/log

   Anyone who runs the generator from a shell without `sudo` must be in
   the `_www` group to read `config.json` and append to the debug log.

7. **gunicorn under launchd:**

       sudo cp /var/lib/planetGen/examples/macos/org.planetgen.gunicorn.plist /Library/LaunchDaemons/
       sudo chown root:wheel /Library/LaunchDaemons/org.planetgen.gunicorn.plist
       sudo chmod 644 /Library/LaunchDaemons/org.planetgen.gunicorn.plist
       sudo launchctl bootstrap system /Library/LaunchDaemons/org.planetgen.gunicorn.plist
       curl -s http://127.0.0.1:8000/api/health

   The plist runs `gunicorn --pythonpath /var/lib/planetGen/src/html
   wsgi:application` as `_www`, with one process and five threads (like
   Apache's daemon), on 127.0.0.1:8000, and restarts it if it exits.
   `--no-control-socket` avoids a gunicorn 25.1+ crash on macOS (its
   control socket thread trips macOS's fork-safety check); remove it only
   for gunicorn older than 25.1.

8. **nginx:**

       cp /var/lib/planetGen/examples/macos/planetgen-nginx.conf $(brew --prefix)/etc/nginx/servers/planetgen.conf

   Edit `server_name` and the certificate paths. To use ports 80 and 443,
   nginx must start as root, and its workers should run as `_www` so they
   can read `static/`. At the top of `$(brew --prefix)/etc/nginx/nginx.conf`:

       user _www _www;

   and remove or change the default `server` on port 8080 there if you
   don't want it. Then:

       sudo nginx -t
       sudo brew services start nginx

   For a certificate: `brew install certbot`, then
   `sudo certbot certonly --standalone -d <name>` with nginx stopped (or
   `--webroot`), and point the config at the files; or use your own. On a
   LAN-only Mac, `mkcert` makes a locally trusted one.

9. Log in at `https://<your site>/login` and change the admin username
   and password. Then run the [server checklist](../server-checklist.md).

The nginx rules are those of the [nginx guide](nginx.md): `/static/` from
disk with nosniff and the `?v=` caching, gzip, a 1 MB body limit, 60 s
proxy timeouts, and everything else to gunicorn. gunicorn listens on
127.0.0.1 only, so only local processes can set the forwarded headers
`proxy_fix` trusts.

## Keeping it running

- **Update:** `sudo /var/lib/planetGen/update.sh`, then the `SIGHUP`
  below. By hand:

      cd /var/lib/planetGen && sudo git pull
      sudo /usr/local/planetgen/venv/bin/pip install --upgrade "/var/lib/planetGen[api]" gunicorn
      sudo /usr/local/planetgen/venv/bin/pip uninstall -y planetGen
      sudo /usr/local/planetgen/venv/bin/python src/migrateDb.py
      # re-apply step 6 if new files came in, then reload:
      sudo launchctl kill SIGHUP system/org.planetgen.gunicorn

  HUP loads the new code and leaves a running Generate job alone. A full
  restart is `sudo launchctl kickstart -k system/org.planetgen.gunicorn`.
- **Monthly orbit update:**
  `sudo examples/maintenance/install-maintenance-timer.sh planetgen`
  installs
  [`org.planetgen.orbits.planetgen.plist`](../../examples/macos/org.planetgen.orbits.planetgen.plist)
  (one per database named) and
  [`org.planetgen.update.plist`](../../examples/macos/org.planetgen.update.plist),
  or install them by hand the same way as the gunicorn one. It runs
  at 03:30 on the 1st, like the Linux timer. A run missed while the Mac
  was switched off is skipped: launchd has no equivalent of systemd's
  `Persistent=true`.
- **Debug log:** `install.sh` writes this `newsyslog` rule to
  `/etc/newsyslog.d/planetgen.conf` (7 copies, rotate past 100 MB,
  compressed); by hand it is:

      /var/log/planetgen.log  _www:_www  660  7  102400  *  Z

- **Sleep:** a Mac serving a site should not sleep (System Settings,
  Energy, "Prevent automatic sleeping"), and should start after a power
  failure.
- **Logs:** `/usr/local/planetgen/log/gunicorn.log`, the nginx logs next
  to it, and the debug log.
- **Stats page:** the admin Stats page shows no memory figure on macOS
  (it reads Linux's `/proc/meminfo`).

## Checks

The [server checklist](../server-checklist.md) applies, with the launchd
commands it lists, plus:

    curl -sI "https://HOST/static/style.css?v=1" | grep -i -e cache-control -e nosniff
    curl -sI https://HOST/ | grep -i strict-transport
