# Apache2 + mod_wsgi on Linux

The reference setup, and the one `install.sh` and `update.sh` are built
for: Debian or Ubuntu, Apache2 with mod_wsgi running the Flask app in its
own daemon process, and Apache serving `/static/` itself. For other web
servers and operating systems, see the [deployment index](README.md).
Apache + mod_wsgi needs no `proxy_fix` setting: Apache owns the client's
connection, so the app already sees the real address and scheme.

## Steps

1. Create the MySQL database and account ([`README.md`](README.md#mysql-accounts)).
2. Clone the repo to `/var/lib/planetGen`, copy `config.json.example` to
   `config.json` and fill in `mysql.*` and `secret_key`
   ([`config.md`](../reference/config.md)).
3. `sudo ./install.sh` from the checkout (it also installs and starts
   `redis-server`; the Generate page cannot start a job without Redis, so
   check it with `redis-cli ping`, which answers `PONG`). Note the admin password it
   prints in step 2/8. After the NLTK step it asks whether to run the
   optional population pass (y/N, default N; see
   [below](#migrating-or-deleting-the-database)).
4. Copy [`planetgen.conf.example`](../../examples/apache/planetgen.conf.example)
   to `/etc/apache2/sites-available/planetgen.conf`, set `ServerName`,
   then `sudo a2ensite planetgen && sudo systemctl reload apache2`.
5. `sudo certbot --apache` for HTTPS (the admin pages need it), then log
   in at `/login` and change the admin username and password.
6. Optional: the monthly maintenance timers in
   [`examples/maintenance/`](../../examples/maintenance/).
7. Run the [server checklist](server-checklist.md).

Later updates: `sudo ./update.sh`. When Apache is running it tests the configuration (`apache2ctl configtest`), reloads
it (a graceful reload loads new modules too; it restarts only when it just installed mod_wsgi), and checks that
`/api/health` answers (`PLANETGEN_HEALTH_URL` names another address). It asks `http://127.0.0.1/api/health` first; if that does not answer (Apache may give the bare address to another virtual host, OPS.77) it asks the site's own address, `base_url` in `config.json` plus `api/health`, pointed at this machine, and passes on a 2xx answer, so keep `base_url` set to the name in your Apache site. A `planetgen-gunicorn` unit, or the launchd
gunicorn daemon on macOS, is sent a reload the same way. If a step fails the update ends non-zero and prints the
command to run by hand; if only the `/api/health` check gets no good answer (the reload itself worked, and the database is migrated by then) it warns, says why ("HTTP 000" means no answer at all, and the warning names the cause) and opens the site; when nothing that serves the site is running it prints the commands and ends zero. A reload
aborts the requests a mod_wsgi daemon is serving after about 4 seconds, where `touch src/html/wsgi.py` does not.

## Deployment files

The files below live in
[`examples/apache/`](../../examples/apache/); nothing in that directory
is copied to the web server's document root.

Most deployments just need `sudo ./install.sh` (first-time setup) or
`sudo ./update.sh` (pulling a later update) from the repo root. Both
call `set-permissions.sh` and `create-cache-dir.sh` (below) as steps. Use
`set-permissions.sh` directly only to re-apply permissions on its own
(e.g. after adding files or a new `config.json` by hand).

| File | Purpose |
|---|---|
| [`planetgen.conf.example`](../../examples/apache/planetgen.conf.example) | Example Apache2 virtual host: `DocumentRoot` at `/var/lib/planetGen/src/html`, a `WSGIScriptAlias` mounting the Flask app (the API under `/api` plus every page, via `src/html/wsgi.py`, see [`api.md`](../reference/api.md#deploying) and [`html-interface.md`](../reference/html-interface.md#flask-pages)) at `/`, with `Alias /static/` for the static files, in its own `WSGIDaemonProcess` (not Apache's embedded/shared-process mode). Directory rules keep Apache out of the checkout: the whole install tree is denied, with no listings, no handlers and `AllowOverride None`; only `wsgi.py` and `static/` are let back in, and source, configuration, state, logs, `.git`, dotfiles and the virtual environment are refused by name wherever they sit (`test_apache_example_restrictions.py` runs the file against a real Apache). `update.sh` never edits your copy of this file: when an update changes the example it prints a note with the `diff` command to compare. Copy to `/etc/apache2/sites-available/`, edit, and enable with `a2ensite` -- the one step `install.sh` deliberately leaves manual. |
| [`set-permissions.sh`](../../examples/apache/set-permissions.sh) | Detects the user/group Apache2 is actually configured to run as (from `/etc/apache2/envvars`, a running `apache2` process, or falling back to the Debian/Ubuntu default `www-data:www-data` if neither is found -- e.g. because apache2 isn't started yet), then makes the whole checkout (`/var/lib/planetGen`, or the directory given as its first argument) and its `state/` directory (made if missing; the second argument moves it) belong to that user and group: directories 750 (the state directory 770, so an operator in the group can run the maintenance command), files 640, and every `*.py` and `*.sh` file at any depth (plus any file already owner-executable) 750. `config.json` is covered, so it is `www-data:www-data` mode 640 (it holds the database password and `secret_key`; see [`config.md`](../reference/config.md)); CLI users who run the generator without `sudo` must be in Apache's group to read it. `install.sh` and `update.sh` both call it on the whole checkout, so a fresh install and an updated one end up identical. Symbolic links are never followed, and it refuses a directory without `src/html`. The virtual environment and the NLTK corpus stay root's, because the root-run scripts execute that Python; the web server is kept out of the checkout by the `<Directory>` rules in `planetgen.conf.example`. Bash, not Python -- a Linux and macOS deployment step (`_www` on macOS), safe to re-run any time as root. |
| [`create-cache-dir.sh`](../../examples/apache/create-cache-dir.sh) | Creates the 3D Galaxy Map's on-disk tile cache (see [`config.md`](../reference/config.md)'s `tile_cache`: `PLANETGEN_TILE_CACHE_DIR`, else `tile_cache.dir`, else `/var/cache/planetgen/tiles`) and the admin Generate page's jobs directory (`PLANETGEN_JOBS_DIR`, else `jobs.dir`, else `/var/lib/planetGen/jobs`), and `chown`s them to Apache's user, detected the same way `set-permissions.sh` does (both source [`apache-identity.sh`](../../examples/apache/apache-identity.sh)). It runs as root, so it imports nothing from the repo: it reads the paths with [`deploy-paths.py`](../../examples/apache/deploy-paths.py), run as `python3 -I` (standard library only). Creates no tile cache when `cache.max_mb` is `0`. Takes the directory as an argument when it's set only by a `SetEnv` in the vhost, which the script can't see. Safe to re-run any time as root. |
| [`setup-debug-log.sh`](../../examples/apache/setup-debug-log.sh) | Creates the debug log (`log_file`, default `/var/log/planetgen.log`) whether or not `debug` is on, owned by Apache's user and group, mode 0660, and re-applies that to an existing one; checks Apache's user and group can write it and the activity log (with [`log-locations.py`](../../examples/apache/log-locations.py)), and for anything it can't set up only warns, with the commands that fix it. Installs a logrotate config for it (`/etc/logrotate.d/planetgen`: daily or past 100 MB, 7 compressed copies) and an hourly size check. Always creates the activity log (`log_dir`, default `/var/log/planetgen/planetgen.log`; folder root and Apache's group, mode 2770, file mode 0660) and its rotation (`/etc/logrotate.d/planetgen-log`: daily or past 100 MB, 30 copies). See [the activity log](../reference/config.md#the-activity-log). Safe to re-run any time as root. |

## planetGen's own Python

planetGen runs from its own virtual environment, `/usr/local/planetgen/venv`
(`PLANETGEN_VENV` changes the path), on Linux and macOS alike. `install.sh`'s
first step ([`scripts/install-python-deps.sh`](../../scripts/install-python-deps.sh))
makes it from the system's `python3` (3.12 or newer; Debian and Ubuntu's
`python3-venv` package is installed when `venv` needs it) and puts every
library `setup.py` asks for, and gunicorn, into it. The system Python and its
packages (apt's `python3-flask` and the rest) are not used and never touched,
so Debian 12+ and Ubuntu 23.04+ marking their Python as externally managed
(PEP 668) does not matter. Anything an older install left in the system
Python's site-packages (`/usr/local/lib/python3.X/dist-packages`) or in a
`/opt/planetgen/venv` is not used any more and can be removed
(`sudo apt purge` of the `python3-*` libraries planetGen once needed, and
`sudo rm -rf /opt/planetgen`).

Everything planetGen runs uses that venv:

- **Apache:** `WSGIDaemonProcess planetgen-api ... python-home=/usr/local/planetgen/venv`
  in the vhost ([`planetgen.conf.example`](../../examples/apache/planetgen.conf.example)).
  A vhost made from an older example needs that one option added;
  `update.sh` prints a warning with the exact line when it is missing.
- **The command line and background jobs:** the `planetgen` command
  (`/usr/local/bin/planetgen` links to the venv's) and
  `/usr/local/planetgen/venv/bin/python -m planetgen.cli.NAME`, which is what
  the orbit-update service runs. A bare `python3 -m planetgen.cli.NAME` finds
  none of the libraries.
- **gunicorn** (nginx and Caddy setups, macOS): `/usr/local/planetgen/venv/bin/gunicorn`.
- **`install.sh` and `update.sh`** themselves, once their Python step has run.

Everything pip installs comes from
[`requirements-server.lock`](../../requirements-server.lock): exact versions of
every library and its dependencies (gunicorn included), each with the sha256
hashes of its files, installed with `--require-hashes`. A new release on PyPI,
or a file replaced under an existing version, never reaches the server until
the lock is regenerated with `scripts/lock-requirements.sh` (after changing
`setup.py`, or with `--upgrade` to take new releases) and committed.

planetGen itself is installed as an editable package of the checkout
(`pip install -e`, into the venv), so the web app, the `planetgen` command and
every `-m planetgen.cli.NAME` tool run the checkout's code, and a `git pull`
takes effect without a reinstall.

`update.sh` reinstalls nothing. It runs the same script with `--check`, which
imports each library with the venv's python and compares its version with
`setup.py`'s floor. It makes the venv again, and installs from the lock, only
when the venv is missing or a library is missing, too old or fails to import.
It prints one line per library and stops the update if anything is still
unusable. Its last step imports the web app as Apache's user, so a library
www-data can't read shows up there rather than as a 500. `sudo ./install.sh` is
still the full reinstall.

### Migrating or deleting the database

When the configured database is behind the current schema, `install.sh`
and `update.sh` first ask:

    Delete all galaxy data in 'planetgen' instead of migrating it? [y/N] (default N in 30s):

`y` wipes every generated sector and system (the same as the Generate
page's Reset; admin logins are kept) and then brings the empty database to
the current schema. Anything else, no answer within 30 seconds, or no
terminal to ask on (the maintenance timer) keeps the data and migrates it.
Nothing is asked when the database is already current. The migration
shows a progress bar with one tick per step, the elapsed time and an
estimate of the time left. `python3 -m planetgen.cli.migrate --status` prints the
current and target versions and how many steps are pending, without
changing anything.

After the database step, both scripts offer the population pass:

    Run the population pass now (species, civilizations, territories)? [y/N] (default N in 30s):

`y` runs `planetgen population`, which names species, dates
civilizations and draws territories from what is already stored.
Anything else, no answer within 30 seconds, or no terminal skips it.
`sudo POPULATION=1 ./update.sh` (or `./install.sh`) runs it without
asking, and `planetgen population` runs it any time.
`install.sh --skip-database` skips both questions.

### What Apache needs to run the venv (and what update.sh checks)

mod_wsgi can run one app in its own virtual environment only in daemon
mode. `update.sh` (and `install.sh`) read your vhost, the enabled site that
holds `WSGIScriptAlias ... src/html/wsgi.py`, and **only warn**: each
warning prints the exact line to add, and neither script edits the vhost.

- **A daemon process group:** `WSGIDaemonProcess planetgen-api ... python-home=...`
  and the script alias with `process-group=planetgen-api` (or a
  `WSGIProcessGroup planetgen-api` line). A script alias with neither runs in
  Apache's embedded mode, which cannot use a venv.
- **`python-home` is the venv's root,** the folder that holds `bin/` and
  `lib/` (`/usr/local/planetgen/venv`), not `bin/python` and not a parent
  folder. A trailing `/` is fine.
- **The main interpreter:** `application-group=%{GLOBAL}` on the script
  alias (or `WSGIApplicationGroup %{GLOBAL}`), because numpy, scipy and numba
  load C code that does not work in a sub-interpreter.
- **Apache's user can read the venv:** it (`www-data` on Debian and Ubuntu,
  `apache` on RHEL and CentOS) must be able to enter every directory down
  to the venv and read the venv itself. The check runs as that user and names
  the first folder that blocks it, with the `chmod o+rx` that opens it. Home
  folders (mode 0700 by default on newer Ubuntu) are the usual cause.
- **The same Python version as mod_wsgi:** see below.

### Which Python and libraries Apache uses

The vhost names the venv with `python-home=/usr/local/planetgen/venv` on
`WSGIDaemonProcess`, so the daemon imports every library from there and none
from the system Python.

The one thing that has to line up is the Python version: mod_wsgi
(`libapache2-mod-wsgi-py3`) embeds the Python it was built for, and the venv
is made from the system's `python3`, so the two must be the same version. Both
scripts check this and warn if they differ. By hand:

    ldd /usr/lib/apache2/modules/mod_wsgi.so | grep libpython   # e.g. libpython3.12
    /usr/local/planetgen/venv/bin/python --version              # must match
    sudo -u www-data /usr/local/planetgen/venv/bin/python -c "import flask; print(flask.__file__)"

If they differ, install the `libapache2-mod-wsgi-py3` that matches, or make the
venv from the Python mod_wsgi was built for (remove the venv, then
`sudo PYTHON=/usr/bin/python3.12 ./update.sh`). A system Python upgrade that
moves to a new minor version needs the same.

To see what the running site really uses, open the admin Stats page
(`/admin/stats`): **Python** shows the daemon's version and prefix, and
**Libraries from** shows the directory it imports Flask from (it should be
under `/usr/local/planetgen/venv`). After an update,
`sudo systemctl reload apache2` restarts the daemon so it picks up anything
newly installed.
