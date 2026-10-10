# Installing planetGen

This guide takes you from nothing to a running planetGen site with a
generated galaxy, on Linux or macOS. It covers the install
scripts, where the example configs live, and what to do after the first
install. For what planetGen is and how to use it, see the
[README](README.md).

## The whole path at a glance

1. [Check the requirements](#requirements).
2. [Get the code](#1-get-the-code) as a git checkout.
3. [Create the MySQL database and account](#2-create-the-database-and-account).
4. [Write `config.json`](#3-write-configjson).
5. [Run the installer](#4-run-the-installer) for your platform.
6. [Set up the web server](#5-set-up-the-web-server) from the example config.
7. [Turn on HTTPS and change the admin login](#6-https-and-the-first-admin-login).
8. [Generate the galaxy](#7-generate-the-galaxy).
9. [Check the server](#8-check-the-server).
10. [Keep it up to date](#keeping-it-up-to-date).

## Requirements

| What | Version | Notes |
|---|---|---|
| Python | 3.9+ on Linux, 3.10+ on macOS (Homebrew's) | Linux uses the system Python. |
| MySQL or MariaDB | MySQL 8.0.16+ or MariaDB 10.4+ | Local or on another machine. |
| A web server | Apache2 + mod_wsgi (reference), nginx, Caddy | Only for the website. The command line needs none. |
| git | any | The checkout must be a git clone so the update scripts can pull. |
| Disk | about 10 GB once the galaxy is planned, more as sectors fill | The plan stores about 60 million bright stars. Most people then fill neighborhoods as they need them rather than the whole galaxy. |

Supported platforms:

| Platform | Installer | Updater | Guide |
|---|---|---|---|
| Debian or Ubuntu (reference) | `install.sh` | `update.sh` | [Apache](docs/deployment/apache.md), [nginx](docs/deployment/nginx.md), [Caddy](docs/deployment/caddy.md) |
| macOS | `install.sh` | `update.sh` | [macOS](docs/deployment/macos.md): Homebrew nginx in front of gunicorn under launchd |
| A rented server | as on Linux | as on Linux | [VPS and PaaS](docs/deployment/paas.md) |

[`docs/deployment/README.md`](docs/deployment/README.md) compares the
platforms in more detail. Windows: see [docs/WINDOWS.md](docs/WINDOWS.md).

> Do not install planetGen with `pip install .`. The installers put its
> libraries where the web server can use them: from apt first on Linux,
> or in a virtual environment built from a hash-checked lock file on
> macOS. They then install the checkout itself as an
> editable package (`pip install -e`), so the site and the command line
> run the checkout's code and a `git pull` takes effect without a
> reinstall.

## 1. Get the code

Clone the repository where the site will run from. The guides use
`/var/lib/planetGen`:

```bash
sudo git clone https://github.com/dwhagar/planetGen.git /var/lib/planetGen
cd /var/lib/planetGen
```

Keep it a git checkout: the update scripts update it by fetching and
resetting to the tracked branch.

On macOS, install the Homebrew packages first:
`brew install python@3.12 nginx mariadb` (or `mysql@8.4` instead of
`mariadb`).

## 2. Create the database and account

planetGen uses one MySQL account for everything: the generator, the
migrations and the website. It creates its own tables and a second
*control* schema (`planetgen_control`: admin logins, sessions, API keys,
audit log), so it needs full rights on every `planetgen%` schema:

```sql
CREATE DATABASE planetgen CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci;
CREATE USER 'planetgen'@'localhost' IDENTIFIED BY 'a long random password';
GRANT ALL PRIVILEGES ON `planetgen%`.* TO 'planetgen'@'localhost';
```

Use `'planetgen'@'%'` (or the web host's address) when MySQL is on
another machine. The admin Generate page's *New galaxy* and *Reset*
actions empty the tables, so a less privileged web account still needs
`DROP`. More in [MySQL accounts](docs/deployment/README.md#mysql-accounts).

## 3. Write `config.json`

```bash
sudo cp config.json.example config.json
sudo ${EDITOR:-nano} config.json
```

At minimum set:

- `mysql.host`, `mysql.user`, `mysql.password`, `mysql.database`;
- `secret_key` to a long random string:
  `python3 -c "import secrets; print(secrets.token_hex(32))"`;
- `site_name` and `base_url` for your site;
- `proxy_fix` when a separate web server proxies to the app (nginx,
  Caddy, macOS): `{"x_for": 1, "x_proto": 1, "x_host": 0}`
  ([why](docs/deployment/README.md#behind-a-reverse-proxy-proxy_fix)).
  Apache with mod_wsgi does not need it.

Every setting, and the `PLANETGEN_*` environment variable that overrides
it, is in [`docs/config.md`](docs/config.md). `config.json` is gitignored,
so updates never overwrite it, and the installers make it readable only
by administrators and the web server's account.

## 4. Run the installer

The installer prints the first admin password **once**, in a box, when
it creates a new database. Write it down.

### Linux (Debian or Ubuntu)

```bash
sudo ./install.sh
```

It runs eight steps and prints each one:

1. Installs the Python libraries. On a managed Python (Ubuntu 24.04+,
   Debian 12+) they come from apt, with pip only for a library apt lacks
   or ships too old, installed system-wide alongside apt's; no virtual
   environment ([Managed Python](docs/deployment/apache.md#managed-python)).
   It also adds a `planetgen` command in `/usr/local/bin`.
2. Creates or migrates the database schema, and prints the admin
   password on a new database.
3. Fetches the NLTK `words` corpus into `/usr/local/share/nltk_data`.
4. Makes the scripts executable.
5. Installs mod_wsgi if needed and enables Apache's `wsgi`, `headers` and
   `deflate` modules.
6. Locks down the code tree and `config.json` for Apache's account.
7. Creates the Galaxy Map tile cache (`/var/cache/planetgen/tiles`) and
   the Generate jobs directory (`/var/lib/planetGen/jobs`).
8. Creates the debug log and its logrotate config when `debug` is on.

After step 3 it also offers to run the optional population pass
(species, civilizations and territories; y/N, defaulting to N after 30
seconds). `sudo POPULATION=1 ./install.sh` runs it without asking, and
`planetgen population` runs it any time later. A new galaxy has no
worlds yet, so there is nothing for it to do on a first install.

`sudo ./install.sh --skip-database` leaves out step 2 and the population
question when the database
isn't ready yet; run `sudo ./update.sh` once it is. The installer never
writes the web server's site file; it prints what to do next.

### macOS

```bash
sudo /var/lib/planetGen/install.sh
```

The same script as on Linux, with macOS's own pieces: a virtual
environment at `/usr/local/planetgen/venv` from Homebrew's `python3`
with the libraries and gunicorn from `requirements-server.lock` (checked
by hash), `_www` ownership, `newsyslog` for the debug log, and gunicorn
installed and started as a launchd daemon on `127.0.0.1:8000`, plus the
`planetgen` command. nginx stays yours (next step). See the [macOS guide](docs/deployment/macos.md).

## 5. Set up the web server

Each platform has a guide and example files. Copy the example, set the
host name, check the paths, and enable it.

**Apache on Linux (the reference setup):**

```bash
sudo cp examples/apache/planetgen.conf.example /etc/apache2/sites-available/planetgen.conf
sudo ${EDITOR:-nano} /etc/apache2/sites-available/planetgen.conf   # set ServerName
sudo a2ensite planetgen
sudo systemctl reload apache2
```

**nginx on macOS:** copy `examples/macos/planetgen-nginx.conf` to
`$(brew --prefix)/etc/nginx/servers/planetgen.conf`, set `server_name`
and the certificate paths, put `user _www _www;` at the top of
`nginx.conf`, then `sudo nginx -t && sudo brew services start nginx`.

Where the examples live:

| Directory | What is in it | Guide |
|---|---|---|
| [`examples/apache/`](examples/apache/) | Apache vhost, and the permission, cache-directory and debug-log scripts the installers call | [apache.md](docs/deployment/apache.md) |
| [`examples/nginx/`](examples/nginx/) | nginx site file | [nginx.md](docs/deployment/nginx.md) |
| [`examples/caddy/`](examples/caddy/) | Caddyfile | [caddy.md](docs/deployment/caddy.md) |
| [`examples/systemd/`](examples/systemd/) | gunicorn service for nginx or Caddy | [nginx.md](docs/deployment/nginx.md) |
| [`examples/macos/`](examples/macos/) | nginx config and launchd plists (gunicorn, orbit update, `update.sh`) | [macos.md](docs/deployment/macos.md) |
| [`examples/maintenance/`](examples/maintenance/) | Monthly maintenance: systemd timers or launchd daemons (`install-maintenance-timer.sh`) | [below](#scheduled-maintenance) |
| [`examples/systems/`](examples/systems/) | Example system files for `planetgen system --system-file` | [system-file-format.md](docs/system-file-format.md) |

## 6. HTTPS and the first admin login

The admin pages only work over HTTPS. With Apache, `sudo certbot
--apache`; each guide says how for its server.

Then open `https://<your site>/login` and sign in as `admin` with the
password the installer printed. You are sent to `/account` to choose
your own username and password; the other admin pages stay locked until
you do. Lost the password? See
[Resetting the admin login](docs/api.md#resetting-the-admin-login).

## 7. Generate the galaxy

A new site has an empty galaxy. Use the browser or the command line.

**In the browser:** as an admin, open **Generate** (`/admin/generate`)
and choose *New galaxy*. It plans the galaxy's shape, places its bright
stars (shown as their own step with a progress bar), and generates a
first neighborhood around a random start, as a background job you can
watch and cancel. *Generate sectors* adds more
later, and so do the Generate buttons on the Galaxy Map and Sector Map.

**On the command line**, from the checkout:

```bash
planetgen plan      # the galaxy's shape and its bright stars; once, before anything else
planetgen galaxy    # a random start and every sector within 100 ly of it
```

The installer adds the `planetgen` command.

**Planning takes a while and needs space.** The plan places every star
of 500 solar luminosities or more across the whole galaxy, about 60
million stars and 10 GB in the database, before any sector is filled,
so the Galaxy Map shows the spiral arms from the start: roughly 20
minutes of drawing plus the database load. For a quick test galaxy, raise the
threshold (`--bright-star-min-luminosity 5000`) or skip the stars
(`--no-bright-stars`); `--bright-stars-only` places them later on the
stored plan. Going the other way, `--bright-star-min-luminosity 100`
places about 220 million (about 35 GB).

Then browse the **Galaxy Map**, the sectors and the systems on the site.
See [Usage](README.md#usage) for the other commands.

## 8. Check the server

Run the [server checklist](docs/server-checklist.md): it confirms the
code, schema, indexes, web server config, cache and permissions are all
right, with commands for each platform.

## Keeping it up to date

```bash
sudo ./update.sh
```

It pulls the latest code (a `git reset --hard` to the tracked
branch; `config.json` and other untracked files are kept), checks the
Python libraries and installs only a missing or too-old one, migrates
the database, re-applies folders and permissions, and checks that the
site imports. It never reinstalls what is already there.

When the update includes a schema migration, it first asks whether to
**delete the galaxy data instead** (y/N, defaulting to N after 30
seconds; admin logins are kept either way). Answer `y` when the release
notes say to regenerate, then generate the galaxy again
([step 7](#7-generate-the-galaxy)). The update never runs the
population pass; run `planetgen population` by hand when wanted. Then
reload the site as the script's last line says (Apache reload on Linux,
a `SIGHUP` to gunicorn on macOS).

Do not use a plain `git pull`: it skips the migration and can drop the
scripts' executable bits.

### Scheduled maintenance

The monthly orbit update (`planetgen.cli.orbits`, at 03:30 on the 1st)
and, optionally, an unattended monthly update (at 03:00):

```bash
sudo examples/maintenance/install-maintenance-timer.sh                      # Linux (systemd) and macOS (launchd)
sudo examples/maintenance/install-maintenance-timer.sh --skip-update-timer  # orbits only
```

An unattended update deploys whatever is on the tracked branch with no review, and
keeps the data when a migration is pending.

## Command line only

The generator does not need the website. Run the installer anyway (it
sets up the libraries, database and NLTK corpus) and skip step 5. Then
run `planetgen`, the command the installer adds.

## Development setup

For working on planetGen itself, use a virtual environment and an
editable install with the test extras:

```bash
python3 -m venv .venv
. .venv/bin/activate
pip install -e ".[test,api]"
pytest -n auto
python3 src/html/wsgi.py      # the site at http://127.0.0.1:5000/
```

Set `admin_cookie_insecure` in `config.json` to log in over plain HTTP
locally. Database tests need a MySQL server and skip without one; see
[`docs/testing.md`](docs/testing.md).
