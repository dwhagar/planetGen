# PlanetGen cloud environment setup

This folder (`examples/cloud-environment/`) holds `setup-planetgen-env.sh` and this README.

`setup-planetgen-env.sh` gives every new Claude cloud session a ready PlanetGen test machine: MySQL 8, Apache 2 with mod_wsgi, Redis, Flask and every other Python dependency, the NLTK words corpus, and a Playwright that matches the Chromium already in the image.

## Where to paste it

1. Open the project's **Project settings**, then the **Cloud environment** menu, then the gear beside the environment, and find **Setup script**.
2. Paste the whole contents of `setup-planetgen-env.sh`.
3. Network access must be **Trusted** (the default) or broader. It needs the Ubuntu archive, PyPI, GitHub and `raw.githubusercontent.com`, all on the Trusted list.
4. Add these under **Environment variables** (belt and braces; the script also writes them to `/etc/profile.d/planetgen-test.sh`, but a variable set here is certain to reach Claude's shell):

```
PLANETGEN_TEST_MYSQL_HOST=127.0.0.1
PLANETGEN_TEST_MYSQL_PORT=3306
PLANETGEN_TEST_MYSQL_USER=root
PLANETGEN_TEST_MYSQL_PASSWORD=planetgen-test
PLANETGEN_TEST_REDIS_URL=redis://127.0.0.1:6379/0
NLTK_ALLOW_PROXIED_URLOPEN=1
PIP_BREAK_SYSTEM_PACKAGES=1
```

`planetgen-test` is a throwaway password for a database that only exists inside the disposable VM; anyone who can use the environment can read these variables, so do not reuse a real password. Change it in the script and here together if you want a different one.

## How the cloud runs a setup script (from the environment docs)

- It runs once, as root on Ubuntu 24.04, before Claude starts, then the **filesystem is snapshotted** and later sessions start from that snapshot without running the script again. It reruns when you edit the script or the network settings, and after about 7 days.
- It must exit 0 and finish in about 5 minutes, or the snapshot is not kept (this one takes about 70 seconds on a fresh machine).
- **Running processes are not kept.** MySQL, Redis and Apache therefore cannot be "already running" in a new session. The script installs and configures them, shuts MySQL down cleanly so the snapshot is tidy, and leaves one command for each session: `planetgen-env-up` (about 3 seconds).
- A script that fails stops the session from starting. The script exits non-zero on any error and prints the failing line; the full log is `/var/log/planetgen-setup.log`.

## Starting the services in each session

Run `planetgen-env-up`. It starts MySQL, Redis and Apache, points the editable install of planetGen at the session's checkout (`pip install -e . --no-deps`), and prints one status line. It is safe to repeat. `planetgen-env-up --check` only reports.

The script also tries to make this automatic two ways, neither of which I could prove from outside a fresh session: it adds a note to `/root/.claude/CLAUDE.md` telling Claude to run it before tests, and it adds a SessionStart hook to `/root/.claude/settings.json`. If a new session shows services down, ask Claude to run `planetgen-env-up`; if that works every time, nothing more is needed.

## What each part does

| Part | What it installs or does |
| --- | --- |
| apt | `mysql-server` (Ubuntu's MySQL 8.0, not MariaDB), `apache2`, `libapache2-mod-wsgi-py3`, `redis-server`. A temporary `policy-rc.d` stops apt starting services during the build. |
| MySQL | `zz-planetgen-test.cnf` (no binary log, no performance schema, relaxed flushing, 500 connections) for a fast throwaway server. Root gets a password with `caching_sha2_password`, MySQL 8's default login, so tests exercise the same path as your real server. A `planetgen` user with an empty password (localhost only) is also created, which is what the test suite falls back to if no `PLANETGEN_TEST_MYSQL_*` variables are set. |
| Apache | `headers`, `deflate` and `wsgi` modules enabled, config checked with `apache2ctl configtest`. It is installed and runs, but it does not serve PlanetGen; your site config stays on your real server. |
| Redis | Preinstalled in the image; installed only if missing. Started per session with no persistence. |
| Python | `pip install --ignore-installed -r requirements-server.lock` (exact pinned versions), then `.[test,api,browser]` from `setup.py` with `playwright==1.56.0` pinned to match the image's Chromium 1194. PlanetGen itself is removed from the snapshot and installed editable per session so it always points at that session's checkout. |
| NLTK | The `words` corpus in `/usr/local/share/nltk_data`, fetched with `NLTK_ALLOW_PROXIED_URLOPEN=1` (the agent proxy blocks it otherwise). |
| Playwright | Never downloads a browser. If the image's Chromium revision ever changes, set `PLAYWRIGHT_VERSION` on the environment to the matching Playwright release. |

## Where the requirements come from

The script reads `requirements-server.lock` and `setup.py` from a checkout if it finds one (`$PLANETGEN_REPO`, the current directory, `/home/user/planetGen` and a few other places). If none exists yet, because setup runs before the clone, it makes a shallow clone of `https://github.com/dwhagar/planetGen` into `/opt/planetgen-seed`. The versions are therefore whatever `main` had when the snapshot was last built; when the lock changes, the script's "rebuild after about 7 days" catches up, or edit any character in the script to force a rebuild.

## Not covered

- Only what the repo's tests need. It does not run `install.sh` or `update.sh`, does not create a PlanetGen database or `config.json`, and never touches your real server.
- Wiki.js and MediaWiki integration tests stay skipped (they need external servers).
- Python is the image's system Python (3.13 at the time of writing; PlanetGen needs 3.12 or newer, the script checks). `mod_wsgi` from Ubuntu is built for Ubuntu's 3.12, so it is not used to serve the app here.
