# tests/test_deploy_scripts.py

"""
The root-run deployment scripts in `examples/apache/` (called by
`install.sh`/`update.sh`):

- `set-permissions.sh` leaves the code tree root-owned (read-only for
  Apache), gives Apache only the runtime `db/` directory, and locks
  `config.json` to `root:<apache group>` 640.
- `create-cache-dir.sh` imports nothing from the repo: it reads its paths
  with `deploy-paths.py` under `python3 -I`.
- `setup-debug-log.sh` makes the debug log 0660, never world-writable.
- `log-locations.py` (OPS.5) sets both logs up wherever they're
  configured, checks the web server's user and group can write them, and
  for anything it can't fix only warns, with the commands that would.

The scripts that only touch a directory they're given run for real
against a throwaway tree when the tests run as root and a `www-data`
user exists (the scripts' fallback Apache identity); otherwise those
tests are skipped.
"""

import grp
import json
import os
import pathlib
import pwd
import shutil
import stat
import subprocess
import sys

import pytest

REPO_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
APACHE_DIR = os.path.join(REPO_DIR, "examples", "apache")
DEPLOY_PATHS = os.path.join(APACHE_DIR, "deploy-paths.py")
SCRIPTS = ["set-permissions.sh", "create-cache-dir.sh", "setup-debug-log.sh", "apache-identity.sh"]


def _read(name):
    with open(os.path.join(APACHE_DIR, name), encoding="utf-8") as f:
        return f.read()


def _code_lines(name):
    """The script without its comment lines."""
    return "\n".join(line for line in _read(name).splitlines() if not line.lstrip().startswith("#"))


def _www_data():
    try:
        return pwd.getpwnam("www-data").pw_uid, grp.getgrnam("www-data").gr_gid
    except KeyError:
        return None


needs_root = pytest.mark.skipif(
    os.geteuid() != 0 or _www_data() is None or os.path.exists("/etc/apache2/envvars"),
    reason="needs root, a www-data user, and no real Apache config to detect",
)


@pytest.mark.parametrize("name", SCRIPTS + ["../../install.sh", "../../update.sh"])
def test_scripts_parse(name):
    subprocess.run(["bash", "-n", os.path.join(APACHE_DIR, name)], check=True)


# --- deploy-paths.py ---------------------------------------------------------

def _deploy_paths(repo_dir, env=None, script=DEPLOY_PATHS, cwd=None):
    clean = {k: v for k, v in os.environ.items() if not k.startswith("PLANETGEN_")}
    clean.update(env or {})
    proc = subprocess.run([sys.executable, "-I", script, str(repo_dir)], env=clean, cwd=cwd or "/",
                          capture_output=True, text=True, check=True)
    return proc.stdout.splitlines()


def test_deploy_paths_defaults_without_a_config(tmp_path):
    assert _deploy_paths(tmp_path) == ["/var/cache/planetgen/tiles", "/var/lib/planetgen/jobs"]


def test_deploy_paths_reads_config_json(tmp_path):
    (tmp_path / "config.json").write_text(json.dumps({
        "tile_cache": {"dir": "/srv/tiles"}, "jobs": {"dir": "/srv/jobs"},
    }))
    assert _deploy_paths(tmp_path) == ["/srv/tiles", "/srv/jobs"]


def test_deploy_paths_env_wins_and_max_mb_0_turns_the_cache_off(tmp_path):
    (tmp_path / "config.json").write_text(json.dumps({"tile_cache": {"dir": "/srv/tiles", "max_mb": 0}}))
    assert _deploy_paths(tmp_path, {"PLANETGEN_JOBS_DIR": "/env/jobs"}) == ["", "/env/jobs"]
    assert _deploy_paths(tmp_path, {"PLANETGEN_TILE_CACHE_MAX_MB": "5", "PLANETGEN_TILE_CACHE_DIR": "/env/t"}) == [
        "/env/t", "/var/lib/planetgen/jobs"]


def test_deploy_paths_matches_the_web_apps_own_answer(tmp_path, monkeypatch):
    """Mirrors `tilecache.configured_cache_dir` / `jobs.configured_jobs_dir`."""
    sys.path.insert(0, os.path.join(REPO_DIR, "src", "html", "lib"))
    import tilecache
    from web import jobs

    config = {"tile_cache": {"dir": "/srv/tiles", "max_mb": 50}, "jobs": {"dir": "/srv/jobs"}}
    (tmp_path / "config.json").write_text(json.dumps(config))
    monkeypatch.setattr(tilecache, "_config", lambda: config["tile_cache"])
    monkeypatch.setattr(jobs, "_config", lambda: config["jobs"])
    for var in ("PLANETGEN_TILE_CACHE_DIR", "PLANETGEN_TILE_CACHE_MAX_MB", "PLANETGEN_JOBS_DIR"):
        monkeypatch.delenv(var, raising=False)
    assert _deploy_paths(tmp_path) == [tilecache.configured_cache_dir(), jobs.configured_jobs_dir()]


def test_deploy_paths_imports_nothing_planted_beside_it(tmp_path):
    """Under -I a `json.py` next to the script (or in the current
    directory) is never imported -- the root-run helper can't be made to
    run code from a directory the web user could write."""
    planted = tmp_path / "planted"
    planted.mkdir()
    script = planted / "deploy-paths.py"
    shutil.copy(DEPLOY_PATHS, script)
    marker = tmp_path / "imported"
    (planted / "json.py").write_text(f"open({str(marker)!r}, 'w').close()\nraise SystemExit(9)\n")
    assert _deploy_paths(tmp_path, script=str(script), cwd=str(planted)) == [
        "/var/cache/planetgen/tiles", "/var/lib/planetgen/jobs"]
    assert not marker.exists()


def test_root_helpers_run_python_isolated():
    code = _code_lines("create-cache-dir.sh")
    assert '"$PYTHON" -I "$APACHE_DIR/deploy-paths.py"' in code
    assert "sys.path" not in code and "import tilecache" not in code and "appconfig" not in code
    assert '"$PYTHON" -I -' in _code_lines("setup-debug-log.sh")
    assert '"$PYTHON" -I "$APACHE_DIR/log-locations.py"' in _code_lines("setup-debug-log.sh")


# --- setup-debug-log.sh ------------------------------------------------------

def test_debug_log_is_0660_never_world_writable():
    code = _code_lines("setup-debug-log.sh")
    assert "0666" not in code
    assert "create 0660 $APACHE_USER $APACHE_GROUP" in code
    helper = _read("log-locations.py")
    assert "0o666" not in helper and "0o777" not in helper
    assert "os.chmod, path, 0o660" in helper


def test_the_log_setup_never_stops_an_install_or_update():
    for script in ("install.sh", "update.sh"):
        with open(os.path.join(REPO_DIR, script), encoding="utf-8") as f:
            text = f.read()
        assert 'setup-debug-log.sh" \\\n    || echo "warning:' in text, script
    # Only running it without root, or with no Python at all, is an error.
    assert _code_lines("setup-debug-log.sh").count("exit 1") == 2


# --- log-locations.py (OPS.5), run as the current user -----------------------

LOG_LOCATIONS = os.path.join(APACHE_DIR, "log-locations.py")


def _me():
    return pwd.getpwuid(os.getuid()).pw_name, grp.getgrgid(os.getgid()).gr_name


def _run_log_locations(tmp_path, user=None, group=None, config=None, **env_overrides):
    """Runs log-locations.py against a fake checkout holding the real
    appconfig.py and, when given, a config.json."""
    repo = tmp_path / "repo"
    package = repo / "src" / "planetgen" / "util"
    package.mkdir(parents=True, exist_ok=True)
    shutil.copy(os.path.join(REPO_DIR, "src", "planetgen", "util", "appconfig.py"), package / "appconfig.py")
    if config is not None:
        (repo / "config.json").write_text(json.dumps(config))
    me_user, me_group = _me()
    env = {k: v for k, v in os.environ.items() if not k.startswith("PLANETGEN_")}
    env.update(env_overrides)
    return subprocess.run([sys.executable, "-I", LOG_LOCATIONS, str(repo), user or me_user, group or me_group],
                          env=env, capture_output=True, text=True, check=False)


@pytest.fixture
def open_dir():
    """A folder every user can pass through (pytest's own tmp_path sits
    in a 0700 folder, which the group check rightly reports)."""
    import tempfile
    path = tempfile.mkdtemp(prefix="planetgen-logs-", dir="/tmp")
    os.chmod(path, 0o755)
    yield pathlib.Path(path)
    shutil.rmtree(path, ignore_errors=True)


def test_log_locations_sets_up_both_logs_from_config(tmp_path, open_dir):
    logs = open_dir / "logs"
    result = _run_log_locations(tmp_path, config={"log_file": str(logs / "debug" / "planetgen.log"),
                                                  "log_dir": str(logs / "activity")})
    assert result.returncode == 0
    assert result.stderr == ""
    assert f'Debug log: {logs / "debug" / "planetgen.log"} (from "log_file" in config.json; debug is off' \
        in result.stdout
    assert f'Activity log: {logs / "activity" / "planetgen.log"} (from "log_dir" in config.json' in result.stdout
    for path in (logs / "debug" / "planetgen.log", logs / "activity" / "planetgen.log"):
        info = os.stat(path)
        assert stat.S_IMODE(info.st_mode) == 0o660, path
        assert (info.st_uid, info.st_gid) == (os.getuid(), os.getgid())
    assert stat.S_IMODE(os.stat(logs / "activity").st_mode) == 0o2770
    # Run again: nothing to change, still quiet.
    again = _run_log_locations(tmp_path, config={"log_file": str(logs / "debug" / "planetgen.log"),
                                                 "log_dir": str(logs / "activity")})
    assert again.returncode == 0 and again.stderr == ""


def test_log_locations_says_where_each_setting_came_from(tmp_path, open_dir):
    result = _run_log_locations(tmp_path, PLANETGEN_LOG_FILE=str(open_dir / "env.log"),
                                PLANETGEN_LOG_DIR=str(open_dir / "envdir"))
    assert "(from PLANETGEN_LOG_FILE;" in result.stdout
    assert "(from PLANETGEN_LOG_DIR;" in result.stdout


def test_an_unusable_log_path_only_warns_with_the_commands(tmp_path):
    blocker = tmp_path / "not-a-folder"
    blocker.write_text("")
    result = _run_log_locations(tmp_path, PLANETGEN_LOG_FILE=str(blocker / "debug.log"),
                                PLANETGEN_LOG_DIR=str(blocker / "activity"))
    assert result.returncode == 0
    user, group = _me()
    err = result.stderr
    assert f"warning: the debug log {blocker / 'debug.log'}" in err
    assert f"{blocker} is a file, not a folder" in err
    assert f"sudo mkdir -p {blocker}" in err
    assert f"sudo chown {user}:{group} {blocker / 'debug.log'}" in err
    assert f"sudo chmod 0660 {blocker / 'debug.log'}" in err
    assert '"log_file" in config.json (or PLANETGEN_LOG_FILE)' in err
    assert f"warning: the activity log {blocker / 'activity' / 'planetgen.log'}" in err
    assert f"sudo chown root:{group} {blocker / 'activity'}" in err
    assert f"sudo chmod 2770 {blocker / 'activity'}" in err
    assert '"log_dir" in config.json (or PLANETGEN_LOG_DIR)' in err
    assert "Traceback" not in err


def test_a_web_server_user_that_does_not_exist_only_warns(tmp_path):
    result = _run_log_locations(tmp_path, user="planetgen-no-such-user", group="planetgen-no-such-group",
                                PLANETGEN_LOG_FILE=str(tmp_path / "d" / "debug.log"),
                                PLANETGEN_LOG_DIR=str(tmp_path / "a"))
    assert result.returncode == 0
    assert "the web server's user 'planetgen-no-such-user' doesn't exist" in result.stderr
    assert "the web server's group 'planetgen-no-such-group' doesn't exist" in result.stderr
    assert f"sudo chown planetgen-no-such-user:planetgen-no-such-group {tmp_path / 'd' / 'debug.log'}" \
        in result.stderr


def test_a_broken_config_only_warns(tmp_path):
    result = _run_log_locations(tmp_path, config={"log_dir": 42})
    assert result.returncode == 0
    assert "couldn't work out where the logs go" in result.stderr
    assert '"log_dir" must be a string' in result.stderr


def _load_log_locations():
    import importlib.util
    spec = importlib.util.spec_from_file_location("log_locations", LOG_LOCATIONS)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_can_write_follows_every_folder_on_the_way(open_dir):
    module = _load_log_locations()
    folder = open_dir / "private"
    folder.mkdir()
    log = folder / "planetgen.log"
    log.write_text("")
    os.chmod(log, 0o660)
    owner, group = os.stat(log).st_uid, os.stat(log).st_gid
    other_uid, other_gid = owner + 4242, group + 4242
    os.chmod(folder, 0o700)
    # Neither a stranger nor the file's group gets through a 0700 folder...
    assert module.can_write(str(log), other_uid, {group}) == (False, str(folder))
    assert module.can_write(str(log), None, {group})[0] is False
    # ...the owner does, and root always does.
    assert module.can_write(str(log), owner, {group}) == (True, None)
    assert module.can_write(str(log), 0, set()) == (True, None)
    os.chmod(folder, 0o2770)
    assert module.can_write(str(log), other_uid, {group}) == (True, None)
    assert module.can_write(str(log), other_uid, {other_gid}) == (False, str(folder))
    os.chmod(folder, 0o755)
    os.chmod(log, 0o640)
    assert module.can_write(str(log), other_uid, {group}) == (False, str(log))
    assert module.can_write(str(folder / "missing.log"), owner, {group}) == (False, str(folder / "missing.log"))


# --- set-permissions.sh / create-cache-dir.sh, run for real ------------------

@pytest.fixture
def deployed(tmp_path):
    """A fake checkout: src/html with code (owned by www-data, as the old
    script left it), a legacy db/ dir and a world-readable config.json."""
    uid, gid = _www_data()
    repo = tmp_path / "planetGen"
    html = repo / "src" / "html"
    (html / "lib").mkdir(parents=True)
    (html / "wsgi.py").write_text("application = None\n")
    (html / "lib" / "fmt.py").write_text("")
    (html / "static").mkdir()
    (html / "static" / "style.css").write_text("")
    db = repo / "db"
    db.mkdir()
    (db / "old.db").write_text("")
    config = repo / "config.json"
    config.write_text("{}")
    os.chmod(config, 0o644)
    for root, dirs, files in os.walk(html):
        for name in dirs + files:
            os.chown(os.path.join(root, name), uid, gid)
    os.chown(html, uid, gid)
    return {"repo": repo, "html": html, "db": db, "config": config, "uid": uid, "gid": gid}


@needs_root
def test_set_permissions_makes_the_code_root_owned_and_config_640(deployed):
    subprocess.run([os.path.join(APACHE_DIR, "set-permissions.sh"), str(deployed["html"]), str(deployed["db"])],
                   check=True, capture_output=True)
    uid, gid = deployed["uid"], deployed["gid"]
    for root, dirs, files in os.walk(deployed["html"]):
        for name in [""] + dirs + files:
            info = os.lstat(os.path.join(root, name))
            assert (info.st_uid, info.st_gid) == (0, gid), os.path.join(root, name)
            assert not info.st_mode & (stat.S_IWGRP | stat.S_IWOTH | stat.S_IRWXO)
    assert stat.S_IMODE(os.stat(deployed["html"] / "wsgi.py").st_mode) == 0o750
    assert stat.S_IMODE(os.stat(deployed["html"] / "static" / "style.css").st_mode) == 0o640
    # Runtime data stays Apache's.
    assert os.stat(deployed["db"] / "old.db").st_uid == uid
    info = os.stat(deployed["config"])
    assert (info.st_uid, info.st_gid, stat.S_IMODE(info.st_mode)) == (0, gid, 0o640)


@needs_root
def test_set_permissions_without_a_config_json(deployed):
    deployed["config"].unlink()
    subprocess.run([os.path.join(APACHE_DIR, "set-permissions.sh"), str(deployed["html"]), str(deployed["db"])],
                   check=True, capture_output=True)
    assert not deployed["config"].exists()


@needs_root
def test_create_cache_dir_gives_apache_only_the_runtime_dirs(tmp_path):
    uid, gid = _www_data()
    env = {k: v for k, v in os.environ.items() if not k.startswith("PLANETGEN_")}
    env.update(PLANETGEN_TILE_CACHE_DIR=str(tmp_path / "cache" / "tiles"),
               PLANETGEN_JOBS_DIR=str(tmp_path / "lib" / "jobs"))
    subprocess.run([os.path.join(APACHE_DIR, "create-cache-dir.sh")], env=env, check=True, capture_output=True)
    for path in (tmp_path / "cache" / "tiles", tmp_path / "lib" / "jobs"):
        info = os.stat(path)
        assert (info.st_uid, info.st_gid, stat.S_IMODE(info.st_mode)) == (uid, gid, 0o750)
    assert os.stat(tmp_path / "cache").st_uid == 0


# --- apache-identity.sh ------------------------------------------------------

DEBIAN_ENVVARS = """\
unset HOME
if [ "${APACHE_CONFDIR##/etc/apache2-}" != "${APACHE_CONFDIR}" ] ; then
\tSUFFIX="-${APACHE_CONFDIR##/etc/apache2-}"
else
\tSUFFIX=
fi
export APACHE_RUN_USER=www-data
export APACHE_RUN_GROUP=www-data
"""


def test_apache_identity_reads_debians_envvars_under_set_u(tmp_path):
    """Debian's envvars reads $APACHE_CONFDIR unset; under the callers'
    `set -u` that used to kill detect_apache_group silently, and with it
    set-permissions.sh, create-cache-dir.sh and setup-debug-log.sh (TEST.62)."""
    envvars = tmp_path / "envvars"
    envvars.write_text(DEBIAN_ENVVARS)
    script = ('set -euo pipefail; source "$1"; '
              'read -r user group < <(detect_apache_group); echo "$user:$group"')
    env = dict(os.environ, APACHE_ENVVARS=str(envvars))
    env.pop("APACHE_CONFDIR", None)
    result = subprocess.run(["bash", "-c", script, "t", os.path.join(APACHE_DIR, "apache-identity.sh")],
                            env=env, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "www-data:www-data"


@pytest.mark.parametrize("name", ["set-permissions.sh", "create-cache-dir.sh", "setup-debug-log.sh"])
def test_an_unknown_apache_identity_is_an_error_not_a_silent_exit(name):
    code = _code_lines(name)
    assert "couldn't work out Apache's user and group" in code
    assert 'read -r APACHE_USER APACHE_GROUP < <(detect_apache_group) || [[ -z "${APACHE_GROUP:-}" ]]' in code
