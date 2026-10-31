# tests/test_deploy_scripts.py

"""
The root-run deployment scripts in `examples/apache/` (called by
`install.sh`/`update.sh`):

- `set-permissions.sh` gives the whole checkout, `config.json` (640) and
  the `state/` directory to Apache's user (www-data:www-data), the same
  way after an install and after an update.
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
import re
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
    assert _deploy_paths(tmp_path) == ["/var/cache/planetgen/tiles", "/var/lib/planetGen/jobs", "/var/lib/planetGen/state"]


def test_deploy_paths_reads_config_json(tmp_path):
    (tmp_path / "config.json").write_text(json.dumps({
        "tile_cache": {"dir": "/srv/tiles"}, "jobs": {"dir": "/srv/jobs"},
    }))
    assert _deploy_paths(tmp_path) == ["/srv/tiles", "/srv/jobs", "/var/lib/planetGen/state"]


def test_deploy_paths_env_wins_and_max_mb_0_turns_the_cache_off(tmp_path):
    (tmp_path / "config.json").write_text(json.dumps({"tile_cache": {"dir": "/srv/tiles"}, "cache": {"max_mb": 0}}))
    assert _deploy_paths(tmp_path, {"PLANETGEN_JOBS_DIR": "/env/jobs"}) == ["", "/env/jobs", "/var/lib/planetGen/state"]
    assert _deploy_paths(tmp_path, {"PLANETGEN_CACHE_MAX_MB": "5", "PLANETGEN_TILE_CACHE_DIR": "/env/t",
                                    "PLANETGEN_STATE_DIR": "/env/state"}) == [
        "/env/t", "/var/lib/planetGen/jobs", "/env/state"]


def test_deploy_paths_matches_the_web_apps_own_answer(tmp_path, monkeypatch):
    """Mirrors `tilecache.configured_cache_dir` / `jobs.configured_jobs_dir`."""
    from planetgen.web.lib import tilecache
    from planetgen.web import jobs

    from planetgen.util import logpaths
    config = {"tile_cache": {"dir": "/srv/tiles"}, "cache": {"max_mb": 50}, "jobs": {"dir": "/srv/jobs"}}
    (tmp_path / "config.json").write_text(json.dumps(config))
    monkeypatch.setattr(logpaths, "CONFIG_PATH", str(tmp_path / "config.json"))
    for var in ("PLANETGEN_TILE_CACHE_DIR", "PLANETGEN_CACHE_MAX_MB", "PLANETGEN_JOBS_DIR"):
        monkeypatch.delenv(var, raising=False)
    from planetgen.admin import maintenance
    monkeypatch.delenv("PLANETGEN_STATE_DIR", raising=False)
    assert _deploy_paths(tmp_path) == [tilecache.configured_cache_dir(), jobs.configured_jobs_dir(), maintenance.state_dir()]


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
        "/var/cache/planetgen/tiles", "/var/lib/planetGen/jobs", "/var/lib/planetGen/state"]
    assert not marker.exists()


def _deploy_paths_module():
    import importlib.util
    spec = importlib.util.spec_from_file_location("deploy_paths", DEPLOY_PATHS)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _job(root, job_id, status=None):
    path = root / job_id
    path.mkdir(parents=True)
    (path / "output.log").write_text("log of " + job_id)
    if status:
        (path / "state.json").write_text(json.dumps({"status": status}))
    return path


def test_the_jobs_default_is_inside_the_checkouts_folder():
    """OPS.19: the old default, /var/lib/planetgen/jobs, differed from the
    checkout's /var/lib/planetGen only by case: two folders."""
    from planetgen.web import jobs
    assert jobs.DEFAULT_JOBS_DIR == _deploy_paths_module().DEFAULT_JOBS_DIR == "/var/lib/planetGen/jobs"


def test_old_jobs_move_to_the_new_default_and_the_old_folders_go(tmp_path, monkeypatch):
    module = _deploy_paths_module()
    old = tmp_path / "planetgen" / "jobs"
    new = tmp_path / "planetGen" / "jobs"
    monkeypatch.setattr(module, "DEFAULT_JOBS_DIR", str(new))
    _job(old, "20261001-010101-aaaa", "succeeded")
    _job(old, "20261002-010101-bbbb", "failed")
    (old / "active").write_text("20261002-010101-bbbb")  # stale: that job finished
    lines = module.move_old_jobs(str(new), str(old))
    assert sorted(os.listdir(new)) == ["20261001-010101-aaaa", "20261002-010101-bbbb"]
    assert (new / "20261001-010101-aaaa" / "output.log").read_text() == "log of 20261001-010101-aaaa"
    assert not old.exists() and not old.parent.exists()
    assert lines[0] == f"Moved 2 Generate job(s) from {old} to {new}."


def test_a_running_old_job_stays_put(tmp_path, monkeypatch):
    module = _deploy_paths_module()
    old = tmp_path / "planetgen" / "jobs"
    new = tmp_path / "planetGen" / "jobs"
    monkeypatch.setattr(module, "DEFAULT_JOBS_DIR", str(new))
    _job(old, "20261001-010101-aaaa", "succeeded")
    _job(old, "20261002-010101-bbbb")  # still running
    (old / "active").write_text("20261002-010101-bbbb")
    lines = module.move_old_jobs(str(new), str(old))
    assert os.listdir(new) == ["20261001-010101-aaaa"]
    assert sorted(os.listdir(old)) == ["20261002-010101-bbbb", "active"]
    assert "20261002-010101-bbbb" in lines[-1]


def test_a_configured_jobs_dir_is_left_alone(tmp_path):
    module = _deploy_paths_module()
    old = tmp_path / "planetgen" / "jobs"
    _job(old, "20261001-010101-aaaa", "succeeded")
    assert module.move_old_jobs(str(tmp_path / "elsewhere"), str(old)) == []
    assert os.listdir(old) == ["20261001-010101-aaaa"]


def test_the_old_galaxy_settings_folder_is_removed_once(tmp_path):
    module = _deploy_paths_module()
    jobs = tmp_path / "jobs"
    (jobs / "galaxy-settings").mkdir(parents=True)
    (jobs / "galaxy-settings" / "a.json").write_text("{}")
    _job(jobs, "20261001-010101-aaaa", "succeeded")
    lines = module.remove_galaxy_settings(str(jobs))
    assert len(lines) == 1 and "galaxy-settings" in lines[0]
    assert os.listdir(jobs) == ["20261001-010101-aaaa"]
    assert module.remove_galaxy_settings(str(jobs)) == []  # nothing left to do


def test_a_galaxy_settings_link_is_removed_not_followed(tmp_path):
    module = _deploy_paths_module()
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    (elsewhere / "keep.txt").write_text("keep")
    jobs = tmp_path / "jobs"
    jobs.mkdir()
    os.symlink(elsewhere, jobs / "galaxy-settings")
    module.remove_galaxy_settings(str(jobs))
    assert not (jobs / "galaxy-settings").exists() and (elsewhere / "keep.txt").exists()


def test_create_cache_dir_moves_old_jobs_before_making_the_folder():
    code = _code_lines("create-cache-dir.sh")
    move = code.index('"$PYTHON" -I "$APACHE_DIR/deploy-paths.py" --move-old-jobs "$JOBS_DIR"')
    assert move < code.index('mkdir -p "$JOBS_DIR"')


def test_root_helpers_run_python_isolated():
    code = _code_lines("create-cache-dir.sh")
    assert '"$PYTHON" -I "$APACHE_DIR/deploy-paths.py"' in code
    assert "sys.path" not in code and "import tilecache" not in code and "logpaths" not in code
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
        warns = 'log_warning "' if script == "update.sh" else 'echo "warning:'      # update.sh speaks through the logging layer (UX.123)
        assert f'setup-debug-log.sh" \\\n    || {warns}' in text, script
    # Only running it without root, or with no Python at all, is an error.
    assert _code_lines("setup-debug-log.sh").count("exit 1") == 2


def test_the_update_never_offers_the_population_pass():
    # OPS.7: an update that just wiped the database mustn't offer to fill
    # it; the installers keep their prompt.
    with open(os.path.join(REPO_DIR, "update.sh"), encoding="utf-8") as f:
        shell = f.read()
    assert "offer_population_pass" not in shell and "POPULATION=" not in shell
    assert "planetgen population" in shell


# --- log-locations.py (OPS.5), run as the current user -----------------------

LOG_LOCATIONS = os.path.join(APACHE_DIR, "log-locations.py")


def _me():
    return pwd.getpwuid(os.getuid()).pw_name, grp.getgrgid(os.getgid()).gr_name


def _run_log_locations(tmp_path, user=None, group=None, config=None, **env_overrides):
    """Runs log-locations.py against a fake checkout holding the real
    logpaths.py and, when given, a config.json."""
    repo = tmp_path / "repo"
    package = repo / "src" / "planetgen" / "util"
    package.mkdir(parents=True, exist_ok=True)
    shutil.copy(os.path.join(REPO_DIR, "src", "planetgen", "util", "logpaths.py"), package / "logpaths.py")
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

needs_www_data = pytest.mark.skipif(
    os.geteuid() != 0 or _www_data() is None,
    reason="needs root and a www-data user (the Apache identity is given to the script through APACHE_ENVVARS)",
)


def _make_checkout(root):
    """A fresh `sudo git clone`: everything root's, with the layout the
    scripts expect. Returns the checkout's path."""
    repo = root / "planetGen"
    html = repo / "src" / "html"
    (html / "lib").mkdir(parents=True)
    (html / "wsgi.py").write_text("application = None\n")
    (html / "lib" / "fmt.py").write_text("")
    (html / "static").mkdir()
    (html / "static" / "style.css").write_text("")
    (repo / "scripts").mkdir()
    (repo / "scripts" / "tool.sh").write_text("#!/bin/sh\n")
    (repo / "scripts" / "helper").write_text("#!/bin/sh\n")
    os.chmod(repo / "scripts" / "helper", 0o755)
    (repo / "install.sh").write_text("#!/bin/sh\n")
    (repo / ".git" / "objects").mkdir(parents=True)
    (repo / ".git" / "HEAD").write_text("ref: refs/heads/main\n")
    os.chmod(repo / ".git" / "HEAD", 0o444)
    (repo / "db").mkdir()
    (repo / "db" / "old.db").write_text("")
    (repo / "config.json").write_text("{}")
    os.chmod(repo / "config.json", 0o644)
    return repo


def _run_permissions(repo, tmp_path, *args):
    envvars = tmp_path / "envvars"
    envvars.write_text("export APACHE_RUN_USER=www-data\nexport APACHE_RUN_GROUP=www-data\n")
    env = dict(os.environ, APACHE_ENVVARS=str(envvars))
    return subprocess.run([os.path.join(APACHE_DIR, "set-permissions.sh"), str(repo), *args],
                          env=env, capture_output=True, text=True, check=False)


def _snapshot(repo):
    """(relative path, owner, group, mode) of everything under the checkout."""
    found = []
    for root, dirs, files in os.walk(repo):
        for name in [""] + dirs + files:
            path = os.path.join(root, name)
            info = os.lstat(path)
            found.append((os.path.relpath(path, repo), info.st_uid, info.st_gid, stat.S_IMODE(info.st_mode)))
    return sorted(found)


@needs_www_data
def test_set_permissions_gives_the_whole_checkout_to_the_web_user(tmp_path):
    uid, gid = _www_data()
    repo = _make_checkout(tmp_path)
    outside = tmp_path / "outside.txt"
    outside.write_text("not planetGen's")
    (repo / "link").symlink_to(outside)
    proc = _run_permissions(repo, tmp_path)
    assert proc.returncode == 0, proc.stderr
    for rel, owner, group, mode in _snapshot(repo):
        if rel == "link":
            continue
        assert (owner, group) == (uid, gid), rel
        assert not mode & stat.S_IRWXO, rel
        if rel != "state":  # the group may write the state folder, nothing else
            assert not mode & stat.S_IWGRP, rel
    modes = {rel: mode for rel, _, _, mode in _snapshot(repo)}
    assert modes["."] == 0o750 and modes["src/html/lib"] == 0o750
    assert modes["src/html/wsgi.py"] == 0o750 and modes["scripts/tool.sh"] == 0o750
    assert modes["scripts/helper"] == 0o750  # keeps its owner execute bit, so git sees no change
    assert modes["src/html/static/style.css"] == 0o640
    assert modes[".git/HEAD"] == 0o640
    assert modes["config.json"] == 0o640
    # The state directory is made, and is the web user's from the start.
    assert (modes["state"], (repo / "state").is_dir()) == (0o770, True)
    # A link is never followed: what it points at is not planetGen's to chown.
    assert os.stat(outside).st_uid == 0


@needs_www_data
def test_set_permissions_takes_a_state_directory_outside_the_checkout(tmp_path):
    uid, gid = _www_data()
    repo = _make_checkout(tmp_path)
    state = tmp_path / "elsewhere" / "state"
    proc = _run_permissions(repo, tmp_path, str(state))
    assert proc.returncode == 0, proc.stderr
    info = os.stat(state)
    assert (info.st_uid, info.st_gid, stat.S_IMODE(info.st_mode)) == (uid, gid, 0o770)
    assert not (repo / "state").exists()


@needs_www_data
def test_a_fresh_install_and_an_update_leave_identical_ownership(tmp_path):
    """install.sh and update.sh both end in set-permissions.sh on the whole
    checkout, so it does not matter which one touched the tree last."""
    for name in ("install.sh", "update.sh"):
        with open(os.path.join(REPO_DIR, name), encoding="utf-8") as f:
            assert '"$SCRIPT_DIR/examples/apache/set-permissions.sh" "$SCRIPT_DIR"' in f.read(), name

    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    fresh = _make_checkout(tmp_path / "a")
    updated = _make_checkout(tmp_path / "b")
    # `updated` was installed earlier, then an update pulled new and changed
    # files as root and left a root-owned bytecode cache behind.
    assert _run_permissions(updated, tmp_path / "b").returncode == 0
    for repo in (fresh, updated):
        (repo / "src" / "html" / "new.py").write_text("")
        (repo / "src" / "html" / "lib" / "fmt.py").write_text("changed")
        os.chmod(repo / "src" / "html" / "lib" / "fmt.py", 0o644)
    (updated / "src" / "html" / "__pycache__").mkdir()
    (updated / "src" / "html" / "__pycache__" / "wsgi.pyc").write_text("")
    (fresh / "src" / "html" / "__pycache__").mkdir()
    (fresh / "src" / "html" / "__pycache__" / "wsgi.pyc").write_text("")
    shutil.rmtree(updated / "state")  # the install tree can lose it; the next run makes it again
    assert _run_permissions(fresh, tmp_path / "a").returncode == 0
    assert _run_permissions(updated, tmp_path / "b").returncode == 0
    assert _snapshot(fresh) == _snapshot(updated)


@needs_www_data
def test_set_permissions_refuses_a_folder_that_is_not_a_checkout(tmp_path):
    (tmp_path / "keep").write_text("")
    proc = _run_permissions(tmp_path, tmp_path)
    assert proc.returncode != 0
    assert "doesn't look like a planetGen checkout" in proc.stderr
    assert os.stat(tmp_path / "keep").st_uid == 0


def test_the_update_trusts_its_own_checkout_though_the_web_user_owns_it():
    """git refuses a repository owned by another user ("dubious ownership"),
    and the checkout now belongs to www-data while update.sh runs as root."""
    with open(os.path.join(REPO_DIR, "update.sh"), encoding="utf-8") as f:
        text = f.read()
    export = 'export GIT_CONFIG_COUNT=1 GIT_CONFIG_KEY_0=safe.directory GIT_CONFIG_VALUE_0="$SCRIPT_DIR"'
    assert export in text
    assert text.index(export) < text.index("git rev-parse --is-inside-work-tree")


@needs_www_data
def test_git_as_root_accepts_a_web_user_checkout_with_the_updates_setting(tmp_path):
    uid, gid = _www_data()
    repo = tmp_path / "repo"
    repo.mkdir()
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    for root, dirs, files in os.walk(repo):
        for name in [""] + dirs + files:
            os.chown(os.path.join(root, name), uid, gid)
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_CONFIG")}
    env["HOME"] = str(tmp_path)  # no global safe.directory
    env["GIT_CONFIG_NOSYSTEM"] = "1"
    plain = subprocess.run(["git", "rev-parse", "--is-inside-work-tree"], cwd=repo, env=env,
                           capture_output=True, text=True)
    if plain.returncode == 0:
        pytest.skip("this git does not check repository ownership")
    env.update(GIT_CONFIG_COUNT="1", GIT_CONFIG_KEY_0="safe.directory", GIT_CONFIG_VALUE_0=str(repo))
    trusted = subprocess.run(["git", "rev-parse", "--is-inside-work-tree"], cwd=repo, env=env,
                             capture_output=True, text=True)
    assert trusted.returncode == 0, trusted.stderr


@needs_root
def test_create_cache_dir_gives_apache_only_the_runtime_dirs(tmp_path):
    uid, gid = _www_data()
    env = {k: v for k, v in os.environ.items() if not k.startswith("PLANETGEN_")}
    env.update(PLANETGEN_TILE_CACHE_DIR=str(tmp_path / "cache" / "tiles"),
               PLANETGEN_JOBS_DIR=str(tmp_path / "lib" / "jobs"),
               PLANETGEN_STATE_DIR=str(tmp_path / "lib" / "state"))
    subprocess.run([os.path.join(APACHE_DIR, "create-cache-dir.sh")], env=env, check=True, capture_output=True)
    for path in (tmp_path / "cache" / "tiles", tmp_path / "lib" / "jobs"):
        info = os.stat(path)
        assert (info.st_uid, info.st_gid, stat.S_IMODE(info.st_mode)) == (uid, gid, 0o750)
    # The state folder (the maintenance switch) is the web user's to write, the group's too, never the world's.
    info = os.stat(tmp_path / "lib" / "state")
    assert (info.st_uid, info.st_gid, stat.S_IMODE(info.st_mode)) == (uid, gid, 0o770)
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


# --- reload_web_if_running (OPS.8, OPS.29) -----------------------------------

def _run_reload(tmp_path, active=("apache2",), restart=0, fail=(), configtest=True, health=True, reply="", curl=None, site=None):
    """Runs the helper with fake `systemctl`, `apache2ctl` and `curl` that log their calls; returns
    `(process, calls)`. `fail` names the systemctl verbs that fail. `curl` is the body of the fake `curl` (it sees the
    address as its last argument) when the test needs different answers for different addresses; `site` is what a fake
    python prints for `site_health_target` (the site's own address and its two --resolve values)."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    log = tmp_path / "calls"

    def fake(name, body):
        path = bin_dir / name
        path.write_text("#!/bin/bash\n" + f'echo "{name} $@" >> {log}\n' + body)
        path.chmod(0o755)

    fake("systemctl", 'if [[ "$1" == is-active ]]; then case " %s " in *" $3 "*) exit 0;; esac; exit 3; fi\n'
         'case " %s " in *" $1 "*) exit 1;; esac\nexit 0\n' % (" ".join(active), " ".join(fail)))
    fake("apache2ctl", "exit %d\n" % (0 if configtest else 1))
    fake("curl", curl or "printf '%%s\\n%%s' '%s' %d\n" % (reply, 200 if health else 503))
    if site is not None:
        fake_python = bin_dir / "fakepython"
        fake_python.write_text(f"#!/bin/bash\necho '{site}'\n")
        fake_python.chmod(0o755)
    script = (f'source "{REPO_DIR}/scripts/deploy-common.sh"; APACHE_NEEDS_RESTART={restart}; '
              'reload_web_if_running')
    env = dict(os.environ, PATH=f"{bin_dir}:{os.environ['PATH']}", HEALTH_WAIT_SECONDS="1", HEALTH_SETTLE_SECONDS="0")
    if site is not None:
        env["PYTHON"] = str(bin_dir / "fakepython")
    proc = subprocess.run(["bash", "-c", script], env=env, capture_output=True, text=True)
    calls = log.read_text().splitlines() if log.exists() else []
    return proc, calls


def test_the_update_reloads_a_running_apache_after_its_config_test_and_checks_the_site_answers(tmp_path):
    proc, calls = _run_reload(tmp_path)
    assert proc.returncode == 0 and "reloaded Apache" in proc.stdout
    assert [c for c in calls if not c.startswith("curl")] == [
        "systemctl is-active --quiet apache2", "apache2ctl configtest", "systemctl reload apache2",
        "systemctl is-active --quiet planetgen-gunicorn"]
    assert any(c.startswith("curl") and "/api/health" in c for c in calls)


def test_the_update_restarts_apache_only_when_mod_wsgi_was_just_installed(tmp_path):
    proc, calls = _run_reload(tmp_path, restart=1)
    assert proc.returncode == 0 and "restarted Apache" in proc.stdout
    assert "systemctl restart apache2" in calls and "systemctl reload apache2" not in calls


def test_the_update_leaves_a_stopped_apache_alone_and_says_nothing_ran(tmp_path):
    proc, calls = _run_reload(tmp_path, active=())
    assert proc.returncode == 2
    assert not any(c.startswith(("apache2ctl", "systemctl reload", "systemctl restart", "curl")) for c in calls)


def test_a_failed_config_test_leaves_the_running_apache_alone(tmp_path):
    proc, calls = _run_reload(tmp_path, configtest=False)
    assert proc.returncode == 1 and "configuration test failed" in proc.stderr
    assert "systemctl reload apache2" not in calls


def test_a_failed_reload_is_a_failure_for_the_caller(tmp_path):
    proc, _calls = _run_reload(tmp_path, fail=("reload",))
    assert proc.returncode == 1 and "failed" in proc.stderr


def test_a_gunicorn_host_with_no_apache_reloads_the_unit(tmp_path):
    proc, calls = _run_reload(tmp_path, active=("planetgen-gunicorn",))
    assert proc.returncode == 0 and "reloaded gunicorn" in proc.stdout
    assert "systemctl reload planetgen-gunicorn" in calls and "apache2ctl configtest" not in calls
    assert any(c.startswith("curl") and ":8000/api/health" in c for c in calls)


def test_a_reload_that_the_site_does_not_answer_after_is_a_warning_not_a_failed_reload(tmp_path):
    """Every reload step worked and only the probe got no good answer: status 3, which update.sh turns into a warning
    and an open site (the migration is done by then), not status 1 (a step failed, the site stays closed)."""
    proc, _calls = _run_reload(tmp_path, health=False)
    assert proc.returncode == 3 and "health check of Apache failed" in proc.stderr


MAINTENANCE_REPLY = '{"status":"ok","maintenance":{"cause":"update","enabled":true,"reason":"Updating planetGen"}}'


def test_a_reload_check_that_finds_the_site_closed_for_the_update_still_passes_and_says_why(tmp_path):
    proc, _calls = _run_reload(tmp_path, reply=MAINTENANCE_REPLY)
    assert proc.returncode == 0
    assert "in maintenance mode (cause: update, Updating planetGen)" in proc.stdout


def test_a_failed_health_check_names_the_check_the_answer_and_the_maintenance_state(tmp_path):
    proc, _calls = _run_reload(tmp_path, health=False, reply='{"maintenance":true,"error":"closed"}')
    assert proc.returncode == 3
    assert "health check of Apache failed" in proc.stderr and "HTTP 503" in proc.stderr
    assert "http://127.0.0.1/api/health" in proc.stderr and "The site is in maintenance mode" in proc.stderr


SITE = "https://starmap.example/api/health starmap.example:80:127.0.0.1 starmap.example:443:127.0.0.1"
OTHER_VHOST = ("for last; do :; done\ncase \"$last\" in\n"
               "  http://127.0.0.1/*) printf '<html>403</html>\\n403';;\n"
               "  https://starmap.example/*) printf '{\"status\":\"ok\"}\\n%s' \"$SITE_CODE\";;\nesac\n")


def test_a_bare_loopback_address_that_belongs_to_another_virtual_host_is_checked_through_the_sites_own_name(tmp_path):
    """OPS.77: Apache refused http://127.0.0.1/api/health with its own 403 (another virtual host answers the bare
    address) while the site's name answered; the check asks the site's name, pointed at this machine."""
    os.environ["SITE_CODE"] = "200"
    try:
        proc, calls = _run_reload(tmp_path, curl=OTHER_VHOST, site=SITE)
    finally:
        del os.environ["SITE_CODE"]
    assert proc.returncode == 0, proc.stderr
    assert "answered HTTP 403" in proc.stdout and "the site answers at https://starmap.example/api/health" in proc.stdout
    site_calls = [c for c in calls if c.startswith("curl") and "starmap.example/api/health" in c]
    assert site_calls and "--resolve starmap.example:443:127.0.0.1" in site_calls[0] and "--resolve starmap.example:80:127.0.0.1" in site_calls[0]


def test_a_failed_check_names_both_addresses_when_neither_answers(tmp_path):
    os.environ["SITE_CODE"] = "503"
    try:
        proc, _calls = _run_reload(tmp_path, curl=OTHER_VHOST, site=SITE)
    finally:
        del os.environ["SITE_CODE"]
    assert proc.returncode == 3
    assert "http://127.0.0.1/api/health answered HTTP 403" in proc.stderr
    assert "own address https://starmap.example/api/health answered HTTP 503" in proc.stderr


def test_the_sites_own_name_is_not_taken_for_an_answer_unless_it_is_a_2xx(tmp_path):
    os.environ["SITE_CODE"] = "301"
    try:
        proc, _calls = _run_reload(tmp_path, curl=OTHER_VHOST, site=SITE)
    finally:
        del os.environ["SITE_CODE"]
    assert proc.returncode == 3 and "answered HTTP 301" in proc.stderr


NO_ANSWER = ("for last; do :; done\ncase \"$last\" in\n"
             "  http://127.0.0.1/*) printf '<html>403</html>\\n403';;\n"
             "  https://starmap.example/*) printf '\\n000'; exit 7;;\nesac\n")


def test_an_https_probe_with_no_http_answer_says_why_and_the_reload_still_counts_as_done(tmp_path):
    """A live update ended with "answered HTTP 000" and the site closed. 000 is no status code: curl got no
    answer. The warning now names the reason, and the update opens the migrated site."""
    proc, _calls = _run_reload(tmp_path, curl=NO_ANSWER, site=SITE)
    assert proc.returncode == 3
    assert "answered HTTP 000 (no HTTP answer at all: the connection was refused" in proc.stderr
    assert "reloaded Apache" in proc.stdout


def test_the_site_probe_does_not_fail_on_a_certificate_the_loopback_check_cannot_judge(tmp_path):
    os.environ["SITE_CODE"] = "200"
    try:
        _proc, calls = _run_reload(tmp_path, curl=OTHER_VHOST, site=SITE)
    finally:
        del os.environ["SITE_CODE"]
    assert any(c.startswith("curl -skL") and "starmap.example/api/health" in c for c in calls)


def test_update_sh_opens_the_site_when_only_the_health_probe_failed():
    update = (pathlib.Path(REPO_DIR) / "update.sh").read_text()
    branch = update.split("elif (( reload_status == 3 ))", 1)[1].split("else", 1)[0]
    assert "log_warning" in branch and "exit" not in branch


def test_the_dirty_tree_warning_leaves_out_the_state_folder_the_update_owns():
    update = (pathlib.Path(REPO_DIR) / "update.sh").read_text()
    assert "status --porcelain -- . ':(exclude)state'" in update
    assert 'UPDATE_STEP="reloading what serves the site' in update


def test_the_update_ends_non_zero_and_prints_the_commands_when_the_reload_failed():
    update = (pathlib.Path(REPO_DIR) / "update.sh").read_text()
    assert "reload_web_if_running || reload_status=$?" in update and "print_reload_commands >&2\n        exit 1" in update


# --- scripts only reach for files that exist --------------------------------

def test_the_log_location_check_runs_against_this_checkout(tmp_path):
    """The inline Python in setup-debug-log.sh, run exactly as the script runs it
    (isolated, from `/`, loading the settings code by its path): it loaded the
    deleted util/appconfig.py and every update printed "couldn't work out where
    the logs go" while nothing in the suite noticed."""
    text = _read("setup-debug-log.sh")
    code = text.split("<<'PY'\n", 1)[1].split("\nPY\n", 1)[0]
    env = {k: v for k, v in os.environ.items() if not k.startswith("PLANETGEN_")}
    result = subprocess.run([sys.executable, "-I", "-", REPO_DIR], input=code, cwd="/", env=env,
                            capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    first, second = result.stdout.splitlines()
    assert first.split(" ", 1)[0] in ("0", "1")
    assert first.split(" ", 1)[1].endswith(".log") and second.endswith(".log")


def _tracked_scripts():
    tracked = subprocess.run(["git", "ls-files", "*.sh", "*.py"], cwd=REPO_DIR, capture_output=True,
                             text=True, check=False).stdout.split()
    return [p for p in tracked if p in ("install.sh", "update.sh") or p.startswith(("scripts/", "examples/"))]


def test_tracked_shell_scripts_are_executable_in_git():
    """A fresh `git clone` gives the mode git tracks, and INSTALL.md says `sudo ./install.sh`: with
    install.sh tracked as 100644 that fails with "Permission denied" on a new server (the executable
    bit only survives a clone if git has it)."""
    listing = subprocess.run(["git", "ls-files", "-s", "*.sh"], cwd=REPO_DIR, capture_output=True,
                             text=True, check=False).stdout.splitlines()
    if not listing:
        pytest.skip("not a git checkout")
    not_executable = []
    for line in listing:
        mode, _sha, _stage, path = line.replace("\t", " ").split(None, 3)
        with open(os.path.join(REPO_DIR, path), encoding="utf-8") as f:
            runnable = f.readline().startswith("#!")
        if runnable and mode != "100755":
            not_executable.append(f"{path} is tracked as {mode}")
    assert not not_executable, "\n".join(not_executable)


def test_scripts_only_name_repo_paths_that_exist():
    """Any `src/...`, `scripts/...`, `examples/...`, `db/...` file a deploy script names, and any
    os.path.join(..., "src", "planetgen", ...) it builds, must exist in the checkout (a moved or
    deleted module left setup-debug-log.sh pointing at nothing)."""
    import re

    scripts = _tracked_scripts()
    if not scripts:
        pytest.skip("not a git checkout")
    literal = re.compile(r"""(?<![\w./-])((?:src|scripts|examples|db)/[\w./-]*\.(?:py|sh|sql|json|conf|plist|service|timer)(?:\.example)?)(?![\w-])""")
    joined = re.compile(r'"src"((?:,\s*"[\w.-]+")+)')
    missing = []
    for rel in scripts:
        with open(os.path.join(REPO_DIR, rel), encoding="utf-8") as f:
            text = f.read()
        for m in literal.finditer(text):
            if "*" not in m.group(1) and not os.path.exists(os.path.join(REPO_DIR, m.group(1))):
                missing.append(f"{rel}: {m.group(1)}")
        for m in joined.finditer(text):
            parts = re.findall(r'"([\w.-]+)"', m.group(1))
            path = os.path.join("src", *parts)
            if not os.path.exists(os.path.join(REPO_DIR, path)):
                missing.append(f"{rel}: {path}")
    assert not missing, "\n".join(missing)


# --- example configs the admin owns -------------------------------------------

def _example_repo(tmp_path):
    repo = tmp_path / "repo"
    (repo / "examples" / "apache").mkdir(parents=True)
    (repo / "examples" / "redis").mkdir(parents=True)
    (repo / "examples" / "apache" / "planetgen.conf.example").write_text("one\n")
    (repo / "examples" / "redis" / "planetgen.conf").write_text("one\n")
    git = lambda *a: subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", *a], cwd=repo,
                                    check=True, capture_output=True, text=True).stdout.strip()
    git("init", "-q")
    git("add", "-A")
    git("commit", "-q", "-m", "one")
    return repo, git


def _run_note(repo, before, after):
    script = (f'SCRIPT_DIR="{repo}"; source "{REPO_DIR}/scripts/deploy-common.sh"; '
              f'cd "{repo}"; note_changed_example_configs "$1" "$2"')
    return subprocess.run(["bash", "-c", script, "t", before, after], capture_output=True, text=True)


def test_the_update_says_when_an_example_config_changed(tmp_path):
    repo, git = _example_repo(tmp_path)
    before = git("rev-parse", "HEAD")
    (repo / "examples" / "apache" / "planetgen.conf.example").write_text("two\n")
    git("commit", "-qam", "two")
    after = git("rev-parse", "HEAD")
    proc = _run_note(repo, before, after)
    assert proc.returncode == 0, proc.stderr
    assert "planetgen.conf.example changed in this update" in proc.stdout
    assert "was not touched" in proc.stdout
    assert "Redis" not in proc.stdout


def test_the_update_stays_quiet_when_no_example_config_changed(tmp_path):
    repo, git = _example_repo(tmp_path)
    before = git("rev-parse", "HEAD")
    (repo / "other.txt").write_text("x")
    git("add", "-A")
    git("commit", "-qm", "other")
    assert _run_note(repo, before, git("rev-parse", "HEAD")).stdout == ""
    assert _run_note(repo, before, before).stdout == ""


def test_the_scripts_never_write_the_admins_apache_or_redis_config():
    """The vhost and the Redis config are the admin's: only the examples in
    the repo change, and update.sh says so (note_changed_example_configs)."""
    for rel in ("install.sh", "update.sh", "scripts/deploy-common.sh", "examples/apache/set-permissions.sh"):
        with open(os.path.join(REPO_DIR, rel), encoding="utf-8") as f:
            code = "\n".join(line for line in f.read().splitlines() if not line.lstrip().startswith("#"))
        for target in ("sites-available/planetgen.conf", "sites-enabled", "/etc/redis"):
            for line in code.splitlines():
                # `sudo ...` and echo lines are instructions printed for the admin to run.
                if target in line and not line.lstrip().startswith(("sudo ", "echo ", "diff ")):
                    assert not re.search(r"(>|\bsed -i\b|\bcp\b|\bmv\b|\btee\b|\binstall\b|\bln\b)[^|]*" + re.escape(target), line), (rel, line)


@needs_www_data
def test_an_update_that_started_before_the_ownership_change_is_told_to_run_again(tmp_path):
    """The update.sh that was running when it pulled the change still hands
    src/html to set-permissions.sh, which now takes the checkout."""
    repo = _make_checkout(tmp_path)
    proc = _run_permissions(repo / "src" / "html", tmp_path)
    assert proc.returncode != 0
    assert "run ./update.sh again" in proc.stderr
    assert os.stat(repo / "src" / "html" / "wsgi.py").st_uid == 0


def test_the_updates_dirty_check_ignores_mode_changes():
    """set-permissions.sh makes every script executable; that is not a local edit to warn about."""
    with open(os.path.join(REPO_DIR, "update.sh"), encoding="utf-8") as f:
        assert 'git -c core.fileMode=false status --porcelain' in f.read()


def test_the_update_works_out_the_cluster_system_masses_after_the_migration():
    """The web daemon is stopped at 60 seconds and building the galaxy shape needs a mean system mass for every stored
    cluster age, which takes longer on a big galaxy: the update does it once and the web app reads the result."""
    update = (pathlib.Path(REPO_DIR) / "update.sh").read_text()
    assert update.index("migrate_or_reset_db\nrecord_version_key\nwarm_shape_means") < update.index("log_stage 6 9")
    assert "run_cli warm_shape" in (pathlib.Path(REPO_DIR) / "scripts" / "deploy-common.sh").read_text()


def test_the_update_starts_the_warm_up_early_and_reloads_after_its_opening_phase():
    """The warm-up starts in the background once everything it needs is in place; the reload and the health check wait
    for its opening phase, then the probe waits for the workers to start."""
    update = (pathlib.Path(REPO_DIR) / "update.sh").read_text()
    common = (pathlib.Path(REPO_DIR) / "scripts" / "deploy-common.sh").read_text()
    assert update.index("\nwarm_galaxy_map\n") < update.index("\n    wait_for_warm_opening\n") < update.index("reload_web_if_running ||")
    assert update.index("check_app_imports") < update.index("\nwarm_galaxy_map\n")
    assert "--opening-done" in common
    assert 'sleep "${HEALTH_SETTLE_SECONDS:-10}"' in common


def test_wait_for_warm_opening_returns_when_the_file_appears_or_the_warm_up_ends(tmp_path):
    common = pathlib.Path(REPO_DIR) / "scripts" / "deploy-common.sh"
    script = f"""
source {common}
WARM_OPENING_FILE={tmp_path}/opening
sleep 30 & WARM_PID=$!
( sleep 1; touch "$WARM_OPENING_FILE" ) &
wait_for_warm_opening
kill $WARM_PID
rm -f "$WARM_OPENING_FILE"
sleep 0.2 & WARM_PID=$!
wait $WARM_PID
WARM_OPENING_WAIT_SECONDS=5 wait_for_warm_opening
"""
    result = subprocess.run(["bash", "-c", script], capture_output=True, text=True, timeout=60)
    assert "The opening phase is built" in result.stdout
    assert "ended before its opening phase" in result.stderr


def test_every_wait_in_the_update_says_what_it_waits_for_and_for_how_long():
    """Boss, 10-30 EDT: a long step or a wait that prints nothing looks like a hang."""
    common = (pathlib.Path(REPO_DIR) / "scripts" / "deploy-common.sh").read_text()
    for words in ("Working out the cluster system masses", "Running the math check", "Importing the web app as",
                  "at most $(( ${WARM_OPENING_WAIT_SECONDS:-1200} / 60 )) minutes", "Giving the workers",
                  "asking again,", "Testing Apache's configuration, then"):
        assert words in common, words
