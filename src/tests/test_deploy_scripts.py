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

The scripts that only touch a directory they're given run for real
against a throwaway tree when the tests run as root and a `www-data`
user exists (the scripts' fallback Apache identity); otherwise those
tests are skipped.
"""

import grp
import json
import os
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


# --- setup-debug-log.sh ------------------------------------------------------

def test_debug_log_is_0660_never_world_writable():
    code = _code_lines("setup-debug-log.sh")
    assert "0666" not in code
    assert 'chmod 0660 "$LOG_FILE"' in code
    assert "create 0660 $APACHE_USER $APACHE_GROUP" in code


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
