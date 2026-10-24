# tests/test_maintenance.py

"""
Maintenance mode (UX.112): the state file (`admin/maintenance.py`), the
command (`cli/maintenance.py`), and the request gate (`web/maintenance_gate.py`)
that turns every page and API call into a "being updated" answer for anyone
but an admin. The gate tests build the real app but need no database: a
refused visitor never reaches one, and the admin check is faked.
"""

import contextlib
import grp
import json
import os
import pwd
import shutil
import subprocess
import tempfile

import pytest

from planetgen.admin import maintenance
from planetgen.api.config import Config
from planetgen.api.version import API_VERSION_HEADER
from planetgen.cli import maintenance as maintenance_cli


@pytest.fixture
def state_file(tmp_path, monkeypatch):
    path = tmp_path / "maintenance.json"
    monkeypatch.setenv(maintenance.FILE_ENV, str(path))
    return path


def test_no_file_means_off(state_file):
    assert maintenance.read() == {"enabled": False, "eta": "", "reason": "", "cause": "", "since": None, "by": None}
    assert not maintenance.is_on()


def test_turn_on_and_off_round_trip(state_file):
    maintenance.turn_on("  about   two hours ", by="boss")
    state = maintenance.read()
    assert state["enabled"] and state["eta"] == "about two hours" and state["by"] == "boss"
    assert state["since"].endswith("Z")
    assert json.loads(state_file.read_text())["enabled"] is True
    maintenance.turn_off(by="boss")
    assert not maintenance.is_on() and maintenance.read()["eta"] == ""


def test_the_eta_has_a_limit(state_file):
    with pytest.raises(ValueError, match="at most 200"):
        maintenance.turn_on("x" * 201)
    assert not state_file.exists()


@pytest.mark.parametrize("text", ["{not json", "[]", ""])
def test_a_file_that_cannot_be_read_counts_as_on(state_file, text):
    state_file.write_text(text)
    assert maintenance.is_on()


def test_a_changed_file_is_seen_without_a_restart(state_file):
    maintenance.turn_on("later")
    assert maintenance.read()["eta"] == "later"
    state_file.write_text(json.dumps({"enabled": True, "eta": "tomorrow, then", "since": None, "by": None}))
    assert maintenance.read()["eta"] == "tomorrow, then"


def test_the_sector_switch_is_cleared_only_by_its_own_cause(state_file):
    maintenance.turn_on("later", by="boss", cause="")
    assert not maintenance.turn_off_for(maintenance.SECTOR_SWITCH) and maintenance.is_on()
    maintenance.turn_on(by="migration", reason=maintenance.SECTOR_SWITCH_REASON, cause=maintenance.SECTOR_SWITCH)
    state = maintenance.read()
    assert state["reason"] == maintenance.SECTOR_SWITCH_REASON and state["cause"] == maintenance.SECTOR_SWITCH
    # An admin adding an expected finish keeps the reason and the cause.
    maintenance.turn_on("about two hours", by="boss")
    state = maintenance.read()
    assert state["eta"] == "about two hours" and state["cause"] == maintenance.SECTOR_SWITCH
    assert maintenance.turn_off_for(maintenance.SECTOR_SWITCH, by="map rebuild")
    assert not maintenance.is_on() and not maintenance.turn_off_for(maintenance.SECTOR_SWITCH)


@pytest.mark.parametrize("version, closes", [(0, False), (61, True), (90, True), (91, False), (92, False)])
def test_migrating_across_the_sector_switch_closes_the_site(state_file, monkeypatch, version, closes):
    from planetgen.cli import migrate

    monkeypatch.setattr(migrate, "schema_status", lambda config: (version, 1))
    monkeypatch.setattr(migrate, "SCHEMA_VERSION", 93)
    lines = []
    assert migrate.close_site_for_the_sector_switch(None, write=lines.append) is closes
    assert maintenance.is_on() is closes
    if closes:
        assert maintenance.read()["cause"] == maintenance.SECTOR_SWITCH and "only admins" in lines[0]


def test_migrating_leaves_an_admins_maintenance_alone(state_file, monkeypatch):
    from planetgen.cli import migrate

    monkeypatch.setattr(migrate, "schema_status", lambda config: (90, 1))
    monkeypatch.setattr(migrate, "SCHEMA_VERSION", 93)
    maintenance.turn_on("later", by="boss")
    assert not migrate.close_site_for_the_sector_switch(None, write=lambda line: None)
    assert maintenance.read()["cause"] == "" and maintenance.read()["eta"] == "later"


def test_the_command_switches_it(state_file, capsys):
    assert maintenance_cli.main(["status"]) == 0
    assert "off" in capsys.readouterr().out
    assert maintenance_cli.main(["on", "--eta", "by Thursday"]) == 0
    out = capsys.readouterr().out
    assert "ON" in out and "by Thursday" in out
    assert maintenance.is_on()
    assert maintenance_cli.main(["off"]) == 0
    assert not maintenance.is_on()
    assert maintenance_cli.main(["on", "--eta", "x" * 300]) == 1


# --- The web user writes the switch (the Permission denied on /var/lib/planetGen, 2026-10-15) --------------------

APACHE_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "examples", "apache")


def _web_user():
    try:
        return pwd.getpwnam("www-data").pw_uid, grp.getgrnam("www-data").gr_gid
    except KeyError:
        return None


needs_root = pytest.mark.skipif(os.geteuid() != 0 or _web_user() is None,
                                reason="needs root and a www-data user to act as the web server's user")


@contextlib.contextmanager
def _as_web_user():
    """Everything inside runs with the web user's file permissions (not root's, which ignore them)."""
    uid, gid = _web_user()
    os.setegid(gid)
    os.seteuid(uid)
    try:
        yield
    finally:
        os.seteuid(0)
        os.setegid(0)


@pytest.fixture
def install_root():
    """A root-owned folder the web user can enter, like /var/lib/planetGen."""
    base = tempfile.mkdtemp(prefix="pg-maint-")
    os.chmod(base, 0o755)
    yield base
    shutil.rmtree(base, ignore_errors=True)


def test_the_default_state_file_is_in_the_state_folder_not_the_install_folder(monkeypatch):
    monkeypatch.delenv(maintenance.FILE_ENV, raising=False)
    monkeypatch.delenv(maintenance.STATE_DIR_ENV, raising=False)
    assert maintenance.state_path() == "/var/lib/planetGen/state/maintenance.json"
    monkeypatch.setenv(maintenance.STATE_DIR_ENV, "/srv/state")
    assert maintenance.state_path() == "/srv/state/maintenance.json"


@needs_root
def test_a_switch_in_a_root_owned_folder_says_what_to_do_and_the_admin_page_warns(install_root, monkeypatch):
    monkeypatch.delenv(maintenance.FILE_ENV, raising=False)
    monkeypatch.setenv(maintenance.STATE_DIR_ENV, install_root)
    maintenance.turn_on(by="migration", cause=maintenance.SECTOR_SWITCH)  # root: the migration
    with _as_web_user():
        assert maintenance.is_on(), "the web user can read it"
        assert "not writable" in maintenance.write_problem()
        with pytest.raises(maintenance.MaintenanceWriteError, match="update.sh") as raised:
            maintenance.turn_off(by="boss")
        assert install_root in str(raised.value) and "Permission denied" in str(raised.value)
        with pytest.raises(maintenance.MaintenanceWriteError):
            maintenance.turn_off_for(maintenance.SECTOR_SWITCH, by="map rebuild")
    assert maintenance.is_on()


@needs_root
def test_the_web_user_turns_it_on_and_off_in_the_folder_create_cache_dir_makes(install_root, monkeypatch):
    """What `update.sh` does: the migration (root) switches it on, `create-cache-dir.sh` gives the folder to
    Apache's user, then the admin page and `planetgen warm` (the web user) switch it off."""
    state = os.path.join(install_root, "state")
    monkeypatch.delenv(maintenance.FILE_ENV, raising=False)
    monkeypatch.setenv(maintenance.STATE_DIR_ENV, state)
    maintenance.turn_on(by="migration", reason=maintenance.SECTOR_SWITCH_REASON, cause=maintenance.SECTOR_SWITCH)
    assert os.stat(state).st_uid == 0, "the migration, run as root, made the folder"
    env = {k: v for k, v in os.environ.items() if not k.startswith("PLANETGEN_")}
    env.update(PLANETGEN_STATE_DIR=state, PLANETGEN_TILE_CACHE_DIR=os.path.join(install_root, "tiles"),
               PLANETGEN_JOBS_DIR=os.path.join(install_root, "jobs"))
    subprocess.run([os.path.join(APACHE_DIR, "create-cache-dir.sh")], env=env, check=True, capture_output=True)
    with _as_web_user():
        assert maintenance.write_problem() is None
        assert maintenance.turn_off_for(maintenance.SECTOR_SWITCH, by="map rebuild")
        assert not maintenance.is_on()
        maintenance.turn_on("about an hour", by="boss")
        maintenance.turn_off(by="boss")
    assert not maintenance.is_on()


def test_the_command_says_when_it_cannot_write_the_switch(state_file, monkeypatch, capsys):
    monkeypatch.setattr(maintenance.tempfile, "mkstemp", _denied)
    assert maintenance_cli.main(["on"]) == 1
    assert "update.sh" in capsys.readouterr().err


def _denied(*args, **kwargs):
    raise PermissionError(13, "Permission denied")


def test_the_map_rebuild_that_cannot_reopen_the_site_fails_loudly(state_file, monkeypatch, caplog, redis_server):
    from planetgen.cli import generate as generate_cli
    from planetgen.cli import warm as warm_cli
    from planetgen.queue import commands
    from planetgen.web import warmup

    maintenance.turn_on(by="migration", cause=maintenance.SECTOR_SWITCH)
    monkeypatch.setattr(maintenance.tempfile, "mkstemp", _denied)
    monkeypatch.setattr(warm_cli, "_site_database", lambda: "x")
    monkeypatch.setattr(warmup, "warm_family", lambda db, family, cap=None, app=None, config=None: {
        "family": family, "warmed": 1, "left_out": 0, "seconds": 0.0, "line": f"{family}: ok"})
    monkeypatch.setattr(warmup, "expand_warm_galaxy_tiles", lambda params: commands.Expansion(
        [({"database": params["database"], "cap": params.get("cap")}, 1.0)], {}))
    monkeypatch.setattr(warmup, "warm_galaxy_tiles_unit", lambda args, report: warmup._warm_unit("galaxy-tiles", args))
    monkeypatch.setattr(warmup, "expand_warm_galaxy_nodes", lambda params: commands.Expansion(
        [({"database": params["database"], "cap": params.get("cap")}, 1.0)], {}))
    monkeypatch.setattr(warmup, "warm_galaxy_nodes_unit", lambda args, report: warmup._warm_unit("galaxy-nodes", args))
    monkeypatch.setattr(warmup, "expand_warm_galaxy_deep_tiles", lambda params: commands.Expansion(
        [({"database": params["database"], "cap": params.get("cap")}, 1.0)], {}))
    monkeypatch.setattr(warmup, "warm_galaxy_deep_tiles_unit",
                        lambda args, report: warmup._warm_unit("galaxy-deep-tiles", args))
    monkeypatch.setattr(warmup, "warm_object_index_unit", lambda args, report: warmup._warm_unit("object-index", args))
    parser, _parsers = generate_cli.build_parser()
    with pytest.raises(SystemExit) as exit_info:
        warm_cli.run_warm(parser.parse_args(["warm"]))
    assert exit_info.value.code == 1
    assert "MaintenanceWriteError" in caplog.text and "update.sh" in caplog.text
    assert maintenance.is_on()


class _Config(Config):
    SESSION_COOKIE_SECURE = False
    SECRET_KEY = "test-secret"


@pytest.fixture
def app(state_file):
    from planetgen.web.app import create_app

    application = create_app(_Config)
    application.config["TESTING"] = True
    return application


@pytest.fixture
def client(app):
    return app.test_client()


@pytest.fixture
def as_admin(monkeypatch):
    from planetgen.api import authz

    monkeypatch.setattr(authz, "_current_admin", lambda: {"username": "boss"})


def test_with_maintenance_off_the_gate_does_nothing(client):
    assert client.get("/maintenance").status_code == 302   # sent home
    assert client.get("/static/style.css").status_code == 200


def test_a_visitor_is_sent_to_the_maintenance_page(client, state_file):
    maintenance.turn_on("about two hours")
    for path in ("/", "/galaxy", "/search?q=a", "/sector/1", "/admin"):
        response = client.get(path)
        assert response.status_code == 302 and response.headers["Location"].endswith("/maintenance"), path
    page = client.get("/maintenance")
    html = page.get_data(as_text=True)
    assert page.status_code == 503 and page.headers["Retry-After"] and page.headers["Cache-Control"] == "no-store"
    assert "being updated" in html and "about two hours" in html and 'content="noindex"' in html
    maintenance.turn_on(reason=maintenance.SECTOR_SWITCH_REASON, cause=maintenance.SECTOR_SWITCH)
    assert maintenance.SECTOR_SWITCH_REASON in client.get("/maintenance").get_data(as_text=True)
    assert client.post("/admin", data={"action": "x"}).status_code == 302


def test_the_page_says_when_there_is_no_estimate(client):
    maintenance.turn_on()
    assert "no estimate" in client.get("/maintenance").get_data(as_text=True)


def test_the_api_answers_503_json(client):
    maintenance.turn_on("about two hours")
    for path in ("/api/galaxy/sectors", "/api/systems"):
        response = client.get(path)
        assert response.status_code == 503, path
        body = response.get_json()
        assert body["maintenance"] is True and body["eta"] == "about two hours" and "updated" in body["error"]
        assert response.headers["Retry-After"] and response.headers[API_VERSION_HEADER]
    assert client.post("/api/systems", json={}).status_code == 503


@pytest.mark.parametrize("path", ["/login", "/login/code", "/logout", "/maintenance", "/static/style.css",
                                  "/api/auth/login", "/api/auth/login/totp", "/api/auth/logout", "/api/auth/me"])
def test_the_sign_in_paths_and_static_files_pass_the_gate(app, path):
    from planetgen.web import maintenance_gate

    maintenance.turn_on()
    with app.test_request_context(path):
        assert maintenance_gate.gate() is None


def test_an_admin_passes_the_gate_everywhere(app, as_admin):
    from planetgen.web import maintenance_gate

    maintenance.turn_on("later")
    for path in ("/", "/galaxy", "/admin", "/api/galaxy/sectors"):
        with app.test_request_context(path):
            assert maintenance_gate.gate() is None, path


def test_a_control_database_that_does_not_answer_counts_as_a_visitor(app, monkeypatch):
    from planetgen.api import authz
    from planetgen.web import maintenance_gate

    def broken():
        raise RuntimeError("no control database")

    monkeypatch.setattr(authz, "_current_admin", broken)
    maintenance.turn_on()
    with app.test_request_context("/api/galaxy/sectors"):
        assert maintenance_gate.gate().status_code == 503


def test_the_pages_own_api_calls_skip_the_gate(app):
    from planetgen.api.limiter import IN_PROCESS_ENVIRON_KEY
    from planetgen.web import maintenance_gate

    maintenance.turn_on()
    with app.test_request_context("/api/galaxy/sectors", environ_base={IN_PROCESS_ENVIRON_KEY: True}):
        assert maintenance_gate.gate() is None


def test_the_maintenance_page_for_an_admin_is_200(client, as_admin):
    maintenance.turn_on("later")
    page = client.get("/maintenance")
    assert page.status_code == 200 and "later" in page.get_data(as_text=True)


def test_the_galaxy_map_redirects_while_maintenance_is_on(client, as_admin):
    maintenance.turn_on()
    response = client.get("/galaxy")
    assert response.status_code == 302 and response.headers["Location"].endswith("/maintenance")
