# tests/test_activity_log.py

"""
The always-on activity log (SEC.28, `planetgen/admin/activity_log.py`) and
the sign-in lines and audit rows SEC.20 writes through it.
"""

import os
import re
import stat

import pytest

from planetgen.admin import activity_log, auth as adminAuth
from planetgen.util import appconfig

LINE = re.compile(
    r'^\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ planetgen\[\d+\]: (?P<cat>[A-Z]+) (?P<action>\S+) '
    r'ip=(?P<ip>\S+) user="(?P<user>(?:[^"\\]|\\.)*)"(?P<rest>.*)$'
)
"""The documented line shape (docs/config.md)."""

FAIL2BAN = re.compile(r'^\s*planetgen\[\d+\]: AUTH (?:login\.failed|login\.locked|password\.failed) '
                      r'ip=(?P<host>[0-9a-fA-F.:]+)(?:\s|$)')
"""SEC.27's example failregex, with `<HOST>` spelled out; fail2ban cuts the
date off the front before matching."""


@pytest.fixture
def log_dir(tmp_path, monkeypatch):
    folder = tmp_path / "logs"
    monkeypatch.setenv("PLANETGEN_LOG_DIR", str(folder))
    monkeypatch.setattr(appconfig, "SYSTEM_ROTATION_FILES", (str(tmp_path / "no-such-rotation"),))
    activity_log.reset()
    yield folder
    activity_log.reset()


def _lines(folder):
    path = folder / "planetgen.log"
    return path.read_text(encoding="utf-8").splitlines() if path.exists() else []


# --- line format --------------------------------------------------------------

def test_line_shape():
    line = activity_log.format_line("AUTH", "login.failed", ip="203.0.113.5", user="admin",
                                   when=1790000000, pid=1234)
    assert line == '2026-09-21T14:13:20Z planetgen[1234]: AUTH login.failed ip=203.0.113.5 user="admin"'
    assert LINE.match(line)
    assert FAIL2BAN.match(line.split(" ", 1)[1]).group("host") == "203.0.113.5"


def test_missing_ip_and_user_are_dashes():
    line = activity_log.format_line("GEN", "generate.start", fields={"command": "galaxy", "db": "planetgen"})
    match = LINE.match(line)
    assert match.group("ip") == "-" and match.group("user") == "-"
    assert match.group("rest") == " command=galaxy db=planetgen"


@pytest.mark.parametrize("address", ["1.2.3.4 user=x", "not-an-address", "", None, "1.2.3.4\n", "::1 ip=9.9.9.9"])
def test_bad_addresses_never_reach_the_line_unquoted(address):
    line = activity_log.format_line("AUTH", "login.failed", ip=address, user="a")
    assert LINE.match(line).group("ip") in ("-", "1.2.3.4")
    assert "user=x" not in line and "9.9.9.9" not in line


def test_ipv6_is_normalized():
    line = activity_log.format_line("AUTH", "login.failed", ip="2001:DB8:0:0::1", user="a")
    assert " ip=2001:db8::1 " in line


@pytest.mark.parametrize("username", [
    'admin" ip=1.1.1.1',
    "admin\n2026-10-01T00:00:00Z planetgen[1]: AUTH login.failed ip=6.6.6.6 user=\"x\"",
    "back\\slash",
    "tab\there\r",
    "line\u2028sep\x85",
    "x" * 1000,
])
def test_crafted_usernames_cant_fake_fields_or_lines(username):
    line = activity_log.format_line("AUTH", "login.failed", ip="10.0.0.1", user=username)
    assert "\n" not in line and "\r" not in line and "\u2028" not in line and "\x85" not in line
    match = LINE.match(line)
    assert match and match.group("ip") == "10.0.0.1" and match.group("rest") == ""
    assert len(line) < 400
    assert FAIL2BAN.match(line.split(" ", 1)[1]).group("host") == "10.0.0.1"


def test_field_values_quoted_only_when_needed():
    line = activity_log.format_line("DB", "sector.update", user="boss", fields={
        "target": "sector:42", "detail": "{'name': 'New \"Home\"'}", "count": 3, "seconds": 1.5,
        "ok": True, "skip": None,
    })
    assert line.endswith(' target=sector:42 detail="{\'name\': \'New \\"Home\\"\'}" count=3 seconds=1.5 ok=yes')


@pytest.mark.parametrize("category, action, fields", [
    ("NOPE", "x", {}), ("AUTH", "Bad Action", {}), ("AUTH", "x", {"ip": "1"}), ("AUTH", "x", {"Bad": 1}),
])
def test_bad_calls_are_refused(category, action, fields):
    with pytest.raises(ValueError):
        activity_log.format_line(category, action, fields=fields)


# --- the file -----------------------------------------------------------------

def test_event_writes_one_line_to_the_file(log_dir):
    line = activity_log.event("AUTH", "logout", user="boss", ip="192.0.2.7")
    assert _lines(log_dir) == [line]
    assert activity_log.log_path() == str(log_dir / "planetgen.log")
    activity_log.event("DB", "migrate", db="planetgen", from_version=43, to_version=44)
    lines = _lines(log_dir)
    assert len(lines) == 2 and lines[1].endswith("DB migrate ip=- user=\"-\" db=planetgen from_version=43 to_version=44")


@pytest.mark.skipif(os.name != "posix", reason="POSIX file modes")
def test_new_file_is_not_world_readable(log_dir):
    activity_log.event("AUTH", "logout", user="boss")
    mode = stat.S_IMODE(os.stat(log_dir / "planetgen.log").st_mode)
    assert mode & 0o007 == 0


def test_app_rotation_without_system_rotation(log_dir):
    activity_log.event("AUTH", "logout")
    handler = activity_log._handler
    assert handler.maxBytes == activity_log.ROTATE_MAX_BYTES
    assert handler.backupCount == activity_log.ROTATE_BACKUP_COUNT


def test_system_rotation_reopens_the_moved_file(log_dir, tmp_path, monkeypatch):
    marker = tmp_path / "planetgen-log"
    marker.write_text("")
    monkeypatch.setattr(appconfig, "SYSTEM_ROTATION_FILES", (str(marker),))
    activity_log.reset()
    activity_log.event("AUTH", "logout", user="before")
    os.rename(log_dir / "planetgen.log", log_dir / "planetgen.log.1")
    activity_log.event("AUTH", "logout", user="after")
    assert len(_lines(log_dir)) == 1 and 'user="after"' in _lines(log_dir)[0]


def test_unwritable_folder_warns_once_and_never_raises(tmp_path, monkeypatch, capsys):
    blocker = tmp_path / "file-not-folder"
    blocker.write_text("")
    monkeypatch.setenv("PLANETGEN_LOG_DIR", str(blocker / "logs"))
    activity_log.reset()
    try:
        assert activity_log.event("AUTH", "logout") is not None
        assert activity_log.event("AUTH", "logout") is not None
        assert activity_log.log_path() is None
    finally:
        activity_log.reset()


def test_existing_but_unwritable_folder_warns(tmp_path, monkeypatch, capsys):
    folder = tmp_path / "logs"
    folder.mkdir()
    monkeypatch.setenv("PLANETGEN_LOG_DIR", str(folder))
    activity_log.reset()

    def refuse(*args, **kwargs):
        raise PermissionError(13, "Permission denied")
    monkeypatch.setattr(activity_log, "_RotatingHandler", refuse)
    monkeypatch.setattr(activity_log, "_WatchedHandler", refuse)
    try:
        activity_log.event("AUTH", "logout")
        activity_log.event("AUTH", "logout")
        err = capsys.readouterr().err
        assert err.count("can't be opened") == 1
    finally:
        activity_log.reset()


# --- configuration ------------------------------------------------------------

def test_default_folders_per_platform():
    assert appconfig.default_log_dir("linux") == "/var/log/planetgen"
    assert appconfig.default_log_dir("darwin") == "/Library/Logs/planetgen"
    assert appconfig.default_log_dir("win32") == os.path.join(
        os.path.dirname(appconfig.CONFIG_PATH), "logs")


def test_log_dir_precedence(monkeypatch):
    monkeypatch.delenv("PLANETGEN_LOG_DIR", raising=False)
    assert appconfig.log_dir_path({"log_dir": ""}) == appconfig.default_log_dir()
    assert appconfig.log_dir_path({"log_dir": "/srv/logs"}) == "/srv/logs"
    monkeypatch.setenv("PLANETGEN_LOG_DIR", "/env/logs")
    assert appconfig.log_dir_path({"log_dir": "/srv/logs"}) == "/env/logs"


def test_never_shares_the_debug_logs_file(monkeypatch, tmp_path):
    monkeypatch.delenv("PLANETGEN_LOG_DIR", raising=False)
    monkeypatch.delenv("PLANETGEN_LOG_FILE", raising=False)
    config = {"log_dir": str(tmp_path), "log_file": str(tmp_path / "planetgen.log")}
    assert appconfig.activity_log_path(config) == str(tmp_path / "planetgen-activity.log")
    config["log_file"] = str(tmp_path / "debug.log")
    assert appconfig.activity_log_path(config) == str(tmp_path / "planetgen.log")


def test_rotation_mode(monkeypatch, tmp_path):
    marker = tmp_path / "planetgen-log"
    monkeypatch.setattr(appconfig, "SYSTEM_ROTATION_FILES", (str(marker),))
    assert appconfig.log_rotation_mode({"log_rotation": "system"}) == "system"
    assert appconfig.log_rotation_mode({"log_rotation": " APP "}) == "app"
    if not appconfig.sys.platform.startswith("win"):
        assert appconfig.log_rotation_mode({"log_rotation": "auto"}) == "app"
        marker.write_text("")
        assert appconfig.log_rotation_mode({}) == "system"


def test_debug_log_also_gets_each_line(log_dir, monkeypatch):
    from planetgen.util import log
    seen = []
    monkeypatch.setattr(log, "trace", lambda message, *args, **kwargs: seen.append(message % args))
    line = activity_log.event("AUTH", "logout", user="boss")
    assert seen == [f"activity: {line}"]


# --- sign-ins through the API (SEC.20) ----------------------------------------

@pytest.fixture
def real_app(mysql_config, log_dir):
    from api.app import create_app
    from api.config import Config

    class RealConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config
        SESSION_COOKIE_SECURE = False
        SECRET_KEY = "test-secret"

    _username, password = adminAuth.bootstrap_control_schema(mysql_config)
    adminAuth._db.get_connection(mysql_config).close()  # the content schema
    application = create_app(RealConfig)
    application.testing = True
    application.first_password = password
    application.mysql_config = mysql_config
    return application


def _login(client, username, password, ip="198.51.100.9"):
    return client.post("/api/auth/login", json={"username": username, "password": password},
                       environ_base={"REMOTE_ADDR": ip})


def test_failed_and_good_logins_are_logged_with_the_address(real_app, log_dir):
    client = real_app.test_client()
    assert _login(client, "admin", "wrong", ip="203.0.113.5").status_code == 401
    assert _login(client, "admin", real_app.first_password).status_code == 200
    assert client.post("/api/auth/logout", environ_base={"REMOTE_ADDR": "198.51.100.9"}).status_code == 200
    events = [LINE.match(line).groupdict() for line in _lines(log_dir)]
    assert [(e["cat"], e["action"], e["ip"], e["user"]) for e in events] == [
        ("AUTH", "login.failed", "203.0.113.5", "admin"),
        ("AUTH", "login.ok", "198.51.100.9", "admin"),
        ("AUTH", "logout", "198.51.100.9", "admin"),
    ]
    assert "wrong" not in "".join(_lines(log_dir)) and real_app.first_password not in "".join(_lines(log_dir))


def test_failures_go_to_the_audit_log_and_the_api(real_app):
    client = real_app.test_client()
    _login(client, "nobody", "wrong", ip="203.0.113.5")
    conn = adminAuth._db.get_control_connection(real_app.mysql_config)
    try:
        rows = adminAuth.recent_login_failures(conn)
    finally:
        conn.close()
    assert [(r["action"], r["username"], r["ip"]) for r in rows] == [("login.failed", "nobody", "203.0.113.5")]

    password = "a-new-long-password-1"
    assert _login(client, "admin", real_app.first_password).status_code == 200
    assert client.post("/api/auth/change-credentials", json={
        "current_password": real_app.first_password, "new_username": "admin", "new_password": password,
    }).status_code == 200
    body = client.get("/api/admin/login-failures").get_json()
    assert [(i["action"], i["username"], i["ip"]) for i in body["items"]] == [
        ("login.failed", "nobody", "203.0.113.5")]
    assert body["items"][0]["created_at"].endswith("Z")
    assert client.get("/api/admin/login-failures?limit=0").status_code == 400


def test_wrong_current_password_is_logged(real_app, log_dir):
    client = real_app.test_client()
    assert _login(client, "admin", real_app.first_password).status_code == 200
    response = client.post("/api/auth/change-credentials", json={
        "current_password": "guess", "new_username": "admin", "new_password": "a-new-long-password-1",
    }, environ_base={"REMOTE_ADDR": "198.51.100.77"})
    assert response.status_code == 400
    assert any(" AUTH password.failed ip=198.51.100.77 " in line for line in _lines(log_dir))


def test_old_failure_rows_are_pruned(real_app):
    conn = adminAuth._db.get_control_connection(real_app.mysql_config)
    try:
        conn.execute("INSERT INTO admin_audit_log (admin_username, action, created_at) "
                     "VALUES ('old', 'login.failed', DATE_SUB(CURRENT_TIMESTAMP, INTERVAL 91 DAY))")
        conn.execute("INSERT INTO admin_audit_log (admin_username, action, created_at) "
                     "VALUES ('boss', 'sector.create', DATE_SUB(CURRENT_TIMESTAMP, INTERVAL 400 DAY))")
        conn.commit()
        adminAuth.record_login_failure(conn, "login.failed", "new", ip="192.0.2.1")
        names = [r["admin_username"] for r in conn.execute("SELECT admin_username FROM admin_audit_log").fetchall()]
    finally:
        conn.close()
    assert "old" not in names and "new" in names and "boss" in names


def test_refused_requests_are_logged(real_app, log_dir):
    client = real_app.test_client()
    client.get("/api/auth/me")  # "is anyone logged in?" -- not logged
    client.get("/api/auth/api-keys", environ_base={"REMOTE_ADDR": "192.0.2.10"})
    client.get("/api/auth/me", headers={"Authorization": "Bearer pg_nope"})
    client.set_cookie("pg_admin_session", "forged")
    client.get("/api/auth/me")
    actions = [LINE.match(line).group("action") for line in _lines(log_dir)]
    assert actions == ["login.required", "apikey.invalid", "session.invalid"]
    assert ' path=/api/auth/api-keys' in _lines(log_dir)[0] and " ip=192.0.2.10 " in _lines(log_dir)[0]


def test_database_writes_are_logged(real_app, log_dir):
    client = real_app.test_client()
    assert _login(client, "admin", real_app.first_password).status_code == 200
    assert client.post("/api/auth/change-credentials", json={
        "current_password": real_app.first_password, "new_username": "boss",
        "new_password": "a-new-long-password-1",
    }).status_code == 200
    response = client.post("/api/sectors", json={"name": "Logged Sector", "edge_ly": 20})
    assert response.status_code == 201, response.get_data(as_text=True)
    sector_id = response.get_json()["id"]
    lines = _lines(log_dir)
    assert any(" AUTH credentials.changed " in line and 'user="boss" old_user=admin' in line for line in lines)
    created = [line for line in lines if " DB sector.create " in line]
    assert len(created) == 1 and f"target=sector:{sector_id}" in created[0] and 'user="boss"' in created[0]
