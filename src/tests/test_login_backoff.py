# tests/test_login_backoff.py

"""Per-username login backoff (`api/loginbackoff.py`, docs/TODO.md item
39): the counter itself with a fake clock, then the login route and the
login page against a fake `adminAuth.authenticate` (no database)."""

import pytest

from api import auth, loginbackoff
from api.loginbackoff import FREE_FAILURES, LoginBackoff


class FakeClock:
    def __init__(self):
        self.now = 1000.0

    def __call__(self):
        return self.now


@pytest.fixture
def clock():
    return FakeClock()


@pytest.fixture
def counter(clock):
    return LoginBackoff(clock=clock)


def test_free_failures_then_doubling_locks(counter, clock):
    for _ in range(FREE_FAILURES):
        assert counter.record_failure("admin") == 0
        assert counter.retry_after("admin") == 0
    assert counter.record_failure("admin") == 1
    assert counter.retry_after("admin") == 1
    clock.now += 1
    assert counter.retry_after("admin") == 0
    assert counter.record_failure("admin") == 2
    clock.now += 2
    assert counter.record_failure("admin") == 4


def test_lock_is_capped(counter):
    for _ in range(FREE_FAILURES + 60):
        seconds = counter.record_failure("admin")
    assert seconds == loginbackoff.MAX_LOCK_SECONDS
    assert counter.retry_after("admin") == int(loginbackoff.MAX_LOCK_SECONDS)


def test_usernames_are_counted_case_and_space_insensitively(counter):
    for name in ["Admin"] * FREE_FAILURES + [" admin "]:
        counter.record_failure(name)
    assert counter.retry_after("ADMIN") == 1
    assert counter.retry_after("someone-else") == 0


def test_success_clears_the_count(counter):
    for _ in range(FREE_FAILURES):
        counter.record_failure("admin")
    counter.record_success("admin")
    assert counter.record_failure("admin") == 0


def test_an_old_failure_is_forgotten(counter, clock):
    for _ in range(FREE_FAILURES):
        counter.record_failure("admin")
    clock.now += loginbackoff.FORGET_AFTER_SECONDS
    assert counter.record_failure("admin") == 0


def test_tracked_usernames_are_bounded_and_locked_ones_kept(counter, monkeypatch):
    monkeypatch.setattr(loginbackoff, "MAX_TRACKED_USERNAMES", 20)
    for _ in range(FREE_FAILURES + 1):
        counter.record_failure("target")
    for i in range(100):
        counter.record_failure(f"junk-{i}")
    assert len(counter._entries) <= 20
    assert counter.retry_after("target") >= 1


# --- The login route and page ----------------------------------------------

@pytest.fixture
def app(monkeypatch):
    from api.app import create_app
    from api.config import Config

    class _Config(Config):
        TESTING = True
        SECRET_KEY = "test"

    monkeypatch.setattr(auth, "get_control_db", lambda: object())
    good = {"id": 1, "username": "admin", "must_change_credentials": 0}

    def authenticate(_conn, username, password):
        if username == "admin" and password == "right":
            return good
        raise auth.adminAuth.AuthError("invalid username or password")

    monkeypatch.setattr(auth.adminAuth, "authenticate", authenticate)
    monkeypatch.setattr(auth.adminAuth, "create_session", lambda _conn, _id: "token")
    application = create_app(_Config)
    application.testing = True
    return application


_addresses = iter(range(1, 10 ** 6))


def _fresh_address():
    """A new client address per request, so the per-IP login limit never
    answers first: these tests are about the per-username count."""
    n = next(_addresses)
    return {"REMOTE_ADDR": f"10.{n // 65536 % 256}.{n // 256 % 256}.{n % 256}"}


def _login(client, username, password):
    return client.post("/api/auth/login", json={"username": username, "password": password},
                       environ_base=_fresh_address())


def test_route_locks_a_username_from_any_address(app):
    client = app.test_client()
    for _ in range(FREE_FAILURES):
        assert _login(client, "admin", "wrong").status_code == 401
    assert _login(client, "admin", "wrong").status_code == 401
    # Locked now: even the right password is refused, unchecked.
    response = _login(client, "admin", "right")
    assert response.status_code == 429
    # Flask-Limiter's headers raise Retry-After to its own window's reset
    # when that is later; never below the lock.
    assert int(response.headers["Retry-After"]) >= 1
    assert response.get_json()["retry_after"] == 1
    # Other usernames are unaffected.
    assert _login(client, "someone", "wrong").status_code == 401


def test_route_locks_unknown_usernames_the_same_way(app):
    client = app.test_client()
    for _ in range(FREE_FAILURES + 1):
        assert _login(client, "no-such-admin", "wrong").status_code == 401
    assert _login(client, "no-such-admin", "wrong").status_code == 429


def test_route_success_resets_the_count(app):
    client = app.test_client()
    for _ in range(FREE_FAILURES):
        _login(client, "admin", "wrong")
    assert _login(client, "admin", "right").status_code == 200
    assert _login(client, "admin", "wrong").status_code == 401
    assert _login(client, "admin", "wrong").status_code == 401


def test_login_page_names_the_wait(app):
    from web import csrf
    client = app.test_client()
    client.get("/login")
    nonce = client.get_cookie(csrf.COOKIE_NAME).value
    with app.test_request_context():
        token = csrf._sign(nonce, "")
    form = {"username": "admin", "password": "wrong", "csrf_token": token}
    for _ in range(FREE_FAILURES + 1):
        assert client.post("/login", data=form, environ_base=_fresh_address()).status_code == 200
    response = client.post("/login", data=form, environ_base=_fresh_address())
    assert response.status_code == 429
    assert "Too many failed logins for this username. Try again in 1 second." in response.get_data(as_text=True)
