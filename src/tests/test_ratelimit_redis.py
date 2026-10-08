# tests/test_ratelimit_redis.py

"""
The request limits and login lockouts count on Redis (SEC.30): two apps
(two worker processes) on one server share the counts, a server that
can't be reached drops to counting in memory instead of failing requests,
and an install with no `ratelimit.storage_uri` set gets the Redis server
of `redis.url`.
"""

import pytest

from planetgen.admin import throttle
from planetgen.api import config as api_config, loginguard
from planetgen.api.config import Config
from planetgen.db.store import MySQLConfig
from planetgen.web.app import create_app

_UNREACHABLE = MySQLConfig(host="127.0.0.1", port=1, user="secretuser", password="x", database="planetgen_x")


def _app(storage, prefix, **settings):
    attributes = dict(MYSQL_CONFIG=_UNREACHABLE, WRITE_MYSQL_CONFIG=_UNREACHABLE,
                      CONTROL_MYSQL_CONFIG=_UNREACHABLE, RATELIMIT_STORAGE_URI=storage,
                      RATELIMIT_KEY_PREFIX=prefix, SESSION_COOKIE_SECURE=False, SECRET_KEY="test-secret")
    attributes.update(settings)
    application = create_app(type("TestConfig", (Config,), attributes))
    application.config["PROPAGATE_EXCEPTIONS"] = False
    return application


def _health(client, ip="10.0.0.1"):
    return client.get("/api/health", environ_base={"REMOTE_ADDR": ip}).status_code


def test_two_apps_on_one_redis_share_the_request_counts(redis_server, key_prefix):
    pages = {"health": "4 per minute"}
    first = _app(redis_server, key_prefix, RATELIMIT_PAGES=pages).test_client()
    second = _app(redis_server, key_prefix, RATELIMIT_PAGES=pages).test_client()
    statuses = [_health(first), _health(second), _health(first), _health(second)]
    assert 429 not in statuses
    assert _health(first) == 429 and _health(second) == 429
    assert _health(second, ip="10.0.0.2") != 429


def test_an_unreachable_redis_counts_in_memory_and_still_answers():
    client = _app("redis://127.0.0.1:1/0", "planetgen-test:down", RATELIMIT_PAGES={"health": "2 per minute"}).test_client()
    statuses = [_health(client) for _ in range(3)]
    assert 500 not in statuses and statuses[-1] == 429


def test_the_login_lockouts_follow_the_rate_limit_storage(redis_server, key_prefix):
    application = _app(redis_server, key_prefix)
    with application.test_request_context("/"):
        store = loginguard.redis_store()
        assert isinstance(store, throttle.RedisStore) and store.key_prefix == f"{key_prefix}:login:"
        loginguard.with_store(lambda s: throttle.record_failure(s, throttle.SCOPE_USER, "someone"))
        assert store.get(throttle.SCOPE_USER, "someone")["failures"] == 1
        assert loginguard.memory_store.get(throttle.SCOPE_USER, "someone") is None
    in_memory = _app("memory://", key_prefix)
    with in_memory.test_request_context("/"):
        assert loginguard.redis_store() is None


def test_lockouts_count_in_memory_while_redis_is_down(capsys):
    application = _app("redis://127.0.0.1:1/0", "planetgen-test:down")
    with application.test_request_context("/"):
        loginguard.with_store(lambda s: throttle.record_failure(s, throttle.SCOPE_USER, "someone"))
        assert loginguard.memory_store.get(throttle.SCOPE_USER, "someone")["failures"] == 1
    assert "Redis can't be used for the login lockouts" in capsys.readouterr().err


def test_no_configured_storage_means_the_redis_server(monkeypatch):
    monkeypatch.delenv("PLANETGEN_RATELIMIT_STORAGE_URI", raising=False)
    monkeypatch.delenv("PLANETGEN_REDIS_URL", raising=False)
    config = {"ratelimit": {"storage_uri": ""}, "redis": {"url": "redis://cache.example:6380/3"}}
    assert api_config.ratelimit_storage_uri(config) == "redis://cache.example:6380/3"
    config["ratelimit"]["storage_uri"] = "memory://"
    assert api_config.ratelimit_storage_uri(config) == "memory://"
    config["ratelimit"]["storage_uri"] = ""
    monkeypatch.setenv("PLANETGEN_REDIS_URL", "redis://other:6379/1")
    assert api_config.ratelimit_storage_uri(config) == "redis://other:6379/1"
    monkeypatch.setenv("PLANETGEN_RATELIMIT_STORAGE_URI", "memory://")
    assert api_config.ratelimit_storage_uri(config) == "memory://"
