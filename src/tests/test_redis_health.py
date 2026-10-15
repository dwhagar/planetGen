"""
ADM.54: Redis's health on the admin Stats page: `queue.redishealth.redis_info` (INFO only, short timeout, never
waits on a dead server), the API's `redis` block, and the page's row.
"""

import os

import pytest

from planetgen.queue import redisclient, redishealth
from planetgen.web import admin_pages

URL = os.environ.get("PLANETGEN_TEST_REDIS_URL")


@pytest.fixture(autouse=True)
def _fresh_clients():
    redisclient.forget_all()
    yield
    redisclient.forget_all()


def test_a_password_is_not_shown():
    assert redishealth.display_url("redis://:secret@10.0.0.5:6380/2") == "redis://10.0.0.5:6380/2"
    assert redishealth.display_url("redis://localhost:6379/0") == "redis://localhost:6379/0"
    assert redishealth.display_url("redis://user:pw@[::1]:6379/0") == "redis://[::1]:6379/0"


def test_a_dead_port_is_reported_not_raised():
    info = redishealth.redis_info("redis://:hush@127.0.0.1:1/0")
    assert info["answering"] is False and info["error"] and "hush" not in info["url"]


@pytest.mark.skipif(not URL, reason="PLANETGEN_TEST_REDIS_URL is not set")
def test_a_live_redis_reports_its_settings():
    info = redishealth.redis_info(URL)
    assert info["answering"] is True and info["error"] is None
    assert info["version"] and info["uptime_seconds"] >= 0
    assert isinstance(info["aof_enabled"], bool) and info["maxmemory_policy"]
    assert info["used_memory_bytes"] > 0 and info["maxmemory_bytes"] >= 0
    assert isinstance(info["last_save_ok"], bool)
    assert isinstance(info["workers"], int) and info["workers"] >= 0


def test_the_page_row_says_what_matters():
    live = {"url": "redis://127.0.0.1:6379/0", "answering": True, "version": "7.2.4", "uptime_seconds": 90000,
            "aof_enabled": False, "maxmemory_policy": "noeviction", "used_memory_bytes": 2 * 1024 * 1024,
            "maxmemory_bytes": 0, "last_save_ok": True, "last_save_seconds_ago": 3600, "workers": 2}
    text = admin_pages.redis_text(live)
    for fragment in ("answering at redis://127.0.0.1:6379/0", "7.2.4", "append-only file off", "noeviction",
                     "no limit", "last save ok", "2 worker(s) alive"):
        assert fragment in text, fragment
    assert "append-only file on" in admin_pages.redis_text({**live, "aof_enabled": True})
    assert "last save failed" in admin_pages.redis_text({**live, "last_save_ok": False})
    assert "no save yet" in admin_pages.redis_text({**live, "last_save_seconds_ago": None})
    assert "workers unknown" in admin_pages.redis_text({**live, "workers": None})
    assert admin_pages.redis_text({"url": "redis://h:1/0", "answering": False, "error": "refused"}) == \
        "not answering at redis://h:1/0: refused"
    assert admin_pages.redis_text(None) == "unknown"
