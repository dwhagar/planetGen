# tests/test_redis_settings_check.py

"""
OPS.48: `scripts/check_redis_settings.py` warns, with the line to add, about Redis settings that don't suit planetGen,
and the example config sets every one of them. It never changes a setting.
"""

import importlib.util
import os

ROOT = os.path.realpath(os.path.join(os.path.dirname(__file__), "..", ".."))
_spec = importlib.util.spec_from_file_location("check_redis_settings", os.path.join(ROOT, "scripts", "check_redis_settings.py"))
check = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(check)

GOOD = {"appendonly": "yes", "appendfsync": "everysec", "maxmemory": "536870912", "maxmemory-policy": "noeviction",
        "bind": "127.0.0.1 -::1"}


def _lines(**changes):
    return [line for _what, line in check.problems(dict(GOOD, **changes))]


def test_the_settings_planetgen_wants_raise_no_warning():
    assert check.problems(GOOD) == []


def test_redis_defaults_warn_about_each_setting_with_the_line_to_add():
    defaults = {"appendonly": "no", "appendfsync": "everysec", "maxmemory": "0", "maxmemory-policy": "noeviction",
                "bind": "* -::*"}
    lines = [line for _what, line in check.problems(defaults)]
    assert "appendonly yes\nappendfsync everysec" in lines
    assert "maxmemory 512mb" in lines
    assert "bind 127.0.0.1 -::1\nprotected-mode yes" in lines
    assert len(lines) == 3


def test_each_setting_warns_alone():
    assert _lines(appendonly="no") == ["appendonly yes\nappendfsync everysec"]
    assert _lines(**{"maxmemory-policy": "allkeys-lru"}) == ["maxmemory-policy noeviction"]
    assert _lines(appendfsync="no") == ["appendfsync everysec"]
    assert _lines(maxmemory="0") == ["maxmemory 512mb"]


def test_only_loopback_addresses_count_as_loopback():
    for bind in ("127.0.0.1", "127.0.0.1 -::1", "::1", "localhost", "127.0.0.2"):
        assert not check.listens_beyond_loopback(bind), bind
    for bind in ("", "*", "0.0.0.0", "* -::*", "127.0.0.1 192.168.1.5", "::"):
        assert check.listens_beyond_loopback(bind), bind


def test_the_example_config_sets_what_the_check_wants():
    with open(os.path.join(ROOT, "examples", "redis", "planetgen.conf"), encoding="utf-8") as f:
        settings = dict(line.split(None, 1) for line in f.read().splitlines() if line and not line.startswith(("#", " ")))
    assert settings["appendonly"] == "yes" and settings["appendfsync"] == "everysec"
    assert settings["maxmemory-policy"] == "noeviction" and settings["maxmemory"] == "512mb"
    assert not check.listens_beyond_loopback(settings["bind"])


def test_an_unreachable_server_only_warns(capsys):
    assert check.main(["check", "redis://127.0.0.1:1/0"]) == 0
    assert "couldn't read the Redis settings" in capsys.readouterr().err
