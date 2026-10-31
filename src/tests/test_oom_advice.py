"""PERF.223: `scripts/oom_advice.py` says how to keep the kernel's out-of-memory killer off MySQL and Redis and offers a
memory limit on Apache's unit (where the workers run); it only prints."""

import importlib.util
import os

ROOT = os.path.realpath(os.path.join(os.path.dirname(__file__), "..", ".."))
_spec = importlib.util.spec_from_file_location("oom_advice", os.path.join(ROOT, "scripts", "oom_advice.py"))
oom = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(oom)


def test_systemctl_output_is_read_into_properties():
    assert oom.parse_show("LoadState=loaded\nOOMScoreAdjust=-800\n\nMemoryMax=infinity") == {
        "LoadState": "loaded", "OOMScoreAdjust": "-800", "MemoryMax": "infinity"}


def test_a_unit_needs_the_drop_in_only_when_loaded_and_above_minus_900():
    assert oom.needs_adjust({"LoadState": "loaded", "OOMScoreAdjust": "0"})
    assert oom.needs_adjust({"LoadState": "loaded", "OOMScoreAdjust": "-800"})          # the old INSTALL.md advice is not enough
    assert not oom.needs_adjust({"LoadState": "loaded", "OOMScoreAdjust": "-900"})
    assert not oom.needs_adjust({"LoadState": "loaded", "OOMScoreAdjust": "-1000"})
    assert not oom.needs_adjust({"LoadState": "not-found", "OOMScoreAdjust": "0"})
    assert not oom.needs_adjust({"LoadState": "loaded", "OOMScoreAdjust": "oops"})


def test_the_advice_gives_each_unprotected_unit_its_drop_in_and_one_restart():
    lines = oom.advice(["mysql", "redis-server"], apache_loaded=False, apache_limited=False, memory_mb=0)
    text = "\n".join(lines)
    assert "OOMScoreAdjust=-900" in text and "/etc/systemd/system/mysql.service.d/oom.conf" in text
    assert "/etc/systemd/system/redis-server.service.d/oom.conf" in text
    assert text.count("systemctl daemon-reload") == 1 and "restart mysql redis-server" in text
    assert "apache2" not in text


def test_apache_gets_an_optional_limit_above_the_scatter_stages_half_of_memory():
    assert oom.apache_limits(15_992) == (8, 10)                      # 55 and 65 percent of a 16 GB host
    assert oom.apache_limits(1024) == (1, 2)
    text = "\n".join(oom.advice([], apache_loaded=True, apache_limited=False, memory_mb=15_992))
    assert "MemoryHigh=8G" in text and "MemoryMax=10G" in text and text.startswith("optional")
    assert oom.advice([], apache_loaded=True, apache_limited=True, memory_mb=15_992) == []      # already limited
    assert oom.advice([], apache_loaded=False, apache_limited=False, memory_mb=15_992) == []


def test_a_protected_host_with_a_limited_apache_prints_nothing():
    assert oom.advice([], True, True, 15_992) == []
