# tests/test_system_load.py

"""
The load line of the admin Queue page (`planetgen/queue/load.py`).
"""

import os

from planetgen.queue import load as systemLoad


def test_the_load_line_on_this_system():
    load = systemLoad.load_average()
    if hasattr(os, "getloadavg"):
        assert load["kind"] == "load" and len(load["values"]) == 3
        assert load["text"].count(" / ") == 2


def test_a_system_without_a_load_average_says_unknown(monkeypatch):
    monkeypatch.delattr(os, "getloadavg", raising=False)
    assert systemLoad.load_average() == {"values": None, "kind": "load", "text": "unknown"}
