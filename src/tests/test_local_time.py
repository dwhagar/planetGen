# tests/test_local_time.py

"""
Times on the pages (docs/TODO.md item 22): the server writes UTC as
`<time datetime="...Z" data-local-time>` with a labelled UTC fallback,
and `static/localtime.js` rewrites it in the viewer's zone.
"""

import datetime
import os
import sys

_SRC_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_SRC_DIR, "html", "lib"))

import pytest  # noqa: E402

from fmt import utc_time_html  # noqa: E402

EXPECTED = '<time datetime="2026-09-30T21:26:41Z" data-local-time>2026-09-30 21:26 UTC</time>'


@pytest.mark.parametrize("value", [
    1790803601,
    1790803601.4,
    datetime.datetime(2026, 9, 30, 21, 26, 41),
    datetime.datetime(2026, 9, 30, 23, 26, 41, tzinfo=datetime.timezone(datetime.timedelta(hours=2))),
    "2026-09-30T21:26:41Z",
    "2026-09-30T21:26:41",
    "2026-09-30 21:26:41",
    "2026-09-30T21:26:41.123",
    "2026-09-30T17:26:41-04:00",
])
def test_every_input_becomes_the_same_utc_markup(value):
    assert utc_time_html(value) == EXPECTED


def test_missing_and_unreadable_values():
    assert utc_time_html(None) == ""
    assert utc_time_html("") == ""
    assert utc_time_html("<b>soon</b>") == "&lt;b&gt;soon&lt;/b&gt;"


def test_the_shell_loads_the_local_time_script(tmp_path):
    base = os.path.join(_SRC_DIR, "html", "web", "templates", "base.html")
    with open(base, encoding="utf-8") as handle:
        assert "static_url('localtime.js')" in handle.read()
