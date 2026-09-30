"""
html/lib/pagination.py regression tests -- the site's one shared pager.

Same `sys.path` setup as `test_navmap.py` (`html/lib` isn't part of the
installed `stellarObjects` package); no database needed.

Run with: pytest src/tests/test_pagination.py
"""
import os
import re
import sys

_SRC_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_SRC_DIR, "html", "lib"))

import pytest  # noqa: E402

from pagination import (  # noqa: E402
    PAGE_SIZE,
    _page_numbers,
    clamp_page,
    fetch_page,
    page_count,
    page_slice,
    parse_page,
    render_pagination,
)


def test_page_size_is_fifty():
    assert PAGE_SIZE == 50


@pytest.mark.parametrize("raw, expected", [
    (None, 1), ("", 1), ("abc", 1), ("0", 1), ("-3", 1), ("1", 1), ("7", 7), (4, 4),
])
def test_parse_page_falls_back_to_page_one(raw, expected):
    assert parse_page(raw) == expected


def test_page_count_and_clamp():
    assert page_count(0) == 1
    assert page_count(50) == 1
    assert page_count(51) == 2
    assert page_count(1234) == 25
    assert clamp_page(99, 120) == 3
    assert clamp_page(0, 120) == 1


def test_page_slice_clamps_a_stale_page_to_the_last_one():
    items = list(range(120))
    rows, page = page_slice(items, 2)
    assert (rows[0], rows[-1], len(rows), page) == (50, 99, 50, 2)
    rows, page = page_slice(items, 9)
    assert (rows, page) == (list(range(100, 120)), 3)
    assert page_slice([], 4) == ([], 1)


def test_fetch_page_asks_for_one_page_and_refetches_past_the_end():
    calls = []

    def fetch(limit, offset):
        calls.append((limit, offset))
        return {"items": list(range(offset, min(offset + limit, 120))), "total": 120}

    envelope, page = fetch_page(fetch, 2)
    assert (page, calls, envelope["items"][0]) == (2, [(50, 50)], 50)

    calls.clear()
    envelope, page = fetch_page(fetch, 10)
    assert page == 3
    assert calls == [(50, 450), (50, 100)]
    assert envelope["items"][0] == 100


def test_page_numbers_window_with_gaps():
    assert _page_numbers(1, 1) == [1]
    assert _page_numbers(1, 5) == [1, 2, 3, 4, 5]
    assert _page_numbers(1, 20) == [1, 2, 3, None, 20]
    assert _page_numbers(10, 20) == [1, None, 8, 9, 10, 11, 12, None, 20]
    # A one-page gap shows that page rather than "...".
    assert _page_numbers(5, 20) == [1, 2, 3, 4, 5, 6, 7, None, 20]
    assert _page_numbers(20, 20) == [1, None, 18, 19, 20]


def test_render_pagination_is_empty_when_everything_fits():
    assert render_pagination("browse.py", {"db": "x"}, "page", 1, 50) == ""


def _posted_pages(html, page_param="sectors_page"):
    return [int(v) for v in re.findall(rf'name="{page_param}" value="(\d+)"', html)]


def test_render_pagination_middle_page():
    html = render_pagination(
        "browse.py", {"db": "galaxy", "sectors_page": 10, "standalone_page": 3}, "sectors_page", 10, 1000,
        anchor="sectors", label="Sector pages",
    )
    assert html.startswith('<nav class="pagination" aria-label="Sector pages">')
    assert "Showing 451&ndash;500 of 1,000" in html
    # First, Prev, the window, Next, Last -- each posting its own page number.
    assert _posted_pages(html) == [1, 9, 1, 8, 9, 11, 12, 20, 11, 20]
    assert '<span class="page-link current" aria-current="page">10</span>' in html
    assert html.count("&hellip;") == 2
    # Every link keeps the other table's page and lands on the table.
    assert html.count('name="standalone_page" value="3"') == 10
    assert html.count('name="db" value="galaxy"') == 10
    assert 'action="browse.py#sectors"' in html


def test_render_pagination_disables_controls_at_either_end():
    first = render_pagination("phenomena.py", {"db": "x"}, "page", 1, 120)
    assert first.count('aria-disabled="true"') == 2
    assert "&laquo; First</span>" in first and "&lsaquo; Prev</span>" in first
    assert "Showing 1&ndash;50 of 120" in first

    last = render_pagination("phenomena.py", {"db": "x"}, "page", 3, 120)
    assert last.count('aria-disabled="true"') == 2
    assert "Next &rsaquo;</span>" in last and "Last &raquo;</span>" in last
    assert "Showing 101&ndash;120 of 120" in last


def test_render_pagination_keeps_repeated_params():
    params = [("db", "x"), ("class", "M"), ("class", "K"), ("stars_page", 2)]
    html = render_pagination("search.py", params, "stars_page", 2, 200)
    assert html.count('name="class" value="M"') == html.count('name="class" value="K"') > 0
    assert _posted_pages(html, "stars_page") == [1, 1, 1, 3, 4, 3, 4]


def test_render_pagination_get_mode_uses_plain_links():
    """GET mode (the Flask pages): plain <a href> links back to the page's
    own path, carrying the other params and the anchor -- no forms."""
    html = render_pagination(
        "/", {"standalone_page": 2}, "sectors_page", 2, 180, anchor="sectors", method="get",
    )
    assert "<form" not in html
    assert 'href="/?standalone_page=2&amp;sectors_page=1#sectors"' in html
    assert 'href="/?standalone_page=2&amp;sectors_page=3#sectors"' in html
    assert 'href="/?standalone_page=2&amp;sectors_page=4#sectors"' in html  # Last
    assert 'aria-current="page">2<' in html


def test_render_pagination_get_mode_escapes_params():
    html = render_pagination("/", {"q": '"><script>'}, "page", 1, 120, method="get")
    assert "<script>" not in html
    assert "q=%22%3E%3Cscript%3E" in html


def test_render_pagination_default_is_still_post():
    html = render_pagination("browse.py", {"db": "x"}, "page", 1, 120)
    assert '<form method="post" action="browse.py' in html
