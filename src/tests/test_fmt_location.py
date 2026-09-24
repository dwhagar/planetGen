"""
html/lib/fmt.py "Location:" rendering tests: `nearest_neighbors_location`
(links built from live `system_detail` neighbor rows) and
`linkify_location` (its fallback, parsing names out of the stored string).
No database needed.

Run with: pytest src/tests/test_fmt_location.py
"""
import os
import sys

_SRC_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_SRC_DIR, "html", "lib"))

from fmt import linkify_location, nearest_neighbors_location  # noqa: E402


def test_every_live_neighbor_is_linked_even_when_the_stored_name_is_stale():
    # The stored string names "Alpha Prime", but that system was saved as
    # "Alpha Prime II" (name uniqueness) -- the live neighbor list still
    # links it, under its real name.
    location = "Voranthis Kelmoor -- nearest: Alpha Prime (4.2 ly), Beta (5.1 ly)"
    neighbors = [
        {"id": 7, "name": "Alpha Prime II", "distance_ly": 4.21},
        {"id": 9, "name": "Beta", "distance_ly": 5.08},
    ]
    html = nearest_neighbors_location("db", location, neighbors)
    assert html.startswith("Voranthis Kelmoor -- nearest: ")
    assert html.count('action="system.py"') == 2
    assert "Alpha Prime II" in html and "(4.2 ly)" in html
    assert "Beta" in html and "(5.1 ly)" in html
    assert 'name="id" value="7"' in html and 'name="id" value="9"' in html


def test_stale_stored_name_is_left_unlinked_by_the_string_fallback():
    location = "Voranthis Kelmoor -- nearest: Alpha Prime (4.2 ly)"
    html = linkify_location("db", location, {"Alpha Prime II": 7})
    assert "system.py" not in html
