"""The data table's in-memory pieces (UX.41): `in_memory` and `parts_of` in `web/lib/datatable.py`."""

from planetgen.web.lib.datatable import State, in_memory, parts_of

ROWS = [
    {"name": "b", "kind": "x", "n": 2}, {"name": "A", "kind": "y", "n": None},
    {"name": "c", "kind": "x", "n": 1}, {"name": "d", "kind": None, "n": 3},
]
SORTS = {"name": lambda r, desc: r["name"].casefold(), "n": lambda r, desc: r["n"]}
FACETS = {"kind": lambda r: r["kind"]}


def _load(sort="name", descending=False, kinds=(), limit=50, offset=0, facets=True):
    state = State(sort, descending, {"kind": list(kinds)}, 1)
    return in_memory(ROWS, state, limit, offset, facets, SORTS, FACETS, lambda r: [r["name"]])


def test_sorts_ignoring_case_and_both_ways():
    assert [c[0] for c in _load().rows] == ["A", "b", "c", "d"]
    assert [c[0] for c in _load(descending=True).rows] == ["d", "c", "b", "A"]


def test_rows_without_a_value_list_last_in_either_direction():
    assert [c[0] for c in _load("n").rows] == ["c", "b", "d", "A"]
    assert [c[0] for c in _load("n", True).rows] == ["d", "b", "c", "A"]


def test_filters_slice_and_count():
    result = _load(kinds=["x"], limit=1, offset=1)
    assert result.total == 2 and [c[0] for c in result.rows] == ["c"]


def test_facet_counts_ignore_their_own_menu_and_skip_rows_without_a_value():
    result = _load(kinds=["y"])
    assert result.facets["kind"] == [{"value": "x", "label": "x", "count": 2}, {"value": "y", "label": "y", "count": 1}]
    assert _load(facets=False).facets is None


def test_parts_keep_links_and_text_and_drop_other_markup():
    assert parts_of('Nearest: <a href="/system/3">A&amp;B</a> (1.0 ly), <b>x</b>') == [
        "Nearest: ", {"text": "A&B", "href": "/system/3"}, " (1.0 ly), ", "x"]
