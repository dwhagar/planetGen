# tests/test_systems_jump_by_value.py

"""
PERF.75: the Systems table's jump-by-value box. `from=` (`start=` in the readers) begins the list at the first name, or
the first sector, from a typed text on, in the sort's direction, the count is of what is left, and the pages by offset
and by key count from there. The web layer reads it into `State.start`, keeps it only under the sorts that have a box,
and draws the box. Runs against a real throwaway database (the fixture of `test_slow_read_paths`).
"""

import pytest

from planetgen.db import query as queryDb
from planetgen.web.lib import datatable
from planetgen.web.lib.datatable import Column, Table, parse_state, state_params, view
from tests.test_api import client  # noqa: F401
from tests.test_slow_read_paths import galaxy  # noqa: F401


def _names(rows):
    return [row["name"] for row in rows]


def _ids(rows):
    return [row["id"] for row in rows]


def test_the_text_after_a_prefix_steps_the_last_character_up():
    assert queryDb._after_prefix("M") == "N"
    assert queryDb._after_prefix("Alpha S") == "Alpha T"
    assert queryDb._after_prefix("z\U0010ffff") == "{"
    assert queryDb._after_prefix("\U0010ffff") is None
    assert queryDb._after_prefix("") is None
    assert queryDb._after_prefix("\ud7ff") == "\ue000"


@pytest.mark.parametrize("descending", [False, True])
def test_a_name_start_keeps_the_names_from_it_on(galaxy, descending):  # noqa: F811
    everything = queryDb.list_systems(galaxy, sort="name", descending=descending)
    names = _names(everything)
    start = sorted(names, key=str.casefold)[len(names) // 2][:2]
    if descending:
        bound = queryDb._after_prefix(start.casefold())
        expected = [row for row in everything if row["name"].casefold() < bound]
    else:
        expected = [row for row in everything if row["name"].casefold() >= start.casefold()]
    assert expected and len(expected) < len(everything)
    assert _ids(queryDb.list_systems(galaxy, sort="name", descending=descending, start=start)) == _ids(expected)
    assert queryDb.count_systems(galaxy, start=start, sort="name", descending=descending) == len(expected)
    # The pages by offset and by key count from the start.
    for index in range(0, len(expected), 4):
        by_offset = queryDb.list_systems(galaxy, limit=4, offset=index, sort="name", descending=descending, start=start)
        assert _ids(by_offset) == _ids(expected[index:index + 4])
        if index:
            by_key = queryDb.list_systems(galaxy, limit=4, offset=999, sort="name", descending=descending,
                                          start=start, after=expected[index - 1]["id"])
            assert _ids(by_key) == _ids(expected[index:index + 4])


@pytest.mark.parametrize("descending", [False, True])
@pytest.mark.parametrize("start", ["Beta", "B", "Delta", "Zeta", "alpha s"])
def test_a_sector_start_keeps_the_sectors_from_it_on_and_the_systems_with_none(galaxy, descending, start):  # noqa: F811
    everything = queryDb.list_systems(galaxy, sort="sector", descending=descending)
    if descending:
        bound = queryDb._after_prefix(start.casefold())
        expected = [row for row in everything if row["sector_name"] is None or row["sector_name"].casefold() < bound]
    else:
        expected = [row for row in everything if row["sector_name"] is None or row["sector_name"].casefold() >= start.casefold()]
    assert expected  # the standalone systems at least
    assert _ids(queryDb.list_systems(galaxy, sort="sector", descending=descending, start=start)) == _ids(expected)
    assert queryDb.count_systems(galaxy, start=start, sort="sector", descending=descending) == len(expected)
    for index in range(0, len(expected) + 4, 3):
        page = queryDb.list_systems(galaxy, limit=4, offset=index, sort="sector", descending=descending, start=start)
        assert _ids(page) == _ids(expected[index:index + 4]), (start, descending, index)


def test_a_start_is_ignored_under_the_sorts_without_a_box_and_when_empty(galaxy):  # noqa: F811
    for sort in ("octant",):
        assert _ids(queryDb.list_systems(galaxy, sort=sort, start="Z")) == _ids(queryDb.list_systems(galaxy, sort=sort))
        assert queryDb.count_systems(galaxy, start="Z", sort=sort) == queryDb.count_systems(galaxy)
    assert _ids(queryDb.list_systems(galaxy, start="")) == _ids(queryDb.list_systems(galaxy))
    assert _ids(queryDb.list_systems(galaxy, start=None, limit=5, sort="sector")) == \
        _ids(queryDb.list_systems(galaxy, limit=5, sort="sector"))


def test_a_start_works_with_the_filters(galaxy):  # noqa: F811
    everything = queryDb.list_systems(galaxy, sort="sector", binary=True)
    expected = [row for row in everything if row["sector_name"] is None or row["sector_name"] >= "Beta"]
    assert _ids(queryDb.list_systems(galaxy, sort="sector", binary=True, start="Beta", limit=50)) == _ids(expected)
    assert queryDb.count_systems(galaxy, binary=True, start="Beta", sort="sector") == len(expected)


def test_the_systems_route_takes_from(client, galaxy):  # noqa: F811
    top = client.get("/api/systems?limit=50&sort=sector").get_json()
    jumped = client.get("/api/systems?limit=50&sort=sector&from=Beta").get_json()
    assert jumped["total"] == len(jumped["items"]) < top["total"]
    assert all(item["sector_name"] is None or item["sector_name"] >= "Beta" for item in jumped["items"])
    by_name = client.get("/api/systems?limit=50&from=" + _names(top["items"])[10][:1]).get_json()
    assert by_name["total"] == len(by_name["items"])
    assert client.get("/api/systems?limit=5&from=").get_json()["total"] == top["total"]


# --- The web layer ------------------------------------------------------------------------------------------

TABLE = Table("things", "Things", [Column("name", "Name"), Column("sector", "Sector"), Column("kind", "Kind")],
              lambda *args, **kwargs: None, prefix="things_", jump={"name": "A letter", "sector": "A sector"})


class _Args(dict):
    def getlist(self, key):
        value = self.get(key)
        return [] if value is None else [value]


def test_the_start_is_read_only_under_a_sort_with_a_box():
    assert parse_state(TABLE, _Args({"things_from": " Be "})).start == "Be"
    assert parse_state(TABLE, _Args({"things_from": "Be", "things_sort": "sector"})).start == "Be"
    assert parse_state(TABLE, _Args({"things_from": "Be", "things_sort": "kind"})).start == ""
    assert parse_state(TABLE, _Args({"from": "x" * 500}), prefix="").start == "x" * datatable.MAX_VALUE_LENGTH
    assert "things_from" in TABLE.owned_params()
    assert "x_from" not in Table("x", "X", [Column("a", "A")], None, prefix="x_").owned_params()


def test_the_start_rides_in_the_links_only_while_the_sort_stays():
    state = parse_state(TABLE, _Args({"things_from": "Be"}))
    assert ("things_from", "Be") in state_params(TABLE, state)
    assert ("things_from", "Be") in state_params(TABLE, state, page=2)
    assert all(pair[0] != "things_from" for pair in state_params(TABLE, state, sort="sector", descending=False))
    assert all(pair[0] != "things_from" for pair in state_params(TABLE, state, sort="kind"))


def test_the_view_gives_the_box_its_hint_fields_and_a_way_back():
    state = parse_state(TABLE, _Args({"things_from": "Be", "things_sort": "sector", "things_order": "desc"}))
    box = view(TABLE, state, datatable.Result([], 0), "/things", "/table/things")["jump"]
    assert box["param"] == "things_from" and box["start"] == "Be" and box["hint"] == "A sector"
    assert box["hints"] == {"name": "A letter", "sector": "A sector"}
    assert ("things_sort", "sector") in box["fields"] and ("things_order", "desc") in box["fields"]
    assert all(name != "things_from" for name, _ in box["fields"])
    assert "from" not in box["clear_href"]
    plain = view(Table("p", "P", [Column("a", "A")], None), parse_state(Table("p", "P", [Column("a", "A")], None), _Args()),
                 datatable.Result([], 0), "/p", "/table/p")
    assert plain["jump"] is None
