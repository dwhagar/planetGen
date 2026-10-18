# tests/test_uncharted_jump_by_value.py

"""
PERF.135: the jump-by-value box of the Uncharted stars and Phenomena tables (PERF.75 gave it to Systems). `from=`
(`start=` in the readers) begins the Uncharted stars at the first one at or after a luminosity, temperature, age, shell
number or star type in the sort's direction, and the Phenomena at a scattered phenomenon's number; the count is of what
is left; a start past the end answers the last page with `past_end`; text that is not a value is a 400 on the API and
the list from the top, with a line saying so, on the web page. Runs against a real throwaway database (the fixtures of
`test_keyset_paging`).
"""

import pytest

from planetgen.db import query as queryDb
from planetgen.physics import constants
from tests.test_api import client  # noqa: F401
from tests.test_keyset_paging import SORTS, _phenomenon, _star, many_stars, phenomena, site, stars  # noqa: F401


def _ids(rows):
    return [row["id"] for row in rows]


def _value(sort, row):
    return {"luminosity": row["luminosity_sol"], "temperature": row["temperature_k"], "age": row["age_gy"],
            "type": row["star_type"], "sector": row["shell_index"]}[sort]


@pytest.mark.parametrize("descending", [False, True])
@pytest.mark.parametrize("sort", ["luminosity", "temperature", "age", "sector"])
def test_a_numeric_start_keeps_the_stars_from_it_on(stars, sort, descending):  # noqa: F811
    everything = queryDb.list_uncharted_systems(stars, sort=sort, descending=descending)
    values = sorted({_value(sort, row) for row in everything})
    chosen = values[len(values) // 2]
    # The box takes the number the column shows: solar luminosities for the luminosity sort.
    text = repr(chosen)
    keep = (lambda value: value <= chosen * (1 + 1e-9)) if descending else (lambda value: value >= chosen * (1 - 1e-9))
    expected = [row for row in everything if keep(_value(sort, row))]
    assert expected and len(expected) < len(everything)
    assert _ids(queryDb.list_uncharted_systems(stars, sort=sort, descending=descending, start=text)) == _ids(expected)
    assert queryDb.count_uncharted_systems(stars, start=text, sort=sort, descending=descending) == len(expected)
    for index in range(0, len(expected), 4):
        by_offset = queryDb.list_uncharted_systems(stars, limit=4, offset=index, sort=sort, descending=descending,
                                                   start=text)
        assert _ids(by_offset) == _ids(expected[index:index + 4])
        if index:
            by_key = queryDb.list_uncharted_systems(stars, limit=4, offset=999, sort=sort, descending=descending,
                                                    start=text, after=expected[index - 1]["id"])
            assert _ids(by_key) == _ids(expected[index:index + 4]), (sort, descending, index)


@pytest.mark.parametrize("descending", [False, True])
def test_a_type_start_keeps_the_types_that_begin_with_it_first(stars, descending):  # noqa: F811
    everything = queryDb.list_uncharted_systems(stars, sort="type", descending=descending)
    start = "B"
    if descending:
        expected = [row for row in everything if row["star_type"] < "C"]
    else:
        expected = [row for row in everything if row["star_type"] >= "B"]
    assert expected and len(expected) < len(everything)
    page = queryDb.list_uncharted_systems(stars, sort="type", descending=descending, start=start)
    assert _ids(page) == _ids(expected)
    assert page[0]["star_type"].startswith("B")
    assert queryDb.count_uncharted_systems(stars, start=start, sort="type", descending=descending) == len(expected)


def test_the_start_is_read_as_the_column_shows_it(stars):  # noqa: F811
    assert queryDb.parse_uncharted_start("luminosity", "2.5") == pytest.approx(2.5 * constants.SOLAR_LUMINOSITY)
    assert queryDb.parse_uncharted_start("temperature", " 6000 ") == 6000.0
    assert queryDb.parse_uncharted_start("age", "0.5") == 0.5
    assert queryDb.parse_uncharted_start("sector", "12") == 12
    assert queryDb.parse_uncharted_start("type", " G2") == "G2"
    assert queryDb.parse_uncharted_start("luminosity", "") is None
    for sort, text in (("luminosity", "bright"), ("temperature", "-3"), ("age", "nan"), ("sector", "1.5"),
                       ("sector", "-1"), ("age", "inf")):
        with pytest.raises(ValueError):
            queryDb.parse_uncharted_start(sort, text)
    with pytest.raises(ValueError):
        queryDb.list_uncharted_systems(stars, sort="age", start="young")
    with pytest.raises(ValueError):
        queryDb.count_uncharted_systems(stars, start="young", sort="age")


@pytest.mark.parametrize("descending", [False, True])
@pytest.mark.parametrize("sort,start", [("luminosity", "3918"), ("sector", "1")])
def test_a_page_from_a_start_is_read_in_index_order_with_no_sort(many_stars, sort, start, descending):  # noqa: F811
    """The start is a range on the sort's index, first page or after a star (MySQL here; the OR form of PERF.75)."""
    condition, params = queryDb._uncharted_start_condition(sort, descending, start)
    columns = queryDb.UNCHARTED_SYSTEM_SORTS[sort].split(",")
    previous = many_stars.execute(
        f"SELECT id, {', '.join(columns)} FROM bright_stars ORDER BY id LIMIT 1 OFFSET 2000").fetchone()
    keyset, keyset_params = queryDb._keyset_condition(columns, [previous[column] for column in columns], descending,
                                                      "id", previous["id"])
    direction = "DESC" if descending else "ASC"
    order = f"ORDER BY {', '.join(f'{c} {direction}' for c in columns)}, id {direction} LIMIT 5"
    for where, args in ((f"star_system_id IS NULL AND {condition}", params),
                        (f"star_system_id IS NULL AND {condition} AND {keyset}", params + keyset_params)):
        plan = many_stars.execute(f"EXPLAIN SELECT id FROM bright_stars WHERE {where} {order}", args).fetchall()
        assert plan[0]["key"] is not None, dict(plan[0])
        assert "filesort" not in (plan[0]["Extra"] or "").lower(), dict(plan[0])


def test_the_last_stars_are_the_end_of_the_sort(stars):  # noqa: F811
    for sort in SORTS:
        for descending in (False, True):
            everything = queryDb.list_uncharted_systems(stars, sort=sort, descending=descending)
            last = queryDb.list_uncharted_systems(stars, limit=5, sort=sort, descending=descending, last=True)
            assert _ids(last) == _ids(everything[-5:])


def test_the_uncharted_route_takes_from_and_answers_the_last_page_past_the_end(client, stars):  # noqa: F811
    top = client.get("/api/uncharted-systems?limit=50&sort=temperature&order=asc").get_json()
    jumped = client.get("/api/uncharted-systems?limit=50&sort=temperature&order=asc&from=12000").get_json()
    assert jumped["total"] == len(jumped["items"]) < top["total"] and not jumped["past_end"]
    assert all(item["temperature_k"] >= 12000 for item in jumped["items"])
    past = client.get("/api/uncharted-systems?limit=5&sort=temperature&order=asc&from=99999").get_json()
    assert past["past_end"] and past["total"] == 5
    assert [item["id"] for item in past["items"]] == [item["id"] for item in top["items"][-5:]]
    assert client.get("/api/uncharted-systems?limit=5&from=").get_json()["total"] == top["total"]
    assert client.get("/api/uncharted-systems?limit=5&sort=age&from=old").status_code == 400
    assert client.get("/api/uncharted-systems?limit=5&sort=age&from=old").get_json()["error"]
    # The stored count is kept for each start (PERF.112): a second read answers the same.
    assert client.get("/api/uncharted-systems?limit=5&sort=temperature&order=asc&from=12000").get_json()["total"] \
        == jumped["total"]


# --- The Phenomena table: the scattered tail from a number -------------------------------------------------------

def test_a_phenomena_start_lists_the_scattered_tail_from_that_number(phenomena):  # noqa: F811
    everything = queryDb.list_phenomena(phenomena, sort="name")
    tail = [row for row in everything if row["scattered"]]
    assert len(tail) == 11
    chosen = tail[4]["scatter_id"]
    expected = [row for row in tail if row["scatter_id"] >= chosen]
    page = queryDb.list_phenomena(phenomena, limit=50, start=str(chosen))
    assert [row["scatter_id"] for row in page] == [row["scatter_id"] for row in expected]
    assert queryDb.count_phenomena(phenomena, start=str(chosen)) == len(expected)
    # By offset and by key the pages count from the start.
    for index in range(0, len(expected), 3):
        by_offset = queryDb.list_phenomena(phenomena, limit=3, offset=index, start=str(chosen))
        assert [row["scatter_id"] for row in by_offset] == [row["scatter_id"] for row in expected[index:index + 3]]
        if index:
            by_key = queryDb.list_phenomena(phenomena, limit=3, offset=999, start=str(chosen),
                                            after=expected[index - 1]["scatter_id"])
            assert [row["scatter_id"] for row in by_key] == [row["scatter_id"] for row in expected[index:index + 3]]
    assert queryDb.count_phenomena(phenomena, start="") == queryDb.count_phenomena(phenomena)
    assert queryDb.count_phenomena(phenomena, types=["nebula"], start=str(chosen)) == 0


def test_a_phenomena_start_must_be_a_whole_number_from_one(phenomena):  # noqa: F811
    assert queryDb.parse_scatter_start(" 7 ") == 7 and queryDb.parse_scatter_start("") is None
    for text in ("seven", "0", "-2", "1.5"):
        with pytest.raises(ValueError):
            queryDb.parse_scatter_start(text)
        with pytest.raises(ValueError):
            queryDb.list_phenomena(phenomena, limit=5, start=text)


def test_the_phenomena_route_takes_from(client, phenomena):  # noqa: F811
    everything = client.get("/api/phenomena?limit=40").get_json()["items"]
    tail = [item for item in everything if item["scattered"]]
    jumped = client.get(f"/api/phenomena?limit=40&from={tail[6]['scatter_id']}").get_json()
    assert [item["scatter_id"] for item in jumped["items"]] == [item["scatter_id"] for item in tail[6:]]
    assert jumped["total"] == len(tail) - 6 and not jumped["past_end"]
    past = client.get("/api/phenomena?limit=4&from=999999").get_json()
    assert past["past_end"] and [item["scatter_id"] for item in past["items"]] == \
        [item["scatter_id"] for item in tail[-4:]] and past["total"] == 4
    assert client.get("/api/phenomena?limit=4&from=abc").status_code == 400
    assert client.get("/api/phenomena?limit=4&from=").get_json()["total"] == len(everything)


# --- The web layer ------------------------------------------------------------------------------------------

def test_the_uncharted_table_passes_the_start_on_and_says_when_it_lay_past_the_end(site):  # noqa: F811
    client, patch, apiclient = site
    seen = []

    def fake(db, limit=None, offset=None, sort=None, descending=False, after=None, start=None):
        seen.append((sort, descending, start))
        if start == "bright":
            raise apiclient.ApiError("bright is not a number", 400)
        return {"items": [_star(41), _star(42)], "total": 2, "past_end": start == "99999"}

    patch.setattr(apiclient, "get_uncharted_systems", fake)
    body = client.get("/table/uncharted_systems?sort=temperature&from=6000&limit=2").get_json()
    assert seen[-1] == ("temperature", False, "6000") and body["note"] == ""
    # Luminosity lists brightest first on the page, which is the API's descending order.
    client.get("/table/uncharted_systems?sort=luminosity&from=2&limit=2")
    assert seen[-1] == ("luminosity", True, "2")
    past = client.get("/table/uncharted_systems?sort=age&from=99999&limit=2").get_json()
    assert "99999" in past["note"] and "last 2" in past["note"]
    bad = client.get("/table/uncharted_systems?sort=age&from=bright&limit=2").get_json()
    assert seen[-1] == ("age", False, "") and "not a value" in bad["note"] and len(bad["rows"]) == 2


def test_the_phenomena_table_passes_the_start_on_and_says_when_it_lay_past_the_end(site):  # noqa: F811
    client, patch, apiclient = site
    seen = []

    def fake(db, limit=None, offset=None, sort=None, descending=False, types=(), descriptors=(), placed=None,
             facets=False, after=None, start=None):
        seen.append(start)
        if start == "abc":
            raise apiclient.ApiError("abc is not a scatter id", 400)
        return {"items": [_phenomenon("Uncharted black hole 9", True, 9)], "total": 1, "past_end": start == "999"}

    patch.setattr(apiclient, "get_phenomena", fake)
    assert client.get("/table/phenomena?from=5&limit=2").get_json()["note"] == ""
    assert seen[-1] == "5"
    assert "last 1" in client.get("/table/phenomena?from=999&limit=2").get_json()["note"]
    assert "not a value" in client.get("/table/phenomena?from=abc&limit=2").get_json()["note"]
    client.get("/table/phenomena?sort=type&from=5&limit=2")
    assert seen[-1] == ""


def test_the_page_draws_the_note_after_the_count(site):  # noqa: F811
    client, patch, apiclient = site
    patch.setattr(apiclient, "get_uncharted_systems",
                  lambda db, limit=None, offset=None, sort=None, descending=False, after=None, start=None:
                  {"items": [_star(41)], "total": 1, "past_end": True})
    html = client.get("/systems?uncharted=1&uncharted_sort=age&uncharted_from=99").get_data(as_text=True)
    assert "Nothing from \u201c99\u201d on; showing the last 1." in html
    assert "from \u201c99\u201d</p>" not in html
    assert 'name="uncharted_from"' in html
