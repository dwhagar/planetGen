"""
MAP.159: the packed Galaxy Map tile (`planetgen/web/lib/tilepack.py`).
static/tilepack.js reads the same bytes (tests/js/tilepack.test.mjs, fed
by `test_js_unit.tile_pack`).
"""

import json
import math

import pytest

from planetgen.web.lib import tilecache, tilepack


def star(star_id, x, y, z, luminosity, **rest):
    return {"id": star_id, "x": x, "y": y, "z": z, "luminosity_sol": luminosity, **rest}


def test_a_tile_comes_back_with_its_stars_in_order_and_close_to_what_went_in():
    tile = {"clouds": [{"id": 1, "type": "nebula"}], "stars": [
        star(5, 100.001, -2000.5, 7.25, 1.5e5, temperature_k=21500, radius_sol=14.0, star_type="B", system_id=9),
        star(2 ** 40 + 1, 100.5, -2001.0, 8.0, 0.0004, star_type="M"),
    ], "generated": [star(6, 1.0, 1.0, 1.0, 1.0, temperature_k=5778, radius_sol=1.0, star_type="G")]}
    back = tilepack.unpack_tile(tilepack.pack_tile(tile))
    assert back["clouds"] == tile["clouds"]
    assert [s["id"] for s in back["stars"]] == [5, 2 ** 40 + 1]
    first, second = back["stars"]
    assert first["x"] == pytest.approx(100.001, abs=0.5 / 65535 * 2 + 1e-9)
    assert first["y"] == pytest.approx(-2000.5, abs=0.5 / 65535 * 2 + 1e-9)
    assert first["luminosity_sol"] == pytest.approx(1.5e5, rel=1e-3)
    assert first["temperature_k"] == pytest.approx(21500, rel=0.02)
    assert first["radius_sol"] == pytest.approx(14.0, rel=0.04)
    assert first["star_type"] == "B" and first["charted"] is True
    assert "temperature_k" not in second and "radius_sol" not in second and second["charted"] is False
    assert back["generated"][0]["star_type"] == "G"


def test_a_star_costs_about_twenty_bytes_where_json_costs_over_a_hundred():
    stars = [star(i, i * 0.37, -i * 0.11, i * 0.02, 1.0 + i, temperature_k=5000 + i, radius_sol=1.0 + i / 100,
                  star_type="G", system_id=i) for i in range(1000)]
    packed = tilepack.pack_tile({"stars": stars})
    as_json = json.dumps(tilecache.trim_tile({"stars": stars}), separators=(",", ":"))
    assert len(packed) < 21 * len(stars) + 200
    assert len(packed) * 5 < len(as_json)


def test_an_empty_tile_and_one_with_no_stars_are_whole():
    for tile in ({}, {"clouds": [], "stars": [], "generated": []}):
        packed = tilepack.pack_tile(tile)
        assert tilepack.is_whole(packed)
        back = tilepack.unpack_tile(packed)
        assert back["stars"] == [] and back["generated"] == []


def test_a_star_alone_in_its_tile_keeps_its_place():
    back = tilepack.unpack_tile(tilepack.pack_tile({"stars": [star(1, 12.5, -3.25, 99.0, 2.0)]}))
    assert (back["stars"][0]["x"], back["stars"][0]["y"], back["stars"][0]["z"]) == (12.5, -3.25, 99.0)


def test_luminosity_codes_hold_the_dimmest_and_brightest_stars():
    luminosities = [1e-6, 3e-5, 1.0, 9.0e5, 2.5e7]
    back = tilepack.unpack_tile(tilepack.pack_tile(
        {"stars": [star(i, float(i), 0.0, 0.0, value) for i, value in enumerate(luminosities)]}))
    for want, got in zip(luminosities, back["stars"]):
        assert got["luminosity_sol"] == pytest.approx(want, rel=1e-3)


def test_a_dark_or_missing_luminosity_is_zero_and_never_a_crash():
    back = tilepack.unpack_tile(tilepack.pack_tile({"stars": [star(1, 0.0, 0.0, 0.0, 0.0), star(2, 1.0, 1.0, 1.0, -3.0)]}))
    assert [s["luminosity_sol"] for s in back["stars"]] == [0.0, 0.0]


def test_the_incomplete_flag_is_read_from_the_bytes():
    assert tilepack.is_incomplete(tilepack.pack_tile({"stars": [], "incomplete": ["stars"]}))
    assert not tilepack.is_incomplete(tilepack.pack_tile({"stars": [star(1, 0.0, 0.0, 0.0, 1.0)]}))


def test_torn_foreign_and_overlong_bytes_are_not_whole():
    packed = tilepack.pack_tile({"stars": [star(i, float(i), 0.0, 0.0, 1.0) for i in range(5)]})
    assert tilepack.is_whole(packed)
    assert not tilepack.is_whole(packed[:-1]) and not tilepack.is_whole(packed + b"x")
    assert not tilepack.is_whole(b'{"stars":[]}') and not tilepack.is_whole(b"")
    for bad in (packed[:-1], b"PGT1", b'{"stars":[]}'):
        with pytest.raises(tilepack.TilePackError):
            tilepack.unpack_tile(bad)


def test_a_response_carries_its_header_and_each_tile_in_order():
    a = tilepack.pack_tile({"stars": [star(1, 0.0, 0.0, 0.0, 1.0)]})
    b = tilepack.pack_tile({})
    body = tilepack.pack_response({"stamp": "s", "edge_pc": 4.0}, [("1/0/0/0", a), ("1/1/0/0", b)])
    header, tiles = tilepack.unpack_response(body)
    assert header == {"stamp": "s", "edge_pc": 4.0} and list(tiles) == ["1/0/0/0", "1/1/0/0"]
    assert len(tiles["1/0/0/0"]["stars"]) == 1
    assert body.endswith(a + b), "the tiles follow the header one after the other"


def test_the_ranges_cover_what_the_generator_makes():
    low, high = tilepack.LUMINOSITY_LOG_RANGE
    assert low <= math.log10(1e-5) and math.log10(2e6) <= high
    low, high = tilepack.TEMPERATURE_LOG_RANGE
    assert low <= math.log10(1500) and math.log10(300000) <= high
    low, high = tilepack.RADIUS_LOG_RANGE
    assert low <= math.log10(0.005) and math.log10(5000) <= high


@pytest.mark.parametrize("ids", [
    [1, 2 ** 40 + 5, 2 ** 63],
    ["FE81000A2B-0000005-000", "FE81000A2B-0000006-000"],
    [1, "mixed", None],
    ["short", "longer one"],
], ids=["numbers", "printed object ids", "mixed", "texts of different lengths"])
def test_an_id_comes_back_as_it_went_in_whatever_it_is(ids):
    stars = [star(i, float(n), 0.0, 0.0, 1.0) for n, i in enumerate(ids)]
    packed = tilepack.pack_tile({"stars": stars, "generated": stars})
    assert tilepack.is_whole(packed) and not tilepack.is_whole(packed[:-1])
    back = tilepack.unpack_tile(packed)
    assert [s["id"] for s in back["stars"]] == ids and [s["id"] for s in back["generated"]] == ids
