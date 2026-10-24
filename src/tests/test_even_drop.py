"""MAP.251: a list that cannot hold everything drops evenly across space, the highest-ranked going last."""

from planetgen.db import query


def _star(n, x, y, z, key, population="old"):
    return {"id": n, "x": x, "y": y, "z": z, "luminosity_sol": key / query.BRIGHT_STAR_POPULATION_WEIGHT.get(population, 1.0),
            "population": population}


def _field():
    """A crowded corner (400 stars) and a sparse far side (100 stars) of one band of brightness, plus 12 bright
    stars of the bands above it, scattered over both."""
    crowd = [_star(n, (n % 20) * 0.5, (n // 20) * 0.5, 0.0, 1000.0 + n * 0.1) for n in range(400)]
    far = [_star(1000 + n, 100.0 + (n % 10) * 5.0, (n // 10) * 5.0, 0.0, 1000.0 + n * 0.1) for n in range(100)]
    bright = [_star(2000 + n, 10.0 * n, 3.0 * n, 0.0, 1.0e5 * (n + 1)) for n in range(12)]
    return crowd, far, bright


def test_the_cut_band_loses_the_same_share_of_its_stars_in_every_region():
    crowd, far, bright = _field()
    listed = query._top_by_key(crowd + far + bright, 112)
    assert len(listed) == 112 and {s["id"] for s in bright} <= {s["id"] for s in listed}
    kept_crowd = sum(1 for s in listed if s["id"] < 1000)
    kept_far = sum(1 for s in listed if 1000 <= s["id"] < 2000)
    # 100 slots over 500 stars of one band: about a fifth of each region, not the brightest 100 wherever they lie.
    assert 55 <= kept_crowd <= 105 and 8 <= kept_far <= 35, (kept_crowd, kept_far)
    plain = query._top_by_key(crowd + far + bright, 112, even=False)
    assert sum(1 for s in plain if 1000 <= s["id"] < 2000) < 5, "the plain cut empties the sparse side"


def test_no_brighter_band_is_dropped_while_a_dimmer_one_is_listed():
    crowd, far, bright = _field()
    everything = crowd + far + bright
    for limit in (5, 12, 40, 112, 300):
        listed = query._top_by_key(everything, limit)
        cut = min(query._decade_band(query._map_key(s)) for s in listed)
        assert {s["id"] for s in everything if query._decade_band(query._map_key(s)) > cut} <= {s["id"] for s in listed}


def test_the_lists_are_nested_a_longer_one_keeps_every_star_of_a_shorter_one():
    crowd, far, bright = _field()
    everything = crowd + far + bright
    shorter = {s["id"] for s in query._top_by_key(everything, 60)}
    longer = {s["id"] for s in query._top_by_key(everything, 140)}
    assert shorter <= longer


def test_a_list_that_fits_is_whole_and_in_rank_order():
    crowd, far, bright = _field()
    everything = crowd + far + bright
    listed = query._top_by_key(everything, 600)
    assert [s["id"] for s in listed] == [s["id"] for s in query._top_by_key(everything, 600, even=False)]
    assert len(listed) == len(everything)


def test_the_halo_and_the_satellite_galaxies_go_first():
    crowd, far, bright = _field()
    halo = [_star(3000 + n, 50.0 * (n % 5), 50.0 * (n // 5), 0.0, 5.0e6 + n, "halo") for n in range(25)]
    satellite = [_star(4000 + n, 50.0 * (n % 5), 50.0 * (n // 5), 10.0, 9.0e6 + n, "satellite") for n in range(25)]
    everything = crowd + far + bright + halo + satellite
    whole = len(crowd + far + bright)
    listed = query._top_by_key(everything, whole)
    assert {s["id"] for s in listed} == {s["id"] for s in crowd + far + bright}, "the disc and bulge fill it, brighter halo or not"
    listed = query._top_by_key(everything, whole + 10)
    assert sum(1 for s in listed if s["population"] == "halo") == 10
    assert not any(s["population"] == "satellite" for s in listed)
    tiers = [query._map_tier(s) for s in listed]
    assert tiers == sorted(tiers)


def test_the_choice_is_the_same_whatever_order_the_stars_arrive_in():
    crowd, far, bright = _field()
    everything = crowd + far + bright
    first = query._top_by_key(everything, 77)
    second = query._top_by_key(list(reversed(everything)), 77)
    assert [s["id"] for s in first] == [s["id"] for s in second]


def test_points_are_cut_the_same_way():
    points = [{"id": f"s{n}", "type": "neutron_star", "size": 0.35 + (n % 3) * 0.01, "x": float(n), "y": 0.0, "z": 0.0}
              for n in range(300)] + [{"id": "big", "type": "black_hole", "size": 1.0, "x": 0.0, "y": 0.0, "z": 0.0}]
    listed = query._even_pick(points, 40, lambda p: (-p["size"], p["id"]), lambda p: int(p["size"] / query.SCATTERED_SIZE_BAND),
                              lambda p: p["id"])
    assert len(listed) == 40 and listed[0]["id"] == "big"
    assert max(int(p["id"][1:]) for p in listed if p["id"] != "big") > 200, "the cut band reaches the whole range"


def test_the_single_brightest_object_is_kept_even_by_a_list_of_one():
    crowd, far, bright = _field()
    everything = crowd + far + bright
    assert [s["id"] for s in query._top_by_key(everything, 1)] == [max(bright, key=lambda s: s["luminosity_sol"])["id"]]
    assert query._top_by_key(everything, 20)[0]["id"] == query._top_by_key(everything, 20, even=False)[0]["id"]
