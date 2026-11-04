# tests/test_search_values.py

"""
UX.226: the numeric range filters and sort keys of the search (`db/search_values.py`, `db.query.search`, `GET /api/search`),
checked against the stored columns of a generated sample galaxy (`tests/scored_site.py`): every range returns exactly the bodies
a plain reading of the stored column picks, in the unit asked for, and every sort puts them in order.
"""

import json
import urllib.error
import urllib.parse
import urllib.request

import pytest

from planetgen.db import column_units, query, search_values as sv
from planetgen.physics import score_filters, unit_registry

from tests.scored_site import scored_galaxy  # noqa: F401 -- fixture


def _rows(config, table, columns):
    conn = query.open_readonly(config)
    try:
        return [dict(row) for row in conn.execute(f"SELECT id, name, {', '.join(columns)} FROM {table}").fetchall()]
    finally:
        conn.close()


def _search(config, ranges=None, sorts=None):
    conn = query.open_readonly(config)
    try:
        tags = {name: set() for name in query.SEARCH_TAG_FACETS}
        return query.search(conn, {}, tags, sorts=sorts, ranges=ranges, limit=300)
    finally:
        conn.close()


def _get(base_url, **params):
    url = f"{base_url}/api/search?{urllib.parse.urlencode(params, doseq=True)}"
    try:
        with urllib.request.urlopen(url) as response:
            return response.status, json.loads(response.read())
    except urllib.error.HTTPError as error:
        return error.code, json.loads(error.read())


# --- the registry ------------------------------------------------------------------------------------------------------------

def test_every_key_has_a_unit_that_converts_to_the_unit_it_is_stored_in():
    assert len(sv.KEYS) > 100
    for key, value in sv.KEYS.items():
        for entity in value.entities:
            stored, default = value.stored_unit(entity), value.default_unit(entity)
            assert sv.convert(1.0, default, stored) > 0.0 or stored == "K", (key, entity)
            if value.column is not None:
                assert column_units.COLUMNS[(sv.TABLE_OF[entity], value.column)] == stored, (key, entity)


def test_a_key_never_shadows_a_score_sort_of_the_same_panel():
    for key, value in sv.KEYS.items():
        if "planet" in value.entities or "moon" in value.entities:
            assert key not in score_filters.SORT_COLUMNS, key


def test_the_lists_of_keys_is_the_one_the_api_describes():
    described = {entry["key"]: entry for entry in sv.describe()}
    assert set(described) == set(sv.KEYS)
    assert described["gravity"]["units"]["planet"] == {"unit": "g0", "symbol": "g₀", "kind": "acceleration"}
    assert described["orbit_distance"]["units"]["planet"]["unit"] == "au" and described["orbit_distance"]["units"]["moon"]["unit"] == "km"
    assert described["luminosity"]["entities"] == ["star"] and described["moon_count"]["entities"] == ["planet"]
    assert "gas_CO2" in described and described["gas_CO2"]["units"]["moon"]["unit"] == "bar"


def test_a_value_turns_into_another_unit_of_its_kind_and_not_into_one_of_another():
    assert sv.convert(1.0, "bar", "kPa") == pytest.approx(100.0)
    assert sv.convert(1.0, "au", "km") == pytest.approx(149597870.7)
    assert sv.convert(2.0, "mSv/yr", "mSv/yr") == 2.0
    assert sv.convert(1.0, "Mearth", "kg") == pytest.approx(5.9722e24, rel=1e-3)
    with pytest.raises(ValueError):
        sv.convert(1.0, "bar", "km")


# --- the request ---------------------------------------------------------------------------------------------------------------

def test_the_parameters_give_the_ranges_in_the_unit_given_or_the_keys_own():
    ranges = sv.ranges_from_args({"planet_min_gravity": "0.8", "planet_max_gravity": "1.2", "planet_max_surface_dose": "2",
                                  "moon_min_surface_pressure": "50", "moon_surface_pressure_unit": "kPa",
                                  "star_min_luminosity": "0.5", "planet_max_radius_km": "7000"})
    assert ranges == {"planet": {"gravity": (0.8, 1.2, "g0"), "surface_dose": (None, 2.0, "mSv/yr")},
                      "moon": {"surface_pressure": (50.0, None, "kPa")},
                      "star": {"luminosity": (0.5, None, "Lsun")}}
    assert sv.ranges_from_args({}) == {}


@pytest.mark.parametrize("args,message", [
    ({"planet_min_gravity": "abc"}, "must be a number"),
    ({"planet_min_gravity": "nan"}, "finite"),
    ({"planet_min_gravity": "2", "planet_max_gravity": "1"}, "must not exceed"),
    ({"planet_min_gravity": "1", "planet_gravity_unit": "km"}, "must be a unit of acceleration"),
    ({"planet_min_gravity": "1", "planet_gravity_unit": "furlong"}, "not a unit"),
])
def test_a_bad_bound_or_unit_is_refused_by_name(args, message):
    with pytest.raises(ValueError, match=message):
        sv.ranges_from_args(args)


def test_a_key_of_another_panel_is_not_read():
    assert sv.ranges_from_args({"star_min_gravity": "1", "planet_min_luminosity": "1"}) == {}


# --- the search ----------------------------------------------------------------------------------------------------------------

PLAIN = (("gravity", "gravity_g", "g0"), ("surface_temperature", "surface_temperature_k", "K"),
         ("surface_pressure", "atmospheric_pressure_pa", "Pa"), ("surface_dose", "surface_dose_msv_yr", "mSv/yr"),
         ("ocean_depth", "ocean_depth_km", "km"), ("year_length", "period_years", "yr"))


@pytest.mark.parametrize("table,panel,entity", (("planets", "planets", "planet"), ("moons", "moons", "moon")))
@pytest.mark.parametrize("key,column,stored", PLAIN)
def test_a_range_returns_exactly_the_bodies_whose_stored_value_is_in_it(scored_galaxy, table, panel, entity, key, column, stored):
    bodies = [b for b in _rows(scored_galaxy.config, table, [column]) if b[column] is not None]
    values = sorted(b[column] for b in bodies)
    if len(values) < 4:
        pytest.skip(f"the sample galaxy has too few {table} with {column}")
    least, most = values[len(values) // 4], values[3 * len(values) // 4]
    unit = sv.KEYS[key].default_unit(entity)
    shown = lambda number: sv.convert(number, stored, unit)
    want = {b["name"] for b in bodies if least <= b[column] <= most}
    answer = _search(scored_galaxy.config, ranges={entity: {key: (shown(least), shown(most), unit)}})
    got = answer["results"][panel]["rows"]
    assert {r["name"] for r in got} == want and want
    for row in got:
        assert row["values"][key]["unit"] == unit
        assert shown(least) * (1 - 1e-9) <= row["values"][key]["value"] <= shown(most) * (1 + 1e-9)
    assert answer["active_ranges"][0]["key"] == key and answer["active_ranges"][0]["unit"] == unit


def test_the_bounds_may_be_given_in_another_unit_of_the_same_kind(scored_galaxy):
    bodies = [b for b in _rows(scored_galaxy.config, "planets", ["atmospheric_pressure_pa"]) if b["atmospheric_pressure_pa"]]
    assert bodies
    middle = sorted(b["atmospheric_pressure_pa"] for b in bodies)[len(bodies) // 2]
    in_bar = _search(scored_galaxy.config, ranges={"planet": {"surface_pressure": (middle / 1e5, None, "bar")}})
    in_kpa = _search(scored_galaxy.config, ranges={"planet": {"surface_pressure": (middle / 1e3, None, "kPa")}})
    names = lambda answer: sorted(r["name"] for r in answer["results"]["planets"]["rows"])
    assert names(in_bar) == names(in_kpa) == sorted(b["name"] for b in bodies if b["atmospheric_pressure_pa"] >= middle * (1 - 1e-12))
    first = in_kpa["results"]["planets"]["rows"][0]["values"]["surface_pressure"]
    assert first["unit"] == "kPa" and first["symbol"] == "kPa"


def test_gravity_and_dose_together_give_the_bodies_that_have_both(scored_galaxy):
    bodies = _rows(scored_galaxy.config, "planets", ["gravity_g", "surface_dose_msv_yr"])
    bodies = [b for b in bodies if b["gravity_g"] is not None and b["surface_dose_msv_yr"] is not None]
    assert len(bodies) > 3
    gravities = sorted(b["gravity_g"] for b in bodies)
    doses = sorted(b["surface_dose_msv_yr"] for b in bodies)
    least, most, limit = gravities[1], gravities[-2], doses[len(doses) // 2]
    want = {b["name"] for b in bodies if least <= b["gravity_g"] <= most and b["surface_dose_msv_yr"] <= limit}
    ranges = {"planet": {"gravity": (least, most, "g0"), "surface_dose": (None, limit, "mSv/yr")}}
    got = _search(scored_galaxy.config, ranges=ranges)["results"]["planets"]["rows"]
    assert {r["name"] for r in got} == want and want
    assert all(set(r["values"]) == {"gravity", "surface_dose"} for r in got)


def test_a_body_with_no_value_is_not_in_a_range(scored_galaxy):
    missing = [b for b in _rows(scored_galaxy.config, "planets", ["ocean_depth_km"]) if b["ocean_depth_km"] is None]
    assert missing, "the sample galaxy has no planet without an ocean depth"
    got = _search(scored_galaxy.config, ranges={"planet": {"ocean_depth": (0.0, None, "km")}})["results"]["planets"]["rows"]
    assert not {b["name"] for b in missing} & {r["name"] for r in got}


@pytest.mark.parametrize("table,panel,entity", (("planets", "planets", "planet"), ("moons", "moons", "moon")))
def test_a_sort_puts_the_smallest_first_and_with_a_minus_the_largest_first_and_no_value_last(scored_galaxy, table, panel, entity):
    column = "gravity_g"
    stored = {b["name"]: b[column] for b in _rows(scored_galaxy.config, table, [column])}
    assert len({v for v in stored.values() if v is not None}) > 2
    for sort, reverse in ((("gravity"), False), (("-gravity"), True)):
        # a sort alone does not open a panel: ask for it with a range that keeps every body that has a value
        rows = _search(scored_galaxy.config, ranges={entity: {"gravity": (0.0, None, "g0")}}, sorts={entity: sort})["results"][panel]["rows"]
        order = [stored[r["name"]] for r in rows]
        assert order == sorted(order, reverse=reverse) and len(order) > 2
        assert [r["values"]["gravity"]["value"] for r in rows] == pytest.approx(order)


def test_a_sort_by_a_gas_reads_the_air_and_a_body_without_the_gas_has_none(scored_galaxy):
    conn = query.open_readonly(scored_galaxy.config)
    try:
        rows = conn.execute("SELECT p.name AS name, g.partial_pa AS pa FROM planet_gases g JOIN planets p ON p.id = g.planet_id "
                            "WHERE g.gas = 'N2'").fetchall()
    finally:
        conn.close()
    if not rows:
        pytest.skip("the sample galaxy has no planet with nitrogen")
    middle = sorted(r["pa"] for r in rows)[len(rows) // 2]
    want = {r["name"] for r in rows if r["pa"] >= middle}
    got = _search(scored_galaxy.config, ranges={"planet": {"gas_N2": (middle / 1e5, None, "bar")}})["results"]["planets"]["rows"]
    assert {r["name"] for r in got} == want
    none = _search(scored_galaxy.config, ranges={"planet": {"gas_N2": (None, 0.0, "bar")}})["results"]["planets"]["rows"]
    assert not {r["name"] for r in none} & {r["name"] for r in rows}


def test_the_number_of_moons_of_a_planet_is_counted(scored_galaxy):
    conn = query.open_readonly(scored_galaxy.config)
    try:
        counts = {r["name"]: r["n"] for r in conn.execute(
            "SELECT p.name AS name, COUNT(m.id) AS n FROM planets p LEFT JOIN moons m ON m.planet_id = p.id GROUP BY p.id, p.name").fetchall()}
    finally:
        conn.close()
    assert max(counts.values()) >= 1
    got = _search(scored_galaxy.config, ranges={"planet": {"moon_count": (1.0, None, "count")}})["results"]["planets"]["rows"]
    assert {r["name"] for r in got} == {name for name, n in counts.items() if n >= 1}
    assert all(r["values"]["moon_count"]["value"] >= 1 for r in got)


def test_a_star_is_filtered_and_sorted_by_its_mass_in_solar_masses_and_its_luminosity(scored_galaxy):
    stars = _rows(scored_galaxy.config, "stars", ["mass_kg", "luminosity_w"])
    assert len(stars) > 2
    msun, lsun = unit_registry.UNITS["Msun"].factor, unit_registry.UNITS["Lsun"].factor
    masses = sorted(s["mass_kg"] / msun for s in stars)
    middle = masses[len(masses) // 2]
    want = {s["name"] for s in stars if s["mass_kg"] / msun >= middle * (1 - 1e-12)}
    answer = _search(scored_galaxy.config, ranges={"star": {"mass": (middle * (1 - 1e-12), None, "Msun")}}, sorts={"star": "-mass"})
    rows = answer["results"]["stars"]["rows"]
    assert {r["name"] for r in rows} == want
    shown = [r["values"]["mass"]["value"] for r in rows]
    assert shown == sorted(shown, reverse=True)
    by_name = {s["name"]: s["luminosity_w"] / lsun for s in stars}
    lit = _search(scored_galaxy.config, ranges={"star": {"luminosity": (None, max(by_name.values()), "Lsun")}})["results"]["stars"]["rows"]
    assert {r["name"] for r in lit} == set(by_name)


def test_a_search_with_no_range_answers_the_keys_and_no_active_range(scored_galaxy):
    answer = _search(scored_galaxy.config)
    assert answer["active_ranges"] == [] and {entry["key"] for entry in answer["values"]} == set(sv.KEYS)


# --- the API -------------------------------------------------------------------------------------------------------------------

def test_the_api_takes_a_gravity_and_a_dose_range_and_answers_with_the_values(scored_galaxy):
    bodies = _rows(scored_galaxy.config, "planets", ["gravity_g", "surface_dose_msv_yr"])
    bodies = [b for b in bodies if b["gravity_g"] is not None and b["surface_dose_msv_yr"] is not None]
    gravities = sorted(b["gravity_g"] for b in bodies)
    doses = sorted(b["surface_dose_msv_yr"] for b in bodies)
    least, most, limit = gravities[1], gravities[-2], doses[len(doses) // 2]
    status, body = _get(scored_galaxy.base_url, planet_min_gravity=least, planet_max_gravity=most, planet_max_surface_dose=limit)
    assert status == 200
    want = {b["name"] for b in bodies if least <= b["gravity_g"] <= most and b["surface_dose_msv_yr"] <= limit}
    rows = body["results"]["planets"]["rows"]
    assert {r["name"] for r in rows} == want and want
    assert rows[0]["values"]["gravity"]["symbol"] == "g₀" and rows[0]["values"]["surface_dose"]["unit"] == "mSv/yr"
    assert [entry["key"] for entry in body["active_ranges"]] == ["gravity", "surface_dose"]


def test_the_api_sorts_by_a_key_and_takes_another_unit(scored_galaxy):
    status, body = _get(scored_galaxy.base_url, planet_sort="-surface_pressure", planet_min_surface_pressure=0, planet_surface_pressure_unit="kPa")
    assert status == 200
    shown = [r["values"]["surface_pressure"]["value"] for r in body["results"]["planets"]["rows"]]
    assert shown == sorted(shown, reverse=True) and body["results"]["planets"]["rows"][0]["values"]["surface_pressure"]["unit"] == "kPa"


@pytest.mark.parametrize("params", [
    {"planet_min_gravity": "abc"}, {"planet_min_gravity": "2", "planet_max_gravity": "1"},
    {"planet_min_gravity": "1", "planet_gravity_unit": "km"}, {"moon_max_surface_dose": "inf"},
])
def test_the_api_refuses_a_bad_range_with_a_400(scored_galaxy, params):
    status, body = _get(scored_galaxy.base_url, **params)
    assert status == 400 and "error" in body
