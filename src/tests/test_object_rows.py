# tests/test_object_rows.py -- DB.74: a row is written from a list of fields (`db/object_rows.py`), one list per kind table.

import types

import pytest

from planetgen.db import models, object_rows, store
from planetgen import tuning
from planetgen.physics import atmosphere, constants, life_classes

KIND_FIELDS = {
    "stars": store.STAR_FIELDS, "planets": store.PLANET_FIELDS, "moons": store.MOON_FIELDS, "asteroid_belts": store.BELT_FIELDS, "comets": store.COMET_FIELDS,
    "black_holes": store.BLACK_HOLE_FIELDS, "neutron_stars": store.NEUTRON_STAR_FIELDS,
    "rogue_planets": store.ROGUE_PLANET_FIELDS, "interstellar_comets": store.INTERSTELLAR_COMET_FIELDS,
    "asteroid_fields": store.ASTEROID_FIELD_FIELDS,
}


class _Connection:
    def __init__(self):
        self.calls = []

    def execute(self, sql, params=()):
        self.calls.append((sql, params))
        return types.SimpleNamespace(lastrowid=len(self.calls))


@pytest.mark.parametrize("table", sorted(KIND_FIELDS))
def test_a_kinds_fields_name_columns_of_its_table_once_each(table):
    columns = object_rows.columns_of(KIND_FIELDS[table])
    assert len(columns) == len(set(columns))
    assert set(columns) <= {column.name for column in models.metadata.tables[table].columns}


@pytest.mark.parametrize("table", sorted(KIND_FIELDS))
def test_a_kinds_fields_fill_every_column_the_table_requires(table):
    required = {column.name for column in models.metadata.tables[table].columns
                if not column.nullable and column.default is None and column.server_default is None
                and not column.primary_key and not column.computed}
    assert required <= set(object_rows.columns_of(KIND_FIELDS[table])), table


def test_the_constructors_read_the_object_and_the_context():
    obj = types.SimpleNamespace(mass=5.0, size=3.0)
    fields = (object_rows.attr("mass_kg", "mass"), object_rows.optional("missing"), object_rows.context("sector_id"),
              object_rows.computed(("a", "b"), lambda o, c: (o.size, c["sector_id"] + 1)))
    assert object_rows.columns_of(fields) == ("mass_kg", "missing", "sector_id", "a", "b")
    assert object_rows.values_of(fields, obj, {"sector_id": 7}) == [5.0, None, 7, 3.0, 8]


def test_a_field_returning_the_wrong_number_of_values_is_refused():
    with pytest.raises(ValueError):
        object_rows.values_of((object_rows.computed(("a", "b"), lambda o, c: (1,)),), None, {})


def test_insert_row_writes_one_plain_insert_and_returns_the_id():
    connection = _Connection()
    fields = (object_rows.attr("name"), object_rows.attr("age_gy", "age"))
    row_id = object_rows.insert_row(connection, "stars", fields, types.SimpleNamespace(name="A", age=2.0))
    assert row_id == 1
    assert connection.calls == [("INSERT INTO stars (name, age_gy) VALUES (?, ?)", ("A", 2.0))]


def test_a_belt_is_written_in_kilometres_with_a_stated_zero_point_and_its_own_mass():
    belt = types.SimpleNamespace(distance=2.0, lower_limit=1.5, upper_limit=3.0, density=0.4, mass_kg=1.5e21,
                                 chondrite_group="CM", get_composition_summary=lambda: "rock")
    connection = _Connection()
    context = {"star_system_id": 5, "star_id": 2, "orbital_index": 1}
    object_rows.insert_row(connection, "asteroid_belts", store.BELT_FIELDS, belt, context)
    columns = object_rows.columns_of(store.BELT_FIELDS)
    row = dict(zip(columns, connection.calls[0][1]))
    assert row["distance_km"] == 2.0 * constants.AU_TO_KM and row["upper_limit_km"] == 3.0 * constants.AU_TO_KM
    assert row["composition_summary"] == "rock" and row["star_id"] == 2 and row["orbital_index"] == 1
    assert [row[name] for name in columns if name.startswith(("position_", "velocity_"))] == [0.0] * 6
    assert row["mass_kg"] == 1.5e21 and row["chondrite_group"] == "CM"


def test_an_anchored_remnant_leaves_its_motion_and_spin_to_its_star():
    hole = types.SimpleNamespace(name="H", mass_class="stellar", mass_solar=10.0, event_horizon_radius_km=30.0,
                                 spin=0.5, has_accretion_disk=False, temperature=0.0, luminosity=0.0, age=1.0,
                                 galactic_orbital_period_gy=0.2, galactic_orbital_phase_deg=10.0,
                                 galactic_min_update_interval_years=5.0, galactic_orbital_speed_kms=200.0, spin_axis_x=0.0, spin_axis_y=0.0,
                                 spin_axis_z=1.0, axial_tilt_deg=0.0)
    columns = object_rows.columns_of(store.BLACK_HOLE_FIELDS)
    anchored = dict(zip(columns, object_rows.values_of(
        store.BLACK_HOLE_FIELDS, hole, {"star_id": 3, "sector_id": None, "placement": None})))
    assert anchored["star_id"] == 3 and anchored["velocity_x_kms"] is None
    assert anchored["galactic_orbital_period_gy"] is None and anchored["spin_axis_z"] is None
    placement = {"center_x_pc": 1000.0, "center_y_pc": 0.0, "center_z_pc": 5.0, "galactic_radius_pc": 1000.0}
    free = dict(zip(columns, object_rows.values_of(
        store.BLACK_HOLE_FIELDS, hole, {"star_id": None, "sector_id": 4, "placement": placement})))
    assert free["galactic_orbital_period_gy"] == 0.2 and free["spin_axis_z"] == 1.0 and free["center_x_pc"] == 1000.0
    assert free["velocity_x_kms"] is not None and free["galactic_radius_pc"] == 1000.0


# -- the load (DB.74 S3): the same lists read a stored row back ----------------------------------------------------

import pytest as _pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw


_STAR_CONTEXT = {"star_system_id": 1, "role": "single", "parent_node_id": None, "slot": None, "offset_au": (0.0, 0.0, 0.0), "flags": 0}


def _stored(fields, obj, context):
    """The row an insert of `obj` writes, as the mapping a `SELECT` hands back."""
    return dict(zip(object_rows.columns_of(fields), object_rows.values_of(fields, obj, context)))


def _same(stored, loaded):
    if isinstance(stored, (list, tuple)):
        assert len(stored) == len(loaded)
        for a, b in zip(stored, loaded):
            _same(a, b)
    elif isinstance(stored, float) or isinstance(loaded, float):
        assert loaded == _pytest.approx(stored, rel=1e-12, abs=1e-12)
    else:
        assert stored == loaded


def _system(seed):
    draw.set_run_seed(seed)
    return StarSystem(system_config=SystemConfig())


_STAR_ROW_ONLY = {"state", "decay_age_yr"}
"""set: Columns the star row holds that `Star.to_dict()` does not carry: they are derived from the star by `stellar_state.set_state` (GEN.785, GEN.786), not serialised with it."""


@_pytest.mark.parametrize("seed", [1, 2, 3])
def test_a_star_is_read_back_as_the_dict_it_was_written_from(seed):
    star = _system(seed).primary_star       # a binary pair is saved as its two stars, not as the proxy that stands for both
    row = _stored(store.STAR_FIELDS, star, _STAR_CONTEXT)
    loaded = store._star_row_to_dict(row, 12.5)
    assert loaded["galactic_orbital_speed_kms"] == 12.5
    original = star.to_dict()
    assert set(loaded) - _STAR_ROW_ONLY == set(original), set(loaded) ^ set(original)
    for key, value in loaded.items():
        if key != "galactic_orbital_speed_kms" and key not in _STAR_ROW_ONLY:
            _same(original[key], value)


def test_a_stored_temperature_that_is_whole_comes_back_as_the_int_the_generator_makes():
    assert store._whole_number(5800.0) == 5800 and isinstance(store._whole_number(5800.0), int)
    assert store._whole_number(3676064.7) == 3676064.7
    assert store._whole_number(None) is None


def test_a_null_reflex_offset_and_a_null_wide_binary_limit_come_back_as_the_objects_defaults():
    star = _system(1).star
    row = _stored(store.STAR_FIELDS, star, _STAR_CONTEXT)
    row.update(reflex_offset_x_km=None, reflex_offset_y_km=None, reflex_offset_z_km=None, wide_binary_a_crit_km=None)
    loaded = store._star_row_to_dict(row, 0.0)
    assert (loaded["reflex_offset_x"], loaded["reflex_offset_y"], loaded["reflex_offset_z"]) == (0.0, 0.0, 0.0)
    assert loaded["a_crit_au"] is None


def test_a_stored_metallicity_is_not_lost_on_the_way_back():
    star = _system(1).star
    row = _stored(store.STAR_FIELDS, star, _STAR_CONTEXT)
    assert store._star_row_to_dict(row, 0.0)["metallicity_fe_h"] == star.metallicity_fe_h is not None


class _Rows(list):
    def fetchall(self):
        return self


class _SpectrumConnection:
    def __init__(self, layers=None, gases=None, scores=None, living=None):
        self.living = living or []
        self.layers = layers or []
        self.gases = gases or {}
        self.scores = scores or {}

    def execute(self, sql, params=()):
        if "_interior_layers" in sql:
            return _Rows(self.layers)
        if "_gases" in sql:  # DB.71: the air is read from planet_gases or moon_gases, one row a gas
            column = "moon_id" if "moon_gases" in sql else "planet_id"
            return _Rows([{column: params[0], "gas": gas, "partial_pa": pa} for gas, pa in self.gases.items()])
        if "_life_class_scores" in sql:  # GEN.889: the class scores are read from planet_life_class_scores or moon_life_class_scores
            column = "moon_id" if "moon_life_class_scores" in sql else "planet_id"
            return _Rows([{column: params[0], "life_class": name, "score": score} for name, score in self.scores.items()])
        if "_life " in sql:  # GEN.1209: the classes that live on the body are read from planet_life or moon_life
            column = "moon_id" if "moon_life" in sql else "planet_id"
            return _Rows([{column: params[0], "life_class": entry["life_class"], "share": entry["share"], "life_stage": entry["stage"],
                           "life_milestone_age_gyr": entry["milestone_age_gyr"], "life_system_age_gyr": entry["system_age_gyr"],
                           "life_animal_grade": entry["animal_grade"], "life_complexity": entry["complexity"],
                           **dict(zip(tuning.LIFE_SCHEDULE_FIELDS, entry["ages"])),
                           **dict(zip(("life_microbial", "life_flora", "life_fauna", "life_sentient"), entry["kinds"]))} for entry in self.living])
        return _Rows([{"spectrum_type": "visible", "value": 1.5}, {"spectrum_type": "non_visible", "value": 3.5}])


def test_planets_moons_belts_and_comets_are_read_back_as_the_dicts_they_were_written_from():
    checked, wanted = set(), {"planet", "moon", "belt", "comet"}
    for seed in range(1, 61):                       # the first seeds that hold a belt move with the generator: read until every kind was seen
        _check_bodies(_system(seed), checked)
        if checked == wanted:
            break
    assert checked == wanted


def _check_bodies(system, checked):
    context = {"star_system_id": 1, "star_id": 2, "orbital_index": 0}
    for body in system.planets:
        if hasattr(body, "moons"):
            row = _stored(store.PLANET_FIELDS, body, context)
            row["id"] = 7
            loaded = store._planet_or_moon_row_to_dict(_SpectrumConnection(body.interior_layers, getattr(body, atmosphere.GAS_FIELD, None), getattr(body, life_classes.SCORES_FIELD, None), getattr(body, life_classes.LIVING_FIELD, None)), row, False)
            original = body.to_dict()
            for key in set(loaded) - {"moons", "reflection_spectrum_visible", "reflection_spectrum_non_visible"}:
                _same(original[key], loaded[key])
            assert loaded["reflection_spectrum_visible"] == [1.5] and loaded["moons"] == []
            checked.add("planet")
            for moon in body.moons:
                row = _stored(store.MOON_FIELDS, moon, dict(context, planet_id=7))
                row["id"] = 8
                loaded = store._planet_or_moon_row_to_dict(_SpectrumConnection(moon.interior_layers, getattr(moon, atmosphere.GAS_FIELD, None), getattr(moon, life_classes.SCORES_FIELD, None), getattr(moon, life_classes.LIVING_FIELD, None)), row, True)
                assert "reflex_offset_x" not in loaded and loaded["is_moon"] is True
                original = moon.to_dict()
                for key in set(loaded) - {"moons", "reflection_spectrum_visible", "reflection_spectrum_non_visible",
                                          "is_moon"}:
                    _same(original[key], loaded[key])
                checked.add("moon")
        else:
            row = _stored(store.BELT_FIELDS, body, context)
            loaded = store._belt_row_to_dict(row, [{"component": "ice", "concentration": 0.5}])
            original = body.to_dict()
            for key in ("distance", "lower_limit", "upper_limit", "density", "mass_kg"):
                _same(original[key], loaded[key])
            assert loaded["body_type"] == "a" and loaded["composition"] == [["ice", 0.5]]
            checked.add("belt")
    for comet in system.comets:
        row = _stored(store.COMET_FIELDS, comet, {"star_system_id": 1, "star_id": 2})
        loaded = store._comet_row_to_dict(row, [{"component": "ice"}])
        original = comet.to_dict()
        for key in set(loaded) - {"composition"}:
            _same(original[key], loaded[key])
        assert loaded["is_active"] is bool(comet.is_active)
        checked.add("comet")


def test_the_load_helpers_invert_their_conversions():
    row = {"d_km": 300000.0, "n_km": None, "flag": 1, "a": 1, "b": 2}
    fields = (object_rows.scaled("d_km", "d", 150000.0), object_rows.scaled("n_km", "n", 150000.0, null=0.0, nullable=True),
              object_rows.each(("a", "b")), object_rows.context("sector_id"))
    assert object_rows.load_dict(fields, row) == {"d": 2.0, "n": 0.0, "a": 1, "b": 2}
    assert object_rows.load_dict(fields, row, only=("a",)) == {"a": 1, "b": 2}
    obj = types.SimpleNamespace(d=2.0, n=None, a=1, b=2)
    assert object_rows.values_of(fields[:3], obj, {}) == [300000.0, None, 1, 2]
