# tests/test_object_rows.py -- DB.74: a row is written from a list of fields (`db/object_rows.py`), one list per kind table.

import types

import pytest

from planetgen.db import models, object_rows, store
from planetgen.physics import constants

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


def test_a_belt_is_written_in_kilometres_with_a_stated_zero_point_and_mass():
    belt = types.SimpleNamespace(distance=2.0, lower_limit=1.5, upper_limit=3.0, density=0.4,
                                 get_composition_summary=lambda: "rock")
    connection = _Connection()
    context = {"star_system_id": 5, "star_id": 2, "orbital_index": 1}
    object_rows.insert_row(connection, "asteroid_belts", store.BELT_FIELDS, belt, context)
    columns = object_rows.columns_of(store.BELT_FIELDS)
    row = dict(zip(columns, connection.calls[0][1]))
    assert row["distance_km"] == 2.0 * constants.AU_TO_KM and row["upper_limit_km"] == 3.0 * constants.AU_TO_KM
    assert row["composition_summary"] == "rock" and row["star_id"] == 2 and row["orbital_index"] == 1
    assert [row[name] for name in columns if name.startswith(("position_", "velocity_"))] == [0.0] * 6
    assert row["mass_kg"] == 0.0


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
