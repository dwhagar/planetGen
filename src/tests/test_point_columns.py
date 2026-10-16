# tests/test_point_columns.py

"""
GEN.335: the stored point of an object is defined once (`db/point_columns.py`). Every table that holds an object is on
exactly one list: a layout with all four parts (a point, a mass, a movement vector and an epoch), or the objects that move
with another's point (a star with its system), or the few that hold no object. A new kind of object cannot ship without
a point, a mass, a vector and an epoch.
"""

import math

import pytest

from planetgen.db import models, point_columns
from planetgen.physics import constants

LAYOUTS, MOVES_WITH, NOT_BODIES = point_columns.LAYOUTS, point_columns.MOVES_WITH, point_columns.NOT_BODIES
TABLES = models.metadata.tables


def _columns(table):
    return [column.name for column in TABLES[table].columns]


def _holds_an_object(table):
    names = _columns(table)
    return "uid" in names or point_columns.holds_a_point(names) or "velocity_x_kms" in names


def test_every_table_holding_an_object_is_on_exactly_one_list():
    lists = (set(LAYOUTS), set(MOVES_WITH), set(NOT_BODIES))
    assert sum(len(found) for found in lists) == len(set().union(*lists)), "a table is on two lists"
    unlisted = sorted(name for name in TABLES if _holds_an_object(name) and not any(name in found for found in lists))
    assert not unlisted, (f"{unlisted} hold an object or a point and are on none of LAYOUTS, MOVES_WITH or NOT_BODIES "
                          f"of db/point_columns.py: give each a point, a mass, a movement vector and an epoch")
    missing = sorted(name for found in lists for name in found if name not in TABLES)
    assert not missing, f"{missing} are listed but are not tables"


@pytest.mark.parametrize("table", sorted(LAYOUTS))
def test_a_layout_names_columns_the_table_has_and_covers_every_part(table):
    layout = LAYOUTS[table]
    assert set(layout.columns()) <= set(_columns(table))
    parts = point_columns.detected_parts(_columns(table))
    if layout.mass_from is not None:
        source, key, column = layout.mass_from
        assert {key, column} <= set(_columns(source))
        parts.add("mass")
    assert parts == set(point_columns.PARTS)


@pytest.mark.parametrize("table", sorted(MOVES_WITH))
def test_an_object_that_moves_with_another_reads_a_layout_and_holds_no_point_of_its_own(table):
    layout, _why = MOVES_WITH[table]
    assert layout in LAYOUTS
    assert point_columns.detected_parts(_columns(table)) & {"point", "vector", "epoch"} == set()


def test_a_table_with_a_point_but_no_mass_no_vector_or_no_epoch_is_seen_to_lack_them():
    assert point_columns.detected_parts(["center_x_pc", "center_y_pc", "center_z_pc"]) == {"point"}
    assert point_columns.detected_parts(["center_x_pc", "center_y_pc"]) == set()
    full = ["position_x_km", "position_y_km", "position_z_km", "mass_kg", "velocity_x_kms", "velocity_y_kms",
            "velocity_z_kms", "epoch_unix"]
    assert point_columns.detected_parts(full) == set(point_columns.PARTS)
    assert point_columns.detected_parts(full[:-1]) == set(point_columns.PARTS) - {"epoch"}
    assert point_columns.detected_parts(["primary_mass_solar"]) == set()


def test_sql_fragments_come_from_the_layout():
    layout = LAYOUTS["bright_stars"]
    assert layout.select("position", "velocity") == (
        "position_x_mpc, position_y_mpc, position_z_mpc, velocity_x_kms, velocity_y_kms, velocity_z_kms")
    assert layout.select("epoch", "due", alias="b") == "b.epoch_unix, b.next_update_due"
    assert layout.assignments("epoch", "due") == "epoch_unix = ?, next_update_due = ?"
    assert layout.values({"epoch_unix": 1.0, "next_update_due": 2.0}, "epoch", "due") == (1.0, 2.0)
    assert LAYOUTS["star_systems"].select("galaxy_position") == "galaxy_x_pc, galaxy_y_pc, galaxy_z_pc"
    with pytest.raises(ValueError):
        layout.select("galaxy_position")


def test_a_galaxy_frame_row_round_trips_through_a_spatial_position():
    layout = LAYOUTS["phenomenon_scatter"]
    row = {"position_x_mpc": 8_000_000.0, "position_y_mpc": -1200.5, "position_z_mpc": 40.25, "velocity_x_kms": 12.5,
           "velocity_y_kms": 225.0, "velocity_z_kms": -3.0, "mass_solar": 8.0, "epoch_unix": 946684800.0,
           "next_update_due": 1.0}
    point = point_columns.spatial_from_row(layout, row)
    assert point.get_coordinates("galactic", "cartesian")[0] == pytest.approx(8_000_000.0 * constants.PARSEC_M / 1000.0)
    assert point.mass_kg == pytest.approx(8.0 * constants.SOLAR_MASS_TO_KG)
    assert point.get_velocity_vector("galactic") == pytest.approx((12500.0, 225000.0, -3000.0))
    back = point_columns.values_from_spatial(layout, point)
    for name in layout.position + layout.velocity + (layout.mass, layout.epoch):
        assert back[name] == pytest.approx(row[name], rel=1e-12)


def test_a_sector_frame_row_with_its_mass_from_the_stars_round_trips():
    layout = LAYOUTS["star_systems"]
    centre = tuple(c * constants.PARSEC_M for c in (8000.0, -500.0, 20.0))
    row = {"position_x_mpc": 1500.0, "position_y_mpc": -250.0, "position_z_mpc": 75.0, "velocity_x_kms": 1.0,
           "velocity_y_kms": 2.0, "velocity_z_kms": 3.0, "epoch_unix": 5.0, "next_update_due": 6.0}
    point = point_columns.spatial_from_row(layout, row, sector_center_m=centre, mass_kg=2.0e30)
    assert point.mass_kg == 2.0e30
    back = point_columns.values_from_spatial(layout, point)
    for name in layout.position + layout.velocity:
        assert back[name] == pytest.approx(row[name], rel=1e-9)
    galaxy = [back[name] for name in layout.galaxy_position]
    assert galaxy == pytest.approx([8000.0 + 1.5, -500.0 - 0.25, 20.0 + 0.075], rel=1e-9)
    assert layout.mass not in back


def test_a_system_frame_row_needs_its_star_and_is_measured_from_it():
    layout = LAYOUTS["planets"]
    row = {"position_x_km": 1.5e8, "position_y_km": -2.0e7, "position_z_km": 3.0e5, "velocity_x_kms": 0.0,
           "velocity_y_kms": 29.8, "velocity_z_kms": 0.1, "mass_kg": 5.97e24, "epoch_unix": 1.0, "next_update_due": 2.0}
    with pytest.raises(ValueError):
        point_columns.spatial_from_row(layout, row)
    star = (2.0e20, -3.0e19, 4.0e17)
    star_velocity = (1.0e5, 2.2e5, -5.0e3)
    point = point_columns.spatial_from_row(layout, row, star_center_m=star, star_velocity_ms=star_velocity)
    assert point.get_coordinates("galactic", "cartesian") == pytest.approx(
        tuple(s + r * 1000.0 for s, r in zip(star, (1.5e8, -2.0e7, 3.0e5))))
    assert point.get_velocity_vector("galactic") == pytest.approx(
        tuple(v + r * 1000.0 for v, r in zip(star_velocity, (0.0, 29.8, 0.1))))
    back = point_columns.values_from_spatial(layout, point)
    for name in layout.position + layout.velocity + ("mass_kg",):
        assert back[name] == pytest.approx(row[name], rel=1e-9, abs=1e-9)
    assert math.isclose(back["epoch_unix"], 1.0)
