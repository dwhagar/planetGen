# tests/test_precision.py -- DB.77 (UN10): each length column declares its reference frame and how far it must hold.

import math

import pytest

from planetgen.db import column_units, point_columns, precision
from planetgen.physics import constants

LENGTH_COLUMNS = {key for key, unit in column_units.COLUMNS.items() if column_units.unit_info(unit).kind == "length"}


def test_every_length_column_declares_its_role_and_nothing_else_is_declared():
    assert set(precision.PRECISION) == LENGTH_COLUMNS


def test_a_positions_frame_is_the_frame_of_its_tables_point():
    for (table, column), declared in precision.PRECISION.items():
        layout = point_columns.LAYOUTS.get(table)
        if declared.role == "position" and layout is not None and column in layout.position:
            assert declared.frame == layout.frame, (table, column)


def test_every_table_that_stores_a_point_in_a_length_column_declares_a_reach_for_it():
    for table, layout in point_columns.LAYOUTS.items():
        for column in layout.position:
            declared = precision.PRECISION[table, column]
            assert declared.role == "position" and precision.reach_m(table, column) > 0.0, (table, column)


@pytest.mark.parametrize("value, unit, spacing", [
    (constants.AU_M, "m", 3.1e-5), (1e4 * constants.AU_M, "m", 0.25), (constants.PARSEC_M, "m", 4.0),
    (8000.0 * constants.PARSEC_M, "m", 32768.0), (100.0 * constants.AU_M, "m", 2.0e-3),
])
def test_the_spacing_of_a_double_is_the_designs_table(value, unit, spacing):
    assert precision.spacing_m(value, unit) == pytest.approx(spacing, rel=0.05)


@pytest.mark.parametrize("unit", ["m", "km", "au", "pc", "mpc", "ly"])
def test_a_unit_keeps_a_metre_up_to_its_reach_and_not_beyond(unit):
    factor = column_units.unit_info(unit).factor
    reach = precision.metre_reach_m(unit)
    assert precision.spacing_m(reach * 0.999 / factor, unit) <= precision.RESOLUTION_M
    assert precision.spacing_m(reach * 1.001 / factor, unit) > precision.RESOLUTION_M


def test_a_finer_unit_keeps_a_metre_farther_than_a_coarser_one_of_the_same_factor_family():
    assert precision.metre_reach_m("au") < precision.metre_reach_m("km") < precision.metre_reach_m("m")
    assert precision.metre_reach_m("pc") < precision.metre_reach_m("km")


def test_planets_and_moons_place_the_projects_ten_metre_floor_body():
    for table, required in precision.REQUIRED_FLOOR_M.items():
        assert precision.floor_body_m(table) <= required, table
    assert precision.floor_body_m("planets") <= precision.FLOOR_BODY_M


EXPECTED_FLOOR_M = {
    "asteroid_belts": 2.4, "comets": 620.0, "moons": 2.4, "planets": 2.4, "stars": 2.7, "scatter_members": 2.7,
    "orbit_nodes": 2.4, "scatter_nodes": 2.4,
    "bright_stars": 9.2e6, "phenomenon_scatter": 9.2e6, "star_systems": 9.0e6, "rogue_planets": 9.0e6,
    "interstellar_comets": 9.0e6, "asteroid_fields": 9.0e6, "black_holes": 9.0e6, "neutron_stars": 9.0e6,
    "nebulae": 9.0e6, "supernova_remnants": 9.0e6, "quasars": 9.0e6, "facilities": 9.0e6,
    "nebula_block_stats": 9.0e6, "sector_path_knots": 9.0e6, "orbital_event_bodies": 9.0e6, "orbital_events": 2.4,
}
"""dict: The smallest body each table places well, metres, to two figures: a change of a column's unit or a reach moves
one of these, and the table in docs/design/unified-unit-quantity.md (section 6.1) is updated with it."""


def test_the_floor_of_every_table_with_a_position_is_the_documented_one():
    tables = {table for (table, _c), declared in precision.PRECISION.items() if declared.reach is not None}
    assert tables == set(EXPECTED_FLOOR_M)
    for table, expected in EXPECTED_FLOOR_M.items():
        assert precision.floor_body_m(table) == pytest.approx(expected, rel=0.05), table


def test_the_stability_ceiling_is_the_hill_sphere_of_the_largest_star_at_the_suns_distance():
    from planetgen.galaxy import enclosed_mass
    from planetgen.generation.star import calculate_hill_sphere
    parsecs = 8000.0
    hill_m = calculate_hill_sphere(parsecs * constants.PARSEC_M, 150.0 * constants.SOLAR_MASS_TO_KG,
                                   enclosed_mass.current_enclosed_mass_kg(parsecs))
    assert precision.REACHES["stability-ceiling"].metres == pytest.approx(hill_m, rel=0.05)


def test_the_galaxy_reach_is_the_largest_edge_the_dials_allow():
    assert precision.REACHES["galaxy-edge"].metres == pytest.approx(262140.0 * constants.PARSEC_M)
    assert math.isfinite(precision.REACHES["sector"].metres)


def test_a_position_or_separation_declares_a_reach_and_no_other_role_does():
    with pytest.raises(ValueError):
        precision.Precision("position")
    with pytest.raises(ValueError):
        precision.Precision("size", "system", "bound-orbit")
    with pytest.raises(ValueError):
        precision.Precision("radius")
