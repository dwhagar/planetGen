# tests/test_column_units.py -- DB.76: every numeric column declares the unit it is stored in, in one list
# (`db/column_units.py`), and a unit in a column's name is the unit it is declared in.

import sqlalchemy as sa

from planetgen.db import column_units, models, point_columns
from planetgen.physics import constants, unit_registry

NUMERIC = (sa.Integer, sa.Float, sa.Numeric)


def _numeric_columns():
    return {(table.name, column.name) for table in models.metadata.tables.values() for column in table.columns
            if isinstance(column.type, NUMERIC)}


def test_every_numeric_column_is_declared_and_every_declaration_is_a_column():
    columns = _numeric_columns()
    declared = set(column_units.COLUMNS)
    assert not columns - declared, f"numeric columns with no unit in db/column_units.py: {sorted(columns - declared)}"
    assert not declared - columns, f"declared columns the schema does not have as numbers: {sorted(declared - columns)}"


def test_a_column_is_declared_once_in_a_unit_that_exists():
    seen = {}
    for table, by_unit in column_units.TABLES.items():
        for unit, columns in by_unit.items():
            known = unit in unit_registry.UNITS or unit in column_units.LOCAL_UNITS
            assert known, f"{table}: unknown unit {unit!r}"
            for column in columns:
                assert (table, column) not in seen, f"{table}.{column} is declared in {seen[(table, column)]} and {unit}"
                seen[(table, column)] = unit
    assert len(seen) == len(column_units.COLUMNS)


def test_the_unit_in_a_columns_name_is_the_unit_it_is_declared_in():
    wrong = {key: (unit, column_units.unit_matching_suffix(key[1], key[0])) for key, unit in column_units.COLUMNS.items()
             if column_units.unit_matching_suffix(key[1], key[0]) not in (None, unit)}
    assert not wrong, f"name and declaration disagree (declared, name says): {wrong}"


def test_the_suffix_exceptions_are_numeric_columns_whose_name_does_end_in_a_suffix():
    for (table, column), reason in column_units.SUFFIX_EXCEPTIONS.items():
        assert (table, column) in column_units.COLUMNS and reason
        assert any(column.endswith(suffix) for suffix, _unit in column_units.SUFFIXES), f"{table}.{column} needs none"


def test_a_column_of_a_length_a_mass_or_a_time_is_declared_in_one_of_that_kind():
    for kind, word in (("length", "_km"), ("mass", "_kg"), ("speed", "_kms"), ("angle", "_deg")):
        found = [key for key in column_units.COLUMNS if key[1].endswith(word)]
        assert found and all(column_units.kind_of(*key) == kind for key in found)


def test_a_unit_is_in_the_registry_or_local_and_never_in_both():
    assert not set(column_units.LOCAL_UNITS) & set(unit_registry.UNITS)
    unused = set(column_units.LOCAL_UNITS) - set(column_units.COLUMNS.values())
    assert not unused, f"local units no column uses: {sorted(unused)}"


def test_the_factors_come_from_the_registry_and_the_exact_definitions():
    assert column_units.to_si("planets", "radius_km") == 1000.0
    assert column_units.to_si("star_systems", "galaxy_x_pc") == constants.PARSEC_M
    assert column_units.to_si("star_systems", "position_x_mpc") == constants.PARSEC_M / 1000.0
    assert column_units.to_si("stars", "mass_kg") == 1.0
    assert column_units.to_si("black_holes", "mass_solar") == constants.SOLAR_MASS_TO_KG
    assert column_units.to_si("stars", "luminosity_w") == 1.0
    assert column_units.unit_info("Gyr").factor == 1e9 * column_units.unit_info("yr").factor
    assert column_units.unit_info("km3").factor == 1.0e9


def test_a_unit_without_an_si_form_refuses_a_factor():
    try:
        column_units.to_si("planets", "ocean_ph")
    except ValueError as error:
        assert "no SI form" in str(error)
    else:
        raise AssertionError("pH has no SI form")


def test_a_solar_suffix_is_a_mass_when_the_name_says_mass_and_a_luminosity_otherwise():
    assert column_units.unit_matching_suffix("primary_mass_solar") == "Msun"
    assert column_units.unit_matching_suffix("bright_star_mass_limit_sol") == "Msun"
    assert column_units.unit_matching_suffix("map_luminosity_sol") == "Lsun"
    assert column_units.unit_matching_suffix("relative_density") is None


def test_a_layouts_units_are_the_registrys():
    for layout in point_columns.LAYOUTS.values():
        assert column_units.unit_info(layout.position_unit).kind == "length"
        assert column_units.unit_info(layout.mass_unit).kind == "mass"
    assert point_columns.LAYOUTS["star_systems"].position_unit == "mpc"
    assert point_columns.LAYOUTS["star_systems"].mass_unit == "kg"
    assert point_columns.LAYOUTS["black_holes"].mass_unit == "Msun"
    assert point_columns.LAYOUTS["stars"].position_unit == "au"
