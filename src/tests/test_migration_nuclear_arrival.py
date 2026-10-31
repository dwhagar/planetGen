# tests/test_migration_nuclear_arrival.py

"""
GEN.1174: migration 0197 adds `nuclear_arrival_age_gy` to `stars` and `bright_stars`, the age MIST's tracks put a pre-main-sequence star on the
main sequence, kept beside the earlier age its radius settles. A row from before v197 keeps NULL there; a star the scatter placed writes it and the
star built from the row keeps it.
"""

import pytest

from planetgen.db import store
from planetgen.generation import bright_stars as brightStars
from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.generation.system import StarSystem
from planetgen.physics import pre_main_sequence
from planetgen.physics import stellar_evolution as evolution
from planetgen.util import draw
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

COLUMN = "nuclear_arrival_age_gy"


def _types(conn, table):
    return {row["COLUMN_NAME"]: row["DATA_TYPE"] for row in conn.execute(
        "SELECT COLUMN_NAME, DATA_TYPE FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ?", (table,)).fetchall()}


def test_a_v196_galaxy_gets_the_column_empty(mysql_config):
    load_old_schema(mysql_config, 194)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("stars", "bright_stars"):
            assert COLUMN not in _types(conn, table)
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("stars", "bright_stars"):
            assert _types(conn, table)[COLUMN] == "double", table
    finally:
        conn.close()


def _pre_main_sequence_params(mass_sol=0.2, age_gy=0.05):
    return evolution.star_params(mass_sol, age_gy, evolution.evolve_star(mass_sol, age_gy, fe_h=0.0))


def test_the_stored_params_of_a_pre_main_sequence_star_carry_both_times():
    params = _pre_main_sequence_params()
    assert params["yerkes_class"] == "PMS"
    assert params["phase_end_age_gy"] == pytest.approx(pre_main_sequence.contraction_time_yr(0.2, 0.0) / 1e9)
    assert params["nuclear_arrival_age_gy"] == pytest.approx(pre_main_sequence.nuclear_arrival_yr(0.2, 0.0) / 1e9)
    assert params["nuclear_arrival_age_gy"] > 3.0 * params["phase_end_age_gy"]
    main_sequence = evolution.star_params(1.0, 4.6, evolution.evolve_star(1.0, 4.6, fe_h=0.0))
    assert main_sequence["nuclear_arrival_age_gy"] is None


def test_a_star_built_from_the_params_keeps_the_nuclear_time():
    params = _pre_main_sequence_params()
    star = Star.from_params(params, SystemConfig(), name="Young")
    assert star.nuclear_arrival_age_gy == params["nuclear_arrival_age_gy"]
    assert star.phase_end_age_gy == params["phase_end_age_gy"]
    assert Star.from_params({**params, "nuclear_arrival_age_gy": None}, SystemConfig(), name="Old").nuclear_arrival_age_gy is None


def test_a_scattered_row_writes_the_column_and_the_star_built_from_it_keeps_it():
    params = _pre_main_sequence_params()
    row = brightStars._row(1, 2, 3, (0.0, 0.0, 0.0), "disk", params, __import__("random").Random(4))
    assert row[store.BRIGHT_STAR_COLUMNS.index(COLUMN)] == params["nuclear_arrival_age_gy"]
    stored = dict(zip(store.BRIGHT_STAR_COLUMNS, row))
    assert brightStars.star_params(stored)["nuclear_arrival_age_gy"] == params["nuclear_arrival_age_gy"]


def test_the_column_is_stored_with_the_star_and_comes_back(mysql_config):
    draw.set_run_seed(9)
    config = SystemConfig()
    config.PLANETS = False
    config.BINARY_SYSTEM = False
    params = _pre_main_sequence_params()
    system = StarSystem(system_config=config, primary_star_params=params)
    conn = store.get_connection(mysql_config)
    try:
        system_id = store.insert_star_system(conn, system, config)
        conn.commit()
        row = conn.execute("SELECT yerkes_class, nuclear_arrival_age_gy FROM stars WHERE star_system_id = ?", (system_id,)).fetchone()
        assert row["yerkes_class"] == "PMS" and row[COLUMN] == pytest.approx(params["nuclear_arrival_age_gy"])
        assert store.load_star_system(conn, system_id).star.nuclear_arrival_age_gy == pytest.approx(params["nuclear_arrival_age_gy"])
    finally:
        conn.close()
