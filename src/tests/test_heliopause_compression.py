"""
A nebula or supernova remnant around a system presses its heliopause in
(`starData.compressed_heliosphere_radius`): the system text says so, and
`query.system_detail` returns the squeezed radius navigation uses.
"""

import pytest

from planetgen.db import query
from planetgen.db import store
from planetgen.physics import constants
from planetgen.generation.star import cloud_pressure_pa, compressed_heliosphere_radius
from tests.test_db_persistence import _placed_nebula, _sector_with_one_system


def test_a_dense_cold_cloud_squeezes_a_sunlike_heliopause_below_an_au():
    squeezed = compressed_heliosphere_radius(85.0, 3000.0, 10.0)
    assert 0.3 < squeezed < 1.0


def test_the_squeeze_follows_the_inverse_square_root_of_pressure():
    low = compressed_heliosphere_radius(100.0, 1000.0, 10.0)
    high = compressed_heliosphere_radius(100.0, 4000.0, 10.0)
    assert high / low == pytest.approx(0.5, rel=0.01)


def test_thin_gas_leaves_the_heliopause_alone():
    assert cloud_pressure_pa(0.1, 8000.0) < constants.ISM_PRESSURE
    assert compressed_heliosphere_radius(100.0, 0.1, 8000.0) == 100.0
    assert compressed_heliosphere_radius(100.0, None, None) == 100.0


def test_hot_remnant_gas_squeezes_through_its_thermal_pressure():
    # A young remnant's million-kelvin gas outpushes the open medium even
    # at a few particles per cm^3.
    assert compressed_heliosphere_radius(100.0, 5.0, 1e6) < 50.0


def test_a_system_inside_a_dense_nebula_reports_its_squeezed_heliopause(mysql_config):
    sector_id = store.save_sector(_sector_with_one_system("Smothered"), config=mysql_config, galaxy_position={
        "center_x_pc": 100.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 100.0,
    })
    conn = store.get_connection(mysql_config)
    try:
        system_id = conn.execute("SELECT id FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()["id"]
        before = query.system_detail(conn, system_id)
        assert before["inside"] is None
        assert before["heliopause_au"] == before["heliopause_open_space_au"] > 0

        with conn:
            nebula_id = _placed_nebula(conn, sector_id, (100.0, 0.0, 0.0), radius_ly=50.0)
            conn.execute("UPDATE nebulae SET density_cm3 = 3000, temperature_k = 10 WHERE id = ?", (nebula_id,))
        detail = query.system_detail(conn, system_id)
        assert detail["inside"]["id"] == nebula_id
        assert detail["inside"]["density_cm3"] == 3000
        assert detail["heliopause_open_space_au"] == pytest.approx(before["heliopause_open_space_au"])
        assert detail["heliopause_au"] == pytest.approx(compressed_heliosphere_radius(
            detail["heliopause_open_space_au"], 3000, 10))
        assert detail["heliopause_au"] < detail["heliopause_open_space_au"] / 10

        system = store.load_star_system(conn, system_id)
        assert system.surrounding_cloud["id"] == nebula_id
        name = detail["inside"]["name"]
        assert f"the gas of {name} around the system presses it in" in system.summary_paragraph()
    finally:
        conn.close()
