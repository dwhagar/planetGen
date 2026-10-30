# tests/test_gen_star_spin.py

"""
GEN.846: every star the scatter places has a rotation period, drawn from the star's own seed by the cascade a built star
uses (`physics.spin.star_rotation_period_hours`), stored in `bright_stars.rotation_period_hours`, kept by the star that is
built from the row, and shown on the uncharted star page and in the API.
"""

import math
import random

import pytest

from planetgen.db import store
from planetgen.generation import bright_stars as brightStars
from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.physics import constants, spin
from tests.test_api import client  # noqa: F401
from tests.test_bright_star_scatter import _row_dict, _scatter
from tests.test_galactic_motion import SCATTER_RADIUS_PC, _insert_scatter_point, _seed_scatter_grid
from tests.test_uncharted_phenomenon import db_app  # noqa: F401

SUN = dict(mass_kg=constants.SOLAR_MASS_TO_KG, radius_km=6.957e5, temperature_k=5772.0, age_gy=4.6, yerkes_class="V")


def test_a_star_always_spins_with_the_same_period_from_its_seed():
    assert brightStars.rotation_period_hours(99, SUN) == brightStars.rotation_period_hours(99, SUN)
    assert brightStars.rotation_period_hours(99, SUN) != brightStars.rotation_period_hours(100, SUN)


def test_the_period_is_the_cascades_for_each_kind_of_star():
    days = [brightStars.rotation_period_hours(seed, SUN) / 24.0 for seed in range(200)]
    assert 10.0 < sorted(days)[100] < 60.0                    # a Sun-like dwarf turns in tens of days
    white_dwarf = dict(mass_kg=1.2e30, radius_km=7000.0, temperature_k=12000.0, age_gy=4.0, yerkes_class="VII")
    assert all(brightStars.rotation_period_hours(seed, white_dwarf) > 0.0 for seed in range(50))
    giant = dict(mass_kg=2.0e30, radius_km=3.0e7, temperature_k=4500.0, age_gy=1.0, yerkes_class="III")
    assert all(brightStars.rotation_period_hours(seed, giant) > 100.0 * 24.0 for seed in range(50))


def test_the_period_never_beats_breakup():
    fast = dict(mass_kg=20 * constants.SOLAR_MASS_TO_KG, radius_km=6.0e6, temperature_k=30000.0, age_gy=0.005,
                yerkes_class="V")
    cap = spin.period_hours_from_speed(fast["radius_km"],
                                       spin.BREAKUP_FRACTION * spin.breakup_speed_kms(fast["mass_kg"], fast["radius_km"]))
    assert all(brightStars.rotation_period_hours(seed, fast) >= cap for seed in range(100))


def test_every_scattered_star_has_a_period_and_the_scatter_draws_what_it_drew_before():
    rows = [_row_dict(row) for row in _scatter()]
    assert rows and all(row["rotation_period_hours"] > 0.0 and math.isfinite(row["rotation_period_hours"]) for row in rows)
    for row in rows:
        assert row["rotation_period_hours"] == brightStars.rotation_period_hours(row["seed"], row)
    rng = random.Random(5)
    before = rng.getstate()
    brightStars._row(0, 0, 0, (0, 0, 0), "young", {
        "type": "G2V", "yerkes_class": "V", "mass_kg": SUN["mass_kg"], "radius_km": SUN["radius_km"],
        "temperature_k": SUN["temperature_k"], "luminosity_w": 3.8e26, "age_gy": 4.6, "lifespan_gy": 10.0,
        "initial_mass_sol": 1.0, "phase_end_age_gy": 10.0, "metallicity_fe_h": 0.0}, rng)
    witness = random.Random()
    witness.setstate(before)
    witness.getrandbits(63)
    assert rng.getstate() == witness.getstate()               # the row took one draw, the seed, and no more


def test_a_star_built_from_a_stored_row_keeps_the_rows_period():
    row = _row_dict(_scatter()[0])
    star = Star.from_params(brightStars.star_params(row), SystemConfig(), name="Test")
    assert star.rotation_period_hours == row["rotation_period_hours"]
    old_row = dict(row, rotation_period_hours=None)           # a row made before v176 draws its own
    star = Star.from_params(brightStars.star_params(old_row), SystemConfig(), name="Test")
    assert star.rotation_period_hours > 0.0


def test_the_uncharted_star_page_and_the_api_give_the_stored_period(db_app, client, mysql_config):  # noqa: F811
    _seed_scatter_grid(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            row_id = _insert_scatter_point(conn, "bright_stars", (SCATTER_RADIUS_PC, 0.0, 0.0))
            conn.execute("UPDATE bright_stars SET rotation_period_hours = 36.5 WHERE id = ?", (row_id,))
            bare = _insert_scatter_point(conn, "bright_stars", (SCATTER_RADIUS_PC + 3.0, 0.0, 0.0))
    finally:
        conn.close()
    assert client.get(f"/api/uncharted-systems/{row_id}").get_json()["rotation_period_hours"] == 36.5
    page = db_app.test_client().get(f"/systems/uncharted/{row_id}").get_data(as_text=True)
    assert "Rotation period" in page
    assert client.get(f"/api/uncharted-systems/{bare}").get_json()["rotation_period_hours"] is None
    assert "Rotation period" not in db_app.test_client().get(f"/systems/uncharted/{bare}").get_data(as_text=True)
