# tests/test_star_metallicity_storage.py

"""
GEN.158: a star's [Fe/H] is stored (`stars.metallicity_fe_h`, `bright_stars.metallicity_fe_h`, migration 0095).
A scattered star draws its own at scatter time from its population and its place on a stream of its own seed;
the fill that builds its system keeps it; a binary's two stars share one; the sector map shows it. The model itself
(`physics/metallicity.py`) is in `test_gen_metallicity.py`.
"""

import statistics

import pytest

from planetgen import tuning
from planetgen.db import store
from planetgen.generation import bright_stars as brightStars
from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.generation.system import StarSystem
from planetgen.util import draw
from planetgen.web.maps import starmap
from tests.test_bright_star_scatter import _row_dict, _scatter


def test_a_scattered_row_carries_a_drawn_fe_h_in_range():
    rows = _scatter(seed=7)
    assert rows and all(len(row) == len(store.BRIGHT_STAR_COLUMNS) for row in rows)
    low, high = tuning.METALLICITY_RANGE_FE_H
    for row in map(_row_dict, rows):
        assert row["metallicity_fe_h"] is not None and low <= row["metallicity_fe_h"] <= high, row


def test_the_same_seed_draws_the_same_fe_h_and_another_draws_another():
    first = [_row_dict(row)["metallicity_fe_h"] for row in _scatter(seed=11)]
    assert first == [_row_dict(row)["metallicity_fe_h"] for row in _scatter(seed=11)]
    assert first != [_row_dict(row)["metallicity_fe_h"] for row in _scatter(seed=12)]


def test_halo_stars_are_far_poorer_in_metals_than_disc_stars():
    kpc_mpc = 1_000_000
    rng = draw.Stream("test-halo-disc")
    disc = [brightStars.spot_fe_h("young", [((8 * kpc_mpc, 0, 0),)], rng)(0) for _ in range(300)]
    halo = [brightStars.spot_fe_h("halo", [((30 * kpc_mpc, 0, 0),)], rng)(0) for _ in range(300)]
    assert statistics.mean(disc) - statistics.mean(halo) > 1.5


def test_a_cluster_star_has_no_value_until_its_cluster_gives_one():
    assert brightStars.star_fe_h("cluster", {"metallicity_fe_h": -0.4}) is None
    assert brightStars.spot_fe_h("cluster", [((1000, 0, 0),)], draw.Stream("t")) is None


def test_the_stored_value_comes_back_in_the_stars_params():
    row = {"star_type": "G2V Yellow Main Sequence Star", "yerkes_class": "V", "mass_kg": 2e30, "radius_km": 7e5,
           "temperature_k": 5772.0, "luminosity_w": 3.8e26, "age_gy": 4.6, "lifespan_gy": 10.0,
           "initial_mass_sol": 1.0, "phase_end_age_gy": 9.0, "metallicity_fe_h": -0.7}
    assert brightStars.star_params(row)["metallicity_fe_h"] == -0.7
    assert brightStars.star_params({**row, "metallicity_fe_h": None})["metallicity_fe_h"] is None


def test_a_star_built_from_a_stored_row_keeps_its_fe_h_and_one_without_draws_its_own():
    draw.set_run_seed(3)
    config = SystemConfig()
    built = Star(config, name="Seed")
    params = {"type": built.type, "yerkes_class": built.yerkes_class, "mass_kg": built.mass,
              "radius_km": built.radius, "temperature_k": built.temperature, "luminosity_w": built.luminosity,
              "age_gy": built.age, "lifespan_gy": built.lifespan, "initial_mass_sol": built.initial_mass_sol,
              "phase_end_age_gy": built.phase_end_age_gy}
    kept = Star.from_params({**params, "metallicity_fe_h": -0.7}, config, name="Kept")
    assert kept.metallicity_fe_h == -0.7
    drawn = Star.from_params({**params, "metallicity_fe_h": None}, config, name="Drawn")
    assert isinstance(drawn.metallicity_fe_h, float)


@pytest.mark.parametrize("wide", [False, True])
def test_both_stars_of_a_binary_share_one_fe_h(wide):
    for seed in range(6):
        draw.set_run_seed(seed)
        cfg = SystemConfig()
        cfg.BINARY_SYSTEM = True
        cfg.WIDE_BINARY = wide
        system = StarSystem(system_config=cfg)
        assert system.primary_star.metallicity_fe_h == system.secondary_star.metallicity_fe_h, seed


def test_the_map_writes_the_fe_h_with_a_true_minus_sign_and_leaves_an_old_row_blank():
    assert starmap._fe_h_text(-0.123) == "−0.12 dex"
    assert starmap._fe_h_text(0.3) == "+0.30 dex"
    assert starmap._fe_h_text(None) is None
    star = {"star_type": "G2V", "temperature_k": 5772.0, "radius_km": 7e5, "luminosity_w": 3.8e26,
            "temp_display": "5772 K", "name": "Sol", "metallicity_fe_h": -0.4}
    system = {"id": "x", "name": "Sol", "quadrant": "I", "location": ""}
    link = lambda kind, **kwargs: "/system/x"  # noqa: E731
    assert starmap._star_data(link, system, star, 0, 0, 0)["metallicity"] == "−0.40 dex"
    assert starmap._star_data(link, system, {**star, "metallicity_fe_h": None}, 0, 0, 0)["metallicity"] is None
