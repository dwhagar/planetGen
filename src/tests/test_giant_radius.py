# tests/test_giant_radius.py

"""
GEN.460: the giant-planet radius from the Fortney tables, the Sestovic inflation and the Chachan heavy-element mass
(`planetgen.physics.giant_radius`): the tables are read at their own grid points, the radius falls with core mass and age,
is held at the edges, inflation starts at its threshold, and the heavy-element mass never reaches the mass.
"""

import math

import pytest

from planetgen.physics import constants, giant_radius as gr, giant_radius_tables as tables
from planetgen.util import draw

MJ_EARTH = constants.JUPITER_MASS_TO_KG / constants.EARTH_MASS_TO_KG


def test_every_table_cell_is_returned_at_its_grid_point():
    for ai, age in enumerate(tables.AGES_GY):
        for di, distance in enumerate(tables.DISTANCES_AU):
            for ci, core in enumerate(tables.CORE_MASSES_EARTH):
                for mi, mass_mj in enumerate(tables.MASSES_MJ):
                    expected = tables.RADII_RJ[ai][di][ci][mi]
                    if expected is None:
                        continue
                    got = gr.radius_rj(mass_mj * MJ_EARTH, core, 1.0 / distance ** 2, age)
                    assert got == pytest.approx(expected, abs=1e-9), (age, distance, core, mass_mj)


def test_the_table_is_in_the_order_physics_requires():
    for ai in range(3):
        for di in range(5):
            for mi in range(len(tables.MASSES_MJ)):
                column = [tables.RADII_RJ[ai][di][c][mi] for c in range(5)]
                valid = [r for r in column if r is not None]
                assert valid == sorted(valid, reverse=True), (ai, di, mi)
    for di in range(5):
        for ci in range(5):
            for mi in range(len(tables.MASSES_MJ)):
                by_age = [tables.RADII_RJ[a][di][ci][mi] for a in range(3)]
                if None not in by_age:
                    assert by_age == sorted(by_age, reverse=True), (di, ci, mi)


def test_a_heavier_core_a_greater_age_and_a_longer_distance_shrink_a_giant():
    base = gr.radius_rj(317.8, 20.0, 1.0, 1.0)
    assert gr.radius_rj(317.8, 60.0, 1.0, 1.0) < base
    assert gr.radius_rj(317.8, 20.0, 1.0, 4.5) < base
    assert gr.radius_rj(317.8, 20.0, 0.01, 1.0) < base
    assert gr.radius_rj(317.8, 20.0, 100.0, 1.0) > base


def test_a_jupiter_is_about_a_jupiter_radius_at_five_au_and_old():
    assert gr.radius_rj(317.8, 40.0, 1.0 / 5.2 ** 2, 4.5) == pytest.approx(1.0, abs=0.05)


def test_the_edges_are_held_in_flux_and_age_and_continued_in_mass():
    assert gr.radius_rj(317.8, 20.0, 1e6, 1.0) == gr.radius_rj(317.8, 20.0, 1e8, 1.0)
    assert gr.radius_rj(317.8, 20.0, 1.0, 20.0) == gr.radius_rj(317.8, 20.0, 1.0, 4.5)
    top = tables.MASSES_MJ[-1] * MJ_EARTH
    assert gr.in_range(top * 1.1) and not gr.in_range(top * 1.5)
    assert abs(gr.radius_rj(top * 1.1, 100.0, 1.0, 4.5) - gr.radius_rj(top, 100.0, 1.0, 4.5)) < 0.05


def test_a_mass_below_the_table_or_a_bad_flux_raises():
    with pytest.raises(gr.GiantRadiusError):
        gr.radius_rj(10.0, 5.0, 1.0, 1.0)
    with pytest.raises(gr.GiantRadiusError):
        gr.radius_rj(300.0, 20.0, 0.0, 1.0)


def test_a_core_as_heavy_as_the_planet_is_held_at_the_heaviest_the_table_has():
    low = tables.MASSES_MJ[0] * MJ_EARTH
    assert gr.radius_rj(low, 17.0, 1.0, 4.5) == gr.radius_rj(low, 10.0, 1.0, 4.5)


def test_inflation_starts_at_the_threshold_flux_and_not_below_a_third_of_a_jupiter():
    mj = MJ_EARTH
    assert gr.inflation_rj(0.2 * mj, 1e7) == 0.0
    assert gr.inflation_rj(0.5 * mj, 10 ** 5.5) == 0.0
    assert gr.inflation_rj(0.5 * mj, 10 ** 6.5) == pytest.approx(0.70 * (6.5 - 5.52))
    assert gr.inflation_rj(1.5 * mj, 10 ** 6.5) == pytest.approx(0.52 * (6.5 - 5.82))
    assert gr.inflation_rj(5.0 * mj, 10 ** 6.5) == pytest.approx(0.22 * (6.5 - 5.2))
    assert gr.inflation_rj(0.5 * mj, 0.0) == 0.0


def test_the_core_the_table_is_read_at_leaves_out_the_envelopes_heavy_elements():
    assert gr.table_core_mass_earth(317.8, 42.0, 0.0) == 42.0
    assert gr.table_core_mass_earth(317.8, 42.0, 0.04) == pytest.approx((42.0 - 0.04 * 317.8) / 0.96)
    assert gr.table_core_mass_earth(100.0, 3.0, 0.1) == 0.0


def test_the_tables_cover_a_jupiter_but_not_an_ice_giant():
    assert gr.covers(317.8, 42.0) and gr.covers(95.0, 22.0) and gr.covers(1000.0, 100.0)
    assert not gr.covers(17.1, 14.6)   # Neptune is mostly heavy elements
    assert not gr.covers(47.0, 30.0)


def test_the_flux_of_the_earth_in_watts_is_the_solar_constant():
    assert gr.flux_w_m2(1.0) == pytest.approx(1361.0, rel=0.01)


def test_a_hot_jupiter_is_bigger_than_a_cool_one():
    cool = gr.radius_km(317.8, 40.0, 1.0, 4.5)
    hot = gr.radius_km(317.8, 40.0, 3000.0, 4.5)
    assert hot > cool * 1.15


def test_the_predicted_heavy_mass_is_the_chachan_law():
    assert gr.predicted_heavy_mass_earth(10.0) == 10.0
    assert gr.predicted_heavy_mass_earth(317.8) == pytest.approx(14.7 + 0.09 * (317.8 - 14.7))
    assert gr.predicted_heavy_mass_earth(14.7) == 14.7


@pytest.mark.parametrize("mass", [17.0, 25.0, 95.0, 317.8, 1000.0, 4000.0])
def test_a_drawn_heavy_mass_leaves_an_envelope_and_is_deterministic(mass):
    values = []
    for seed in range(1, 40):
        with draw.bound(seed):
            value = gr.heavy_mass_earth(mass)
        assert gr.MIN_HEAVY_FRACTION * mass <= value <= gr.MAX_HEAVY_FRACTION * mass
        values.append(value)
    assert len(set(values)) > 5
    with draw.bound(3):
        again = gr.heavy_mass_earth(mass)
    assert again == values[2]


def test_the_scatter_is_about_the_papers_sixty_six_percent():
    mass = 600.0
    predicted = gr.predicted_heavy_mass_earth(mass)
    values = []
    for seed in range(1, 600):
        with draw.bound(seed):
            values.append(gr.heavy_mass_earth(mass))
    mean = sum(values) / len(values)
    spread = math.sqrt(sum((v - mean) ** 2 for v in values) / len(values))
    assert mean == pytest.approx(predicted, rel=0.15)
    assert 0.35 * predicted < spread < 0.7 * predicted  # a normal cut at the bounds is narrower than 0.66
