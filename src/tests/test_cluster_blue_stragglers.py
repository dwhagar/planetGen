"""
GEN.326: the blue stragglers of the globular clusters (`planetgen/generation/cluster_blue_stragglers.py`, rows BN1 to BN4 of
docs/design/research/gap-filling-plan.md 6.331): a count that follows the core mass to the 0.38 power within Poisson bounds, each in the core,
none outside a globular cluster. No database or Redis: the draws are seeded.
"""

import math

import pytest

from planetgen import tuning
from planetgen.db import store
from planetgen.galaxy import components
from planetgen.generation import cluster_blue_stragglers as blue, component_pass
from planetgen.physics import constants, king, component_profiles

EDGE = 4.0
COLUMN = {name: index for index, name in enumerate(store.BRIGHT_STAR_COLUMNS)}


def _cluster(mass=2.0e5, half=5.0, w0=7, name="NGC Test", center=(8000.0, 0.0, 0.0), cluster_id=3):
    return components.Component(
        kind="cluster", name=name, source="generated", center_pc=center, velocity_kms=(0.0, 0.0, 0.0), epoch_unix=0.0,
        mass_solar=mass, half_mass_radius_pc=half, profile="king", age_gy=12.0, metallicity_fe_h=-1.5, seed=1,
        king_w0=w0, tidal_radius_pc=king.tidal_radius_pc(w0, half), id=cluster_id)


def _rows(component, seed=5, skip=frozenset()):
    return list(blue.blue_straggler_rows(component, EDGE, seed, skip, component_pass.positions_mpc, component_pass._address))


def test_the_count_goes_as_the_core_mass_to_the_exponent_and_is_the_normalization_at_1e4():
    small, big = _cluster(mass=1.0e5), _cluster(mass=1.0e6)
    ratio = blue.expected_count(big) / blue.expected_count(small)
    assert ratio == pytest.approx((blue.core_mass_solar(big) / blue.core_mass_solar(small)) ** tuning.CLUSTER_BLUE_STRAGGLER_CORE_MASS_EXPONENT)
    assert blue.expected_count(_cluster()) == pytest.approx(
        tuning.CLUSTER_BLUE_STRAGGLERS_AT_1E4_CORE_MASS * (blue.core_mass_solar(_cluster()) / 1.0e4) ** tuning.CLUSTER_BLUE_STRAGGLER_CORE_MASS_EXPONENT)


def test_a_clusters_blue_stragglers_are_tens_to_a_few_hundred():
    """The design note's figure: 'typical clusters holding tens to a few hundred'."""
    for mass in (3.0e4, 2.0e5, 1.0e6, 3.0e6):
        assert 10.0 <= blue.expected_count(_cluster(mass=mass)) <= 400.0


def test_the_number_drawn_is_poisson_about_the_expected_count():
    component = _cluster()
    mean = blue.expected_count(component)
    counts = [len(_rows(component, seed=seed)) for seed in range(60)]
    average = sum(counts) / len(counts)
    assert abs(average - mean) < 4.0 * math.sqrt(mean / len(counts)) + 1.0     # a few stars short for the pre-main-sequence draws
    assert min(counts) < max(counts)


def test_each_blue_straggler_sits_in_the_core_in_its_own_cell_and_is_the_clusters():
    component = _cluster()
    core = king.scale_radius_pc(component.king_w0, component.half_mass_radius_pc)
    rows = _rows(component)
    assert rows
    for row in rows:
        point = [row[COLUMN[name]] / 1000.0 for name in ("position_x_mpc", "position_y_mpc", "position_z_mpc")]
        assert math.dist(point, component.center_pc) <= core + 0.01
        assert row[COLUMN["cluster_id"]] == component.id
        assert row[COLUMN["population"]] == "cluster"


def test_a_blue_straggler_is_a_young_main_sequence_star_of_the_mass_range_brighter_than_the_turn_off():
    rows = _rows(_cluster(mass=1.0e6))
    low, high = tuning.BLUE_STRAGGLER_MASS_RANGE_SOL
    assert rows
    for row in rows:
        assert low <= row[COLUMN["initial_mass_sol"]] <= high
        assert 0.0 <= row[COLUMN["age_gy"]] <= tuning.BLUE_STRAGGLER_LIFETIME_GY
        assert row[COLUMN["yerkes_class"]] in ("V", "VI", "IV")
        assert row[COLUMN["luminosity_w"]] / constants.SOLAR_LUMINOSITY > 1.0       # a turn-off star of 0.85 solar masses gives under 1 solar luminosity


def test_the_same_seed_gives_the_same_blue_stragglers_and_another_seed_others():
    component = _cluster()
    assert _rows(component, seed=9) == _rows(component, seed=9)
    assert _rows(component, seed=9) != _rows(component, seed=10)


def test_a_filled_sector_keeps_none_of_them():
    component = _cluster()
    rows = _rows(component)
    addresses = {tuple(row[:3]) for row in rows}
    assert addresses
    assert _rows(component, skip=frozenset(addresses)) == []


def test_a_young_open_cluster_and_a_satellite_draw_none():
    open_cluster = components.Component(
        kind="cluster", name="Open", source="generated", center_pc=(100.0, 0.0, 0.0), velocity_kms=(0.0, 0.0, 0.0), epoch_unix=0.0,
        mass_solar=2.0e3, half_mass_radius_pc=3.0, profile="plummer", age_gy=0.1, metallicity_fe_h=0.0, seed=2,
        tidal_radius_pc=component_profiles.cut_radius_pc("plummer", 3.0), id=4)
    satellite = components.Component(
        kind="satellite", name="Cloud", source="generated", center_pc=(6000.0, 0.0, 0.0), velocity_kms=(0.0, 0.0, 0.0), epoch_unix=0.0,
        mass_solar=2.0e8, half_mass_radius_pc=900.0, profile="exponential", age_gy=6.0, metallicity_fe_h=-0.8, seed=3,
        tidal_radius_pc=component_profiles.cut_radius_pc("exponential", 900.0), id=5)
    for component in (open_cluster, satellite):
        assert blue.expected_count(component) == 0.0
        assert _rows(component) == []


def test_the_blue_stragglers_are_the_passes_own_rows_drawn_with_the_first_slice_only_and_the_stars_of_the_pass_are_unchanged():
    """A cluster drawn in slices has them once; its other stars come from a stream of their own, the blue stragglers' draw takes none of theirs."""
    component = _cluster()
    mine = {tuple(row) for row in _rows(component, seed=7)}
    assert mine
    whole = {tuple(row) for row in component_pass.scatter_blue_stragglers([component], EDGE, 7)}
    first = {tuple(row) for row in component_pass.scatter_blue_stragglers([component], EDGE, 7, slices={("cluster", component.key): (0, 2)})}
    second = {tuple(row) for row in component_pass.scatter_blue_stragglers([component], EDGE, 7, slices={("cluster", component.key): (1, 2)})}
    assert whole == first == mine
    assert not second
