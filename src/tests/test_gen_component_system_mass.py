# tests/test_gen_component_system_mass.py -- GEN.598: a cluster's or satellite's systems are its mass over the mean mass of
# its own systems, worked out from its age with the same mass function, remnants and multiplicity the galaxy's are (GEN.545).

import math
import statistics

import pytest

from planetgen.galaxy import components
from planetgen.generation import multistar, system_mass
from planetgen.physics import imf, remnant_mass, remnant_progenitors, stellar_evolution
from planetgen.util import draw


def _component(age_gy, kind="cluster", mass=2.0e5):
    return components.Component(
        kind=kind, name="Test", source="generated", center_pc=(0.0, 0.0, 0.0), velocity_kms=(0.0, 0.0, 0.0), epoch_unix=0.0,
        mass_solar=mass, half_mass_radius_pc=5.0, profile="plummer", age_gy=age_gy, metallicity_fe_h=-0.5, seed=1,
        tidal_radius_pc=50.0)


@pytest.mark.parametrize("age_gy", [0.03, 0.3, 3.0, 12.0])
def test_a_components_mean_is_a_seeded_draw_of_systems_of_its_age(age_gy):
    rng = draw.Stream(f"gen598-systems-{age_gy}")
    found = []
    for _ in range(4000):
        mass = imf.sample(rng=rng)
        if stellar_evolution.evolve_star(mass, age_gy, rng) is None:
            kind = "black-hole" if remnant_progenitors.fails(mass) else "neutron-star"
            found.append(remnant_mass.draw_mass_sol(kind, rng=rng))
            continue
        found.append(multistar.draw_member_tree(multistar.Placed("star", mass), age_gy * 1.0e9, rng=rng).system_mass_solar)
    error = statistics.stdev(found) / math.sqrt(len(found))
    assert abs(statistics.fmean(found) - _component(age_gy).mean_system_mass_solar) < 4.0 * error + 0.02 * statistics.fmean(found)


def test_the_systems_are_the_mass_over_the_mean_and_a_younger_component_holds_fewer():
    old, young = _component(12.0), _component(0.03)
    assert old.mean_system_mass_solar == pytest.approx(system_mass.mean_system_mass_at_age_solar(12.0))
    assert old.expected_systems == pytest.approx(old.mass_solar / old.mean_system_mass_solar)
    assert young.mean_system_mass_solar > old.mean_system_mass_solar
    assert young.expected_systems < old.expected_systems


def _by_turn_off(age_gy):
    """The mean system mass the long way: the turn-off by the full search, then the sums over every star."""
    import numpy as np

    number, others, placed_states, companion_states = system_mass._placed_states()
    turn_off = system_mass._turn_off_sol(age_gy)
    placed = float(np.sum(number * system_mass._present_mass(placed_states, turn_off)) / np.sum(number))
    further = system_mass._present_mass(companion_states, turn_off).mean(axis=1)
    return placed + float(np.sum(number * others * further) / np.sum(number))


@pytest.mark.parametrize("age_gy", [0.0036, 0.05, 1.03, 2.0, 5.2, 12.0, 20.0])
def test_the_mean_by_the_age_table_is_the_mean_by_the_full_turn_off_search(age_gy):
    """The sorted table with running sums gives what the pass over every star gave, to rounding (a live galaxy's
    page view waited minutes for the long way once per open cluster)."""
    assert system_mass.mean_system_mass_at_age_solar(age_gy) == pytest.approx(_by_turn_off(age_gy), rel=1e-10)


def test_a_second_mean_for_the_same_age_does_no_stellar_evolution(monkeypatch):
    system_mass.mean_system_mass_at_age_solar(7.123456)
    calls = []
    monkeypatch.setattr(stellar_evolution, "remnant_formation_time_gy", lambda *args: calls.append(args) or 1.0)
    system_mass.mean_system_mass_at_age_solar(7.123456)
    assert calls == []


def test_the_turn_off_search_that_stops_early_lands_among_the_same_masses_as_the_full_one():
    births = system_mass._age_table()[0]
    for age_gy in (0.004, 0.3, 1.03, 2.2, 9.0):
        early = stellar_evolution.youngest_progenitor_sol(
            age_gy * 1.0e9, settled=lambda low, high: births.searchsorted(low, side="right") == births.searchsorted(high, side="right"))
        full = stellar_evolution.youngest_progenitor_sol(age_gy * 1.0e9)
        assert births.searchsorted(early, side="right") == births.searchsorted(full, side="right")


def test_means_the_update_stored_are_read_back_without_any_stellar_evolution(tmp_path, monkeypatch):
    """`update.sh` works the cluster ages out once (`save_means`); a web process then looks them up (the daemon is stopped
    at 60 seconds and a galaxy has tens of thousands of ages)."""
    monkeypatch.setattr(system_mass, "_store_path", lambda: str(tmp_path / system_mass.STORE_FILE))
    system_mass._stored_means.cache_clear()
    ages = (0.123457, 4.567891)
    worked_out = {age: system_mass.mean_system_mass_at_age_solar(age) for age in ages}
    seen = []
    assert system_mass.save_means(ages, progress=lambda done, total: seen.append((done, total))) == 2
    assert seen == [(1, 2), (2, 2)]
    assert (tmp_path / system_mass.STORE_FILE).stat().st_mode & 0o044 == 0o044           # the web app only reads it
    system_mass.mean_system_mass_at_age_solar.cache_clear()
    monkeypatch.setattr(stellar_evolution, "youngest_progenitor_sol", lambda *args, **kwargs: pytest.fail("worked out again"))
    try:
        for age, mean in worked_out.items():
            assert system_mass.mean_system_mass_at_age_solar(age) == mean
    finally:
        system_mass._stored_means.cache_clear()
        system_mass.mean_system_mass_at_age_solar.cache_clear()


def test_means_stored_by_another_release_are_not_used(tmp_path, monkeypatch):
    monkeypatch.setattr(system_mass, "_store_path", lambda: str(tmp_path / system_mass.STORE_FILE))
    (tmp_path / system_mass.STORE_FILE).write_text('{"release": "0.0.0", "means": {"1.5": 99.0}}')
    system_mass._stored_means.cache_clear()
    system_mass.mean_system_mass_at_age_solar.cache_clear()
    try:
        assert system_mass.mean_system_mass_at_age_solar(1.5) != 99.0
    finally:
        system_mass._stored_means.cache_clear()
        system_mass.mean_system_mass_at_age_solar.cache_clear()
