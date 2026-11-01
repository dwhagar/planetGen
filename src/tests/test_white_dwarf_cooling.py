# tests/test_white_dwarf_cooling.py

"""
GEN.742: a white dwarf's temperature, luminosity, radius and crystallised core come from cooling sequences by mass and
atmosphere (`datasets.white_dwarf_models.cooling_state`), not from a power law of the cooling age.

The tests read the downloaded [R0308] sequences from the machine's cache folder (like the MIST isochrones, a generator that
makes a white dwarf needs them) and skip on a machine that has not downloaded them; the last two tests are about the files being absent.
"""

import math
import random

import numpy as np
import pytest

from planetgen import tuning
from planetgen.datasets import fetch, white_dwarf_models
from planetgen.physics import constants, stellar_composition, stellar_evolution


@pytest.fixture(autouse=True)
def _needs_the_sequences(request):
    if "without_the_sequences" not in request.node.name and fetch.locate(white_dwarf_models.DATASET).missing:
        pytest.skip("the white-dwarf sequences are not in this machine's cache folder")


def _luminosity_from_radius_and_temperature(state):
    radius_m = state.radius_solar * constants.SOLAR_RADIUS_M
    return 4 * math.pi * radius_m ** 2 * constants.STEFAN_BOLTZMANN_CONSTANT * state.teff_k ** 4 / constants.SOLAR_LUMINOSITY


def test_a_white_dwarf_cools_and_its_core_crystallises():
    states = [white_dwarf_models.cooling_state(0.6, age) for age in (0.001, 0.01, 0.1, 1.0, 3.0, 8.0)]
    assert [state.teff_k for state in states] == sorted((state.teff_k for state in states), reverse=True)
    assert [state.luminosity_sol for state in states] == sorted((state.luminosity_sol for state in states), reverse=True)
    assert states[0].crystallised_fraction == 0.0 and states[-1].crystallised_fraction > 0.9


def test_the_luminosity_obeys_stefan_boltzmann_for_the_radius_and_temperature():
    for mass, age in ((0.5, 0.3), (0.62, 2.0), (0.95, 5.0), (1.2, 0.5)):
        state = white_dwarf_models.cooling_state(mass, age)
        assert state.luminosity_sol == pytest.approx(_luminosity_from_radius_and_temperature(state), rel=0.03)


def test_a_heavier_white_dwarf_is_smaller_and_both_layers_are_read():
    assert white_dwarf_models.cooling_state(1.1, 1.0).radius_solar < white_dwarf_models.cooling_state(0.6, 1.0).radius_solar
    thick = white_dwarf_models.cooling_state(0.6, 3.0, "thick")
    thin = white_dwarf_models.cooling_state(0.6, 3.0, "thin")
    assert thick.radius_solar != thin.radius_solar and thick.radius_solar == pytest.approx(thin.radius_solar, rel=0.1)


def test_a_mass_between_two_grid_masses_lies_between_them():
    low, mid, high = (white_dwarf_models.cooling_state(m, 1.0) for m in (0.6, 0.65, 0.7))
    assert min(low.radius_solar, high.radius_solar) <= mid.radius_solar <= max(low.radius_solar, high.radius_solar)
    assert min(low.teff_k, high.teff_k) <= mid.teff_k <= max(low.teff_k, high.teff_k)


def test_one_older_than_the_last_model_is_held_at_it_and_says_so():
    held = white_dwarf_models.cooling_state(1.2, 30.0)
    assert held.extrapolated
    assert held.teff_k == pytest.approx(white_dwarf_models.sequence(1.2).teff_k[-1], rel=1e-6)
    assert not white_dwarf_models.cooling_state(1.2, 0.5).extrapolated


def test_one_heavier_than_the_grid_reads_the_heaviest_sequence_and_says_so():
    heavy = white_dwarf_models.cooling_state(1.35, 1.0)
    top = white_dwarf_models.cooling_state(white_dwarf_models.CRYSTAL_MASS_SOL, 1.0)
    assert heavy.mass_clamped and not top.mass_clamped
    assert heavy.teff_k == top.teff_k and heavy.radius_solar == top.radius_solar


def test_a_non_positive_mass_or_age_is_refused():
    with pytest.raises(ValueError):
        white_dwarf_models.cooling_state(0.0, 1.0)
    with pytest.raises(ValueError):
        white_dwarf_models.cooling_state(0.6, 0.0)


def test_an_evolved_star_older_than_its_cooling_time_is_a_white_dwarf_from_the_sequences():
    state = stellar_evolution.evolve_star(2.0, 9.0, random.Random(5), fe_h=0.0)
    assert state["yerkes_class"] == "VII"
    envelope = "thick" if state["white_dwarf_atmosphere"] == "DA" else "thin"
    mass = state["mass_sol"]
    # it cooled for less than its whole age, and its radius is the sequence's at its own temperature
    assert white_dwarf_models.cooling_state(mass, 9.0, envelope).teff_k <= state["temperature_k"] <= \
        white_dwarf_models.cooling_state(mass, tuning.WD_MIN_COOLING_AGE_GY, envelope).teff_k
    assert state["radius_km"] / (constants.SOLAR_RADIUS_M / 1000.0) == pytest.approx(
        white_dwarf_models.radius_solar(mass, state["temperature_k"], envelope), rel=0.02)
    radius_m = state["radius_km"] * 1000.0
    watts = 4 * math.pi * radius_m ** 2 * constants.STEFAN_BOLTZMANN_CONSTANT * state["temperature_k"] ** 4
    assert state["luminosity_sol"] == pytest.approx(watts / constants.SOLAR_LUMINOSITY, rel=0.03)
    assert 0.0 <= state["crystallised_fraction"] <= 1.0


def test_the_atmosphere_follows_the_published_share_of_hydrogen_atmospheres():
    rng = random.Random(11)
    kinds = [stellar_evolution.evolve_star(2.0, 9.0, rng, fe_h=0.0)["white_dwarf_atmosphere"] for _ in range(600)]
    assert kinds.count("DA") / len(kinds) == pytest.approx(tuning.STELLAR_WD_DA_SHARE, abs=0.06)


def test_without_the_sequences_a_run_that_disregards_missing_data_takes_the_labelled_power_law(monkeypatch, tmp_path):
    monkeypatch.setenv(fetch.CACHE_DIR_ENV, str(tmp_path))
    with fetch.disregard_missing_data():
        state = stellar_evolution.evolve_star(2.0, 9.0, random.Random(5), fe_h=0.0)
        assert ("white_dwarf_models", "the Mestel power law of the white-dwarf cooling") in fetch.stand_ins_used()
    assert "white_dwarf_atmosphere" not in state
    assert state["radius_km"] == stellar_evolution.white_dwarf_radius_km(state["mass_sol"])


def test_without_the_sequences_a_run_refuses_a_white_dwarf(monkeypatch, tmp_path):
    monkeypatch.setenv(fetch.CACHE_DIR_ENV, str(tmp_path))
    with pytest.raises(fetch.DataMissing):
        stellar_evolution.evolve_star(2.0, 9.0, random.Random(5), fe_h=0.0)


def test_a_white_dwarf_heavier_than_the_sequences_still_has_a_radius_that_keeps_shrinking():
    """The initial-final mass relation reaches 1.33 solar masses and the cooling sequences stop at 1.30."""
    radii = [stellar_composition._grid_radius_solar(mass, 20000.0, "thick") for mass in (1.25, 1.30, 1.31, 1.33)]
    assert all(radius is not None and 0.003 < radius < 0.05 for radius in radii)
    assert radii == sorted(radii, reverse=True)
