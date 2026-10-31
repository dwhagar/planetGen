# tests/test_stellar_state.py

"""GEN.786 and GEN.785: the definitions of a star's evolutionary state and the age since it stopped burning (`physics/stellar_state.py`)."""

import types

import pytest

from planetgen.physics import pre_main_sequence, stellar_evolution, stellar_state


def test_every_state_fits_the_state_column():
    assert all(len(state) <= 16 for state in stellar_state.STATES)
    assert len(set(stellar_state.STATES)) == len(stellar_state.STATES)


@pytest.mark.parametrize("yerkes, phase, state", [
    ("protostar", None, "protostar"), (pre_main_sequence.PRE_MAIN_SEQUENCE, None, "pre-main-seq"), ("V", None, "main-sequence"), ("VI", None, "main-sequence"),
    ("IV", None, "subgiant"), ("III", None, "red-giant"), ("II", None, "red-giant"), ("III", "horizontal-branch", "helium-burning"),
    ("III", "asymptotic-giant", "asymptotic-giant"), ("IB", None, "supergiant"), ("IA+", None, "hypergiant"), ("W", None, "wolf-rayet"),
    ("VII", None, "white-dwarf"), ("D", None, "white-dwarf")])
def test_a_class_and_phase_name_one_state(yerkes, phase, state):
    assert stellar_state.evolutionary_state(yerkes, phase) == state


def test_a_brown_dwarf_is_a_brown_dwarf_whatever_class_it_was_given():
    assert stellar_state.evolutionary_state("V", brown_dwarf=True) == "brown-dwarf"


def test_a_class_the_table_lacks_is_refused():
    with pytest.raises(ValueError):
        stellar_state.evolutionary_state("X")


@pytest.mark.parametrize("star_type, phase", [
    ("G8III Yellow Horizontal Branch Giant Star", "horizontal-branch"), ("K0III Orange Core Helium Burning Giant Star", "core-helium-burning"),
    ("M2III Red Asymptotic Giant Star", "asymptotic-giant"), ("F7II Yellow-White Classical Cepheid Star", "classical-cepheid"),
    ("G2V Yellow Main Sequence Star", None), ("K1VII Orange White Dwarf Star", None), (None, None)])
def test_the_phase_is_read_back_from_the_type_text(star_type, phase):
    assert stellar_state.phase_of(star_type) == phase


@pytest.mark.parametrize("yerkes, phase, state", [
    ("III", "core-helium-burning", "helium-burning"), ("II", "classical-cepheid", "cepheid"), ("II", "asymptotic-giant", "asymptotic-giant"),
    ("II", None, "red-giant")])
def test_the_phases_of_the_tracks_name_a_state(yerkes, phase, state):
    assert stellar_state.evolutionary_state(yerkes, phase) == state


def test_every_state_the_evolution_gives_has_a_word():
    """A sweep of `evolve_star` over masses and ages: every class and phase it returns maps to a state in `STATES`."""
    pytest.importorskip("planetgen.datasets.isochrones")
    from planetgen.datasets import fetch, isochrones
    if fetch.locate(isochrones.DATASET).missing:
        pytest.skip("the MIST isochrones are not in this machine's cache folder")
    seen = set()
    for mass in (0.3, 0.8, 1.0, 1.5, 3.0, 8.0, 15.0, 40.0):
        for age in (0.001, 0.1, 1.0, 5.0, 9.0, 11.0, 13.0):
            state = stellar_evolution.evolve_star(mass, age, fe_h=0.0)
            if state is not None:
                seen.add(stellar_state.evolutionary_state(state["yerkes_class"], state.get("phase")))
    assert seen <= set(stellar_state.STATES) and "main-sequence" in seen and len(seen) >= 4


def test_only_a_white_dwarf_has_an_age_since_it_stopped_burning():
    assert stellar_state.decay_age_yr("white-dwarf", 2.5) == pytest.approx(2.5e9)
    for state in stellar_state.STATES:
        if state != "white-dwarf":
            assert stellar_state.decay_age_yr(state, 2.5) is None
    assert stellar_state.decay_age_yr("white-dwarf", None) is None


def _star(yerkes, star_type, age, initial_mass=1.0, fe_h=0.0):
    return types.SimpleNamespace(yerkes_class=yerkes, type=star_type, age=age, initial_mass_sol=initial_mass, metallicity_fe_h=fe_h)


def test_set_state_gives_a_white_dwarf_its_cooling_age_and_a_burning_star_none():
    dwarf = _star("VII", "K1VII Orange White Dwarf Star", 10.0)
    stellar_state.set_state(dwarf)
    assert dwarf.state == "white-dwarf"
    assert dwarf.decay_age_yr == pytest.approx(1.0e9 * stellar_evolution.cooling_age_gy(1.0, 10.0, 0.0))
    assert 0.0 < dwarf.decay_age_yr < 10.0e9
    sun = _star("V", "G2V Yellow Main Sequence Star", 4.6)
    stellar_state.set_state(sun)
    assert (sun.state, sun.decay_age_yr) == ("main-sequence", None)


def test_a_white_dwarf_of_the_evolution_has_the_cooling_age_the_evolution_cooled_it_for():
    """`evolve_star` and `cooling_age_gy` use one age: a progenitor that ended at 1 Gy leaves a dwarf 9 Gy old at 10 Gy."""
    end = stellar_evolution.post_main_sequence.timeline(1.0, 0.0).end_gy
    assert stellar_evolution.cooling_age_gy(1.0, end + 3.0, 0.0) == pytest.approx(3.0)
    assert stellar_evolution.cooling_age_gy(1.0, end - 1.0, 0.0) == stellar_evolution.tuning.WD_MIN_COOLING_AGE_GY
