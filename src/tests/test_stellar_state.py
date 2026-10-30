# tests/test_stellar_state.py

"""GEN.786 and GEN.785: the definitions of a star's evolutionary state and the age since it stopped burning (`physics/stellar_state.py`)."""

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
    assert stellar_state.decay_age_yr("white-dwarf", 10.0, 7.5) == pytest.approx(2.5e9)
    assert stellar_state.decay_age_yr("white-dwarf", 7.0, 7.5) == 0.0
    for state in stellar_state.STATES:
        if state != "white-dwarf":
            assert stellar_state.decay_age_yr(state, 10.0, 7.5) is None
    assert stellar_state.decay_age_yr("white-dwarf", 10.0, None) is None
