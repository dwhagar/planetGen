"""
The runs of the mantle box that a later age takes up (PERF.289, `planetgen.physics.interior_heat._Run`): the history of the air asks for the heat of one body at
a hundred ages, and each is the result of the box run from the formation; a run made for one age holds, up to the step that landed on its age, the steps of every
younger age, so the later ages take it up and the result is the one a run from the formation gives, bit for bit, in whatever order the ages come.
"""

import random

import pytest

from planetgen.physics import constants, interior_heat as ih, radiogenic, tidal_heating

from tests.test_interior_heat import EARTH, EARTH_ROCK, MID, _bulk_silicate_earth

IO_LIKE = tidal_heating.Forcing(1.9e27, 4.2e8, 0.05)
ICE_AGES = [1.0e-4, 2.0e-4, 0.03, 0.5, 1.0, 2.5, 4.0, 4.4, 4.567, 4.568, 6.0, 9.5]


def _from_formation(age, mass, radius, forcing, lid, surface_k):
    return ih.heat(radiogenic.heat_terms(_bulk_silicate_earth(), EARTH_ROCK.aluminium), MID, mass, radius, age, EARTH_ROCK, forcing, lid, surface_k)


def _taking_up(runs, age, mass, radius, forcing, lid, surface_k):
    return ih.heat(radiogenic.heat_terms(_bulk_silicate_earth(), EARTH_ROCK.aluminium), MID, mass, radius, age, EARTH_ROCK, forcing, lid, surface_k, runs=runs)


@pytest.mark.parametrize("mass,radius,forcing,lid,surface_k", [
    (EARTH, 6371.0, None, ih.STAGNANT, ih.SURFACE_TEMPERATURE_K),
    (EARTH, 6371.0, None, ih.MOBILE, 288.0),
    (0.1 * EARTH, 3000.0, None, ih.STAGNANT, 180.0),
    (3.0 * EARTH, 9000.0, None, ih.STAGNANT, 400.0),
    (EARTH, 6371.0, IO_LIKE, ih.STAGNANT, 288.0),
    (0.015 * EARTH, 1822.0, tidal_heating.Forcing(1.9e27, 4.2e8, 0.0041), ih.STAGNANT, 110.0),
], ids=["earth", "earth-mobile-lid", "small-cold", "super-earth-warm", "tidal-heated-earth", "small-tidal"])
def test_a_run_taken_up_at_another_age_gives_the_heat_a_run_from_the_formation_gives(mass, radius, forcing, lid, surface_k):
    """Ascending, descending and shuffled ages; every `Heat` equal to the last bit, the core's state included."""
    ages = sorted({*ICE_AGES, *(round(random.Random(3).random() * 9.0, 6) for _ in range(30))})
    for order in (ages, ages[::-1], random.Random(11).sample(ages, len(ages))):
        runs = {}
        for age in order:
            assert repr(_taking_up(runs, age, mass, radius, forcing, lid, surface_k)) == repr(_from_formation(age, mass, radius, forcing, lid, surface_k)), age


def test_an_age_takes_up_a_run_made_for_an_older_one_and_a_younger_one_and_does_not_make_the_shared_steps_again():
    """The run for 4.5 Gy holds the steps of 3 Gy up to the one that landed on it: the second age resumes from a kept state, past the start."""
    runs = {}
    _taking_up(runs, 4.5, EARTH, 6371.0, None, ih.STAGNANT, ih.SURFACE_TEMPERATURE_K)
    (run,) = runs.values()
    steps = len(run.times)
    found = run.resume(3.0)
    assert found is not None and 0 < found[0] < steps
    again = _taking_up(runs, 6.0, EARTH, 6371.0, None, ih.STAGNANT, ih.SURFACE_TEMPERATURE_K)
    assert repr(again) == repr(_from_formation(6.0, EARTH, 6371.0, None, ih.STAGNANT, ih.SURFACE_TEMPERATURE_K))
    assert len(run.times) > steps          # the run now reaches the older age


def test_the_heat_a_body_is_asked_for_at_many_ages_is_the_same_through_the_cache_that_remembers_its_runs():
    """`heat_flow` keeps the runs of each body: the answers at ascending ages are the ones a run from the formation gives."""
    for age in ICE_AGES:
        got = ih.heat_flow(_bulk_silicate_earth(), MID, EARTH, 6371.0, age, EARTH_ROCK, None, ih.STAGNANT, 288.0)
        assert repr(got) == repr(_from_formation(age, EARTH, 6371.0, None, ih.STAGNANT, 288.0))
