# tests/test_scatterer_choice.py

"""
GEN.1456: the scatterer of a thrown or swept planet is the planet of the highest Safronov number among all those whose reach holds the body, not the neighbour the clearing
recorded alone (system-interactions.md: every body acts on every other). A body outside every other planet's reach, or a system with one planet, is scattered as before.
"""

import pytest

from planetgen.physics import planet_architecture as pa
from tests.test_architecture_losses import architecture_of, lost_of, seed

LIGHT = seed(2.0, 5.0, 5.0, False)


def _thrown(planets, pair):
    (loss,) = [loss for loss in pa.losses(architecture_of(planets, lost_of(thrown=(pair,)))) if loss.kind == "thrown"]
    return loss


def test_the_neighbour_the_clearing_recorded_scatters_when_it_is_the_strongest_in_reach():
    neighbour = seed(2.4, 60.0, 20.0, True)
    loss = _thrown([neighbour], (LIGHT, neighbour))
    assert loss.perturber_au == 2.4 and loss.perturber_mass_kg == neighbour.mass_kg


def test_a_heavier_giant_whose_reach_holds_the_body_scatters_it_instead():
    neighbour = seed(1.4, 20.0, 10.0, True)
    giant = seed(2.2, 318.0, 20.0, True)
    assert pa._reaches(giant, LIGHT.distance_au, 1.0), "the test needs the body inside the giant's reach"
    assert not pa._reaches(neighbour, LIGHT.distance_au, 1.0), "and the recorded neighbour (the clearing's wider spacing) outside its own"
    loss = _thrown([neighbour, giant], (LIGHT, neighbour))
    assert loss.perturber_mass_kg == giant.mass_kg and loss.ejected > _thrown([neighbour], (LIGHT, neighbour)).ejected


def test_a_giant_that_does_not_reach_the_body_changes_nothing():
    neighbour = seed(1.4, 20.0, 10.0, True)
    far = seed(40.0, 318.0, 20.0, True)
    assert not pa._reaches(far, LIGHT.distance_au, 1.0)
    assert _thrown([neighbour, far], (LIGHT, neighbour)) == _thrown([neighbour], (LIGHT, neighbour))


def test_a_tie_keeps_the_recorded_neighbour():
    neighbour = seed(2.1, 60.0, 20.0, True)
    twin = seed(2.1, 60.0, 20.0, True)
    assert pa._reaches(twin, LIGHT.distance_au, 1.0)
    assert _thrown([neighbour, twin], (LIGHT, neighbour)).perturber_au == neighbour.distance_au
