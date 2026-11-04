# tests/test_ejection_events.py

"""GEN.1111 part 2: the events a rogue is built from, drawn from a donor's record."""

import math
import random

import pytest

from planetgen.generation import ejection_events
from planetgen.physics import constants, rogue_kick
from planetgen.physics.planet_architecture import Loss
from tests.test_gen_builder_internals import make_system

EARTH_KG = constants.EARTH_MASS_TO_KG


class Stream:
    """A seeded stand-in for `draw` (`random` only)."""

    def __init__(self, seed):
        self._random = random.Random(seed)

    def random(self):
        return self._random.random()


def _loss(ejected=1.0, count=1, **over):
    fields = dict(kind="thrown", count=count, mass_kg=count * 100.0 * EARTH_KG, core_kg=count * 10.0 * EARTH_KG, distance_au=5.0, formation_au=5.0,
                  time_yr=1e7, perturber_au=3.0, perturber_mass_kg=3000.0 * EARTH_KG, theta=8.0, hill_radius_au=0.35, ice_rich=True,
                  ejected=ejected, collided=0.0, bound=1.0 - ejected)
    fields.update(over)
    return Loss(**fields)


@pytest.fixture(scope="module")
def donor():
    return make_system().star


def test_a_loss_that_always_ends_ejected_gives_an_event_with_every_field_set(donor):
    [event] = ejection_events.events("uid-1", donor, [_loss(1.0)], Stream(1))
    assert event.donor_uid == "uid-1" and event.index == 0
    assert event.mass_kg == 100.0 * EARTH_KG and event.core_kg == 10.0 * EARTH_KG and event.ice_rich
    assert 0.0 < event.delay_gy <= donor.age and event.age_gy == donor.age
    assert event.since_ejection_gy == pytest.approx(donor.age - event.delay_gy)
    assert len(event.kick_km_s) == 3 and len(event.orbital_velocity_km_s) == 3 and event.moon_limit_km > 0.0


def test_a_loss_that_never_ends_ejected_gives_none_but_still_uses_an_index(donor):
    events = ejection_events.events("uid-2", donor, [_loss(0.0), _loss(1.0)], Stream(2))
    assert [e.index for e in events] == [1], "the index counts the bodies of every loss, ejected or not"


def test_the_events_of_a_donor_are_a_pure_function_of_its_stream(donor):
    again = [ejection_events.events("u", donor, [_loss(0.6), _loss(0.6)], Stream(7)) for _ in range(2)]
    assert again[0] == again[1]
    assert ejection_events.events("u", donor, [_loss(0.6), _loss(0.6)], Stream(8)) != again[0]


def test_a_share_of_one_half_ejects_about_half(donor):
    stream = Stream(11)
    kept = sum(len(ejection_events.events("u", donor, [_loss(0.5)], stream)) for _ in range(2000))
    assert 900 < kept < 1100


def test_the_orbital_velocity_is_the_circular_speed_at_the_distance_and_the_kick_follows_the_scatter_speed(donor):
    circular = rogue_kick.circular_speed_km_s(donor.mass, 5.0)
    speeds, orbits = [], []
    stream = Stream(5)
    for _ in range(300):
        [event] = ejection_events.events("u", donor, [_loss(1.0)], stream)
        orbits.append(math.sqrt(sum(c * c for c in event.orbital_velocity_km_s)))
        speeds.append(math.sqrt(sum(c * c for c in event.kick_km_s)))
    assert all(o == pytest.approx(circular) for o in orbits)
    scale = rogue_kick.scatter_speed_km_s(3000.0 * EARTH_KG, 3.0, 100.0 * EARTH_KG)
    assert 0.2 * scale < sorted(speeds)[len(speeds) // 2] < 5.0 * scale, "the median kick is of the scatter speed's order"


def test_directions_are_isotropic_so_no_axis_is_favoured(donor):
    stream = Stream(21)
    sums = [0.0, 0.0, 0.0]
    for _ in range(1500):
        [event] = ejection_events.events("u", donor, [_loss(1.0)], stream)
        for axis in range(3):
            sums[axis] += event.orbital_velocity_km_s[axis]
    circular = rogue_kick.circular_speed_km_s(donor.mass, 5.0)
    for total in sums:
        assert abs(total) < 4.0 * circular * math.sqrt(1500 / 3.0), "the sum of an isotropic set is within four standard deviations of zero"
