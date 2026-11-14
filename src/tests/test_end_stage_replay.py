# tests/test_end_stage_replay.py

"""
GEN.1547 (built ahead by Lane 07 for Lane 05): every draw of an end stage comes from a stream keyed to the body's birth, so remaking or moving a body never re-rolls it.

`draw.body_stream` keys a stream on the body's mass, the distance it formed at and its star's mass. The end stages (the air, the haze, the aerosols, the surface
water, the spin, the oddities, the mineral rows) still make streams whose keys hold the body's CURRENT distance, its star's age or its surface temperature, so a hair's
move re-rolls them (reports/build-order-air-once-2026-11-11.md, Round 3: 14 of the 65 bodies that changed class on the same dice).

* The helper tests pass today: same inputs same draws, a different body a different stream, and the seed string of `envelope-start` and `tidal-q-prime` is unchanged.
* `test_a_remade_body_makes_the_same_streams_when_it_moved_by_a_hair` replays `planets.remake_end_stages` on a body and on a copy moved by 0.05 percent in distance,
  0.1 percent in star age and 0.2 K in temperature, records the seed of every `draw.Stream` each makes, and requires the same set. It is `xfail(strict=True)`: today it
  names the offenders (refinish-rocky, haze, aerosols, body-strength, ore:*, isotopes:*, oddity:W19/W20 keyed on the current distance, age or temperature). Lane 05
  switches the call sites (`planets.py` lines 352, 524, 816; `aerosols.py` 449 and 662; `oddities.body_key`; the deposits; `isotope_record.py` 449;
  `giant_draw.py` 226; `kozai_bodies.py` 342) to `draw.body_stream_of` and removes the mark when the offender list is empty.
* `test_the_first_pass_air_does_not_move_the_run_stream` is the first-pass air, which draws from the run stream today (`planets.calculate_atmospheric_conditions`).
"""

import copy
import math
from collections import Counter

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.physics import envelope_capture, planets
from planetgen.util import draw
from tests import class_support

SOLAR_L_W = 3.828e26
EARTH_MASS_KG = 5.972e24


def test_a_body_stream_depends_on_the_birth_numbers_alone():
    a = draw.body_stream("haze", 5.97e24, 1.0, 1.989e30)
    b = draw.body_stream("haze", 5.97e24, 1.0, 1.989e30)
    assert [a.random() for _ in range(5)] == [b.random() for _ in range(5)]
    for other in (draw.body_stream("ore", 5.97e24, 1.0, 1.989e30), draw.body_stream("haze", 5.98e24, 1.0, 1.989e30), draw.body_stream("haze", 5.97e24, 1.1, 1.989e30),
                  draw.body_stream("haze", 5.97e24, 1.0, 2.0e30), draw.body_stream("haze", 5.97e24, 1.0, 1.989e30, slot="b")):
        assert other.random() != draw.body_stream("haze", 5.97e24, 1.0, 1.989e30).random()


def test_the_seed_string_is_the_one_envelope_start_and_tidal_q_prime_used():
    assert draw.body_stream("envelope-start", 5.97e24, 1.25, 1.989e30).random() == draw.Stream(f"envelope-start:{5.97e24!r}:{1.25!r}:{1.989e30!r}").random()


@pytest.fixture(scope="module")
def rocky():
    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    with draw.bound("replay-star"):
        star = Star(config)
    distance = math.sqrt(star.luminosity / SOLAR_L_W)
    with draw.bound("replay-planet"):
        return class_support.draw_planet(star, distance, 1, config, mass=EARTH_MASS_KG, moon_count=0)


def test_tidal_q_prime_is_unchanged_by_the_move_to_body_stream(rocky):
    stream = draw.Stream(f"tidal-q-prime:{float(rocky.mass)!r}:{float(planets.formation_distance_au(rocky))!r}:{float(rocky.star.mass)!r}")
    low, high = envelope_capture.super_puff.TIDAL_Q_PRIME_RANGE
    with draw.bound(stream):
        expected = math.exp(math.log(low) + draw.random() * (math.log(high) - math.log(low)))
    assert envelope_capture.tidal_q_prime(rocky) == expected


def _streams_made(planet):
    """The seed strings of every `draw.Stream` that `remake_end_stages` makes for `planet` (the integer seeds follow from them and are left out)."""
    made = []
    real = draw.Stream.__init__

    def recording(self, seed):
        made.append(seed)
        real(self, seed)

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(draw.Stream, "__init__", recording)
        with draw.bound("replay-outer"):
            planets.remake_end_stages(planet, planet)
    return Counter(s for s in made if isinstance(s, str))


@pytest.mark.slow
def test_a_remade_body_makes_the_same_streams_when_it_moved_by_a_hair(rocky):
    still = copy.deepcopy(rocky)
    moved = copy.deepcopy(rocky)
    moved.distance *= 1.0005
    moved.star.age *= 1.001
    moved.surface_temperature += 0.2
    before, after = _streams_made(still), _streams_made(moved)
    # `abundance` is the star's makeup (keyed on its mass and galactic distance, never its age); it is made once per age in a cache, so only the moved copy makes the stream again
    # `ice-giant-heat` is keyed on the body's mass alone (`ice_giant_heat.SEED_DIGITS`) and made once per mass in the same cache, so only the first copy makes it
    offenders = sorted({s.split(":")[0] for s in (before - after) | (after - before)} - {"abundance", "ice-giant-heat"})
    assert not offenders, f"streams keyed on the body's present state: {offenders}"


@pytest.mark.slow
def test_the_first_pass_air_does_not_move_the_run_stream(rocky):
    planet = copy.deepcopy(rocky)
    run = draw.Stream("replay-run")
    with draw.bound(run):
        before = draw.getstate()
        planets.calculate_atmospheric_conditions(planet)
        after = draw.getstate()
    assert before == after
