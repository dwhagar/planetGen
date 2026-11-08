# tests/test_super_puff_branch.py

"""
GEN.1343 (W15): a small share of the sub-Neptunes inside 100 days start as super-puffs, a thick envelope of 7 to 45 percent
(`envelope_capture.draw_puff_start`), the XUV and the age term of the law fade it (`envelope_capture.kept`, the radius of
`super_puff.radius_earth`), and a planet that comes out under 0.1 g/cc in the period range of the Kepler super-puffs is typed a super-puff
(`planet_type`). The module's own laws are tested in test_super_puff.py, the sub-Neptune envelope in test_envelope_capture.py.
"""

import math
from types import SimpleNamespace

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.physics import constants, envelope_capture, planet_type, super_puff
from planetgen.util import draw

EARTH = constants.EARTH_MASS_TO_KG


@pytest.fixture()
def young_sun():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    with draw.bound("super-puff-branch:star"):
        star = Star(cfg)
    star.age = 0.5
    return star


def _body(mass_earth, distance_au, star):
    return SimpleNamespace(star=star, mass=mass_earth * EARTH, distance=distance_au, eccentricity=0.0, is_moon=False)


def test_the_period_follows_keplers_third_law(young_sun):
    body = _body(5.0, 1.0, young_sun)
    expected = 365.25 * math.sqrt(1.0 / (young_sun.mass / constants.SOLAR_MASS_TO_KG))
    assert envelope_capture.period_days(body) == pytest.approx(expected, rel=1e-3)


def test_only_a_sub_neptune_inside_the_period_range_can_start_as_a_super_puff():
    assert envelope_capture.puff_eligible(5.0, 50.0)
    assert not envelope_capture.puff_eligible(1.5, 50.0)           # lighter than Lee and Chiang's class
    assert not envelope_capture.puff_eligible(12.0, 50.0)          # heavier than the module's top
    assert not envelope_capture.puff_eligible(5.0, 400.0)          # outside the period range
    with draw.bound("super-puff-branch:ineligible"):
        assert envelope_capture.draw_puff_start(1.5, 50.0) == 0.0 and envelope_capture.draw_puff_start(5.0, 400.0) == 0.0


def test_the_share_that_starts_as_a_super_puff_is_the_occurrence_row_and_its_envelope_is_in_the_range():
    with draw.bound("super-puff-branch:draws"):
        starts = [envelope_capture.draw_puff_start(5.0, 50.0) for _ in range(20000)]
    puffs = [value for value in starts if value > 0.0]
    assert len(puffs) / len(starts) == pytest.approx(super_puff.START_SHARE, rel=0.2)
    low, high = super_puff.PUFF_ENVELOPE
    assert all(low <= value <= high for value in puffs)
    assert max(puffs) > 0.3 > min(puffs) < 0.15         # it spans the range, not one edge


def test_a_cold_young_super_puff_keeps_its_envelope_and_a_hot_one_loses_it(young_sun):
    cold = envelope_capture.kept(_body(5.0, 0.4, young_sun), 0.41)       # about 90 days, 6 Earth fluxes
    hot = envelope_capture.kept(_body(5.0, 0.05, young_sun), 0.41)
    assert cold is not None
    fraction, radius = cold
    assert fraction > 0.3 and radius > 8.0
    assert super_puff.density_g_cm3(5.0, radius) < super_puff.DENSITY_RANGE_G_CM3[1]
    assert hot is not None and hot[1] > radius          # the XUV loss starts at 1 Gyr: at 0.5 Gyr the hotter planet is only the more inflated
    young_sun.age = 4.5
    hot_old, cold_old = envelope_capture.kept(_body(5.0, 0.05, young_sun), 0.41), envelope_capture.kept(_body(5.0, 0.4, young_sun), 0.41)
    assert hot_old is None or hot_old[1] < 0.5 * cold_old[1]
    assert cold_old is not None and cold_old[0] > 0.3


def test_the_age_term_contracts_a_super_puff(young_sun):
    body = _body(5.0, 0.4, young_sun)
    _, young = envelope_capture.kept(body, 0.41)
    young_sun.age = 4.5
    _, old = envelope_capture.kept(body, 0.41)
    assert old < young


def test_a_planet_that_passes_the_gate_is_typed_a_super_puff(young_sun):
    body = _body(5.0, 0.4, young_sun)
    body.envelope_mass_fraction, body.radius = 0.4, 9.0 * constants.EARTH_RADIUS_KM
    body.body_type = "t"
    assert planet_type.planet_type(body) == planet_type.SUPER_PUFF
    body.radius = 3.0 * constants.EARTH_RADIUS_KM                       # a dense sub-Neptune
    assert planet_type.planet_type(body) == planet_type.TRANSITIONAL
    body.radius, young_sun.age = 9.0 * constants.EARTH_RADIUS_KM, 8.0
    assert planet_type.planet_type(body) == planet_type.SUPER_PUFF      # no age in the gate (GEN.1347): the obliquity tide keeps an old one inflated
    assert planet_type.planet_type({"envelope_mass_fraction": 0.4, "body_type": "t"}) == planet_type.TRANSITIONAL    # a stored row has no star to judge by


def _tilted(mass_earth, distance_au, star, tilt):
    body = _body(mass_earth, distance_au, star)
    body.axial_tilt_deg = tilt
    return body


def test_the_obliquity_tide_inflates_a_hot_tilted_planet_and_not_an_upright_or_untilted_one(young_sun):
    plain = envelope_capture.kept(_body(5.0, 0.08, young_sun), 0.41)
    assert plain is not None and plain[0] > 0.2
    assert envelope_capture.kept(_body(5.0, 0.08, young_sun), 0.41, tides=True)[1] == pytest.approx(plain[1])        # no tilt yet
    assert envelope_capture.kept(_tilted(5.0, 0.08, young_sun, 0.0), 0.41, tides=True)[1] == pytest.approx(plain[1])    # upright
    tilted = envelope_capture.kept(_tilted(5.0, 0.08, young_sun, 60.0), 0.41, tides=True)
    assert tilted[0] == pytest.approx(plain[0])
    assert 1.02 * plain[1] < tilted[1] <= super_puff.TIDAL_MAX_RATIO * plain[1] * (1 + 1e-9)


def test_a_cold_planet_has_no_obliquity_tide_whatever_its_tilt(young_sun):
    cold = _tilted(5.0, 0.8, young_sun, 60.0)          # about 1.5 Earth fluxes
    assert envelope_capture.kept(cold, 0.41, tides=True)[1] == pytest.approx(envelope_capture.kept(cold, 0.41)[1])


def test_the_tidal_quality_is_the_planets_own_and_in_the_range(young_sun):
    body = _body(5.0, 0.12, young_sun)
    value = envelope_capture.tidal_q_prime(body)
    assert super_puff.TIDAL_Q_PRIME_RANGE[0] <= value <= super_puff.TIDAL_Q_PRIME_RANGE[1]
    assert envelope_capture.tidal_q_prime(body) == value
    with draw.bound("super-puff-branch:q"):
        first = draw.random()
        envelope_capture.tidal_q_prime(body)
    with draw.bound("super-puff-branch:q"):
        assert draw.random() == first          # the run's stream does not move
