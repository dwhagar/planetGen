# tests/test_cycle_dial.py

"""GEN.1720: the carbon-cycle dial of a Class M abundance tier (`abundance_tier.cycle_share`, `continental_weathering.draw_cycle`): Science draws as it always did, a tier draws the two constants from the part of their
published ranges nearest the Earth's balance, and a world there holds its CO2 where the life classes ask for it."""

import random
from types import SimpleNamespace

import pytest

from planetgen import tuning
from planetgen.generation import abundance_tier
from planetgen.generation.config import SystemConfig
from planetgen.physics import continental_weathering, surface_air, surface_pressure, waterworld_carbon
from planetgen.util import draw

Q_EARTH = surface_pressure.HEAT_FLUX_EARTH_W_M2


def _body(tier, zone="e", body_type="t", is_moon=False):
    config = SystemConfig()
    config.ABUNDANCE_TIER = tier
    distance = {"e": 1.0, "h": 0.3, "c": 3.0, None: 0.3}[zone]       # outside the habitable zone 0.95 to 1.37 AU widened by half its width: 0.74 to 1.58 AU (GEN.1807)
    return SimpleNamespace(system_config=config, zone=zone, body_type=body_type, is_moon=is_moon, habitable_zone=(0.95, 1.37), distance=distance, star=None)


def test_science_makes_no_draw_and_the_constants_are_those_the_law_always_drew():
    assert abundance_tier.cycle_share(_body(None)) == 0.0 == abundance_tier.cycle_share(_body("science"))
    with draw.bound(draw.Stream(7)):
        design, expected = {}, {}
        got = continental_weathering.draw_cycle(design)
        after = draw.getstate()
    with draw.bound(draw.Stream(7)):
        want = (continental_weathering.bracket(expected, "weathering_lambda_scale", continental_weathering.LAMBDA_SCALE),
                continental_weathering.bracket(expected, "outgassing_exponent", surface_pressure.OUTGASSING_EXPONENT))
        assert draw.getstate() == after
    assert got == want and "cycle_nudged" not in design


def test_a_tier_moves_the_ambient_stream_exactly_as_far_as_science_does():
    """The nudge is rolled on a stream of its own: after the two constants the ambient stream stands where Science's draw leaves it, so the rest of the planet (the air recipe and the other design constants) is the same planet."""
    with draw.bound(draw.Stream(11)):
        continental_weathering.draw_cycle({})
        science = draw.getstate()
    for share in (0.5, 1.0):
        with draw.bound(draw.Stream(11)):
            continental_weathering.draw_cycle({}, share)
            assert draw.getstate() == science, share


def test_a_tier_draws_both_constants_inside_the_nudged_ranges_and_the_published_ones():
    low_m, high_m = tuning.ABUNDANCE_CYCLE_OUTGASSING_RANGE[:2]
    low_l, high_l = tuning.ABUNDANCE_CYCLE_LAMBDA_RANGE[:2]
    assert surface_pressure.OUTGASSING_EXPONENT[0] <= low_m < high_m <= surface_pressure.OUTGASSING_EXPONENT[1]       # inside the source range
    assert continental_weathering.LAMBDA_SCALE[0] <= low_l < high_l <= continental_weathering.LAMBDA_SCALE[1]
    with draw.bound(draw.Stream(3)):
        for _ in range(40):
            design = {}
            scale, exponent = continental_weathering.draw_cycle(design, 1.0)
            assert design["cycle_nudged"] is True
            assert low_l <= scale <= high_l and low_m <= exponent <= high_m
        shares = [continental_weathering.draw_cycle({}, 0.5) for _ in range(400)]
    nudged = sum(low_m <= exponent <= high_m and low_l <= scale <= high_l for scale, exponent in shares)
    assert 150 < nudged < 330        # half of them, the rest over the whole ranges (a few of those fall inside too)


def test_the_dial_reaches_rocky_planets_in_the_ecosphere_only():
    assert abundance_tier.cycle_share(_body("abundant")) == tuning.ABUNDANCE_TIER_DIALS["abundant"][5] > 0.0
    for changes in ({"zone": "h"}, {"zone": None}, {"body_type": "j"}, {"is_moon": True}):
        assert abundance_tier.cycle_share(_body("abundant", **changes)) == 0.0


def test_a_setting_built_by_hand_is_outside_every_tier():
    field = {name: field for name, field in surface_air.Setting.__dataclass_fields__.items()}["cycle_share"]
    assert field.default == 0.0


@pytest.mark.parametrize("heat_over_earth", [0.6, 1.0, 1.5, 2.0])
def test_a_planet_drawn_from_the_nudged_ranges_holds_its_co2_between_the_photosynthesis_floor_and_the_class_limit(heat_over_earth):
    """At 290 K and the Earth's absorbed starlight, a world whose heat flow is 0.6 to 2 times the Earth's holds 2 Pa to 1 kPa of CO2 for every corner of the nudged ranges: the band the life classes ask
    for (the floor of photosynthesis and the class limit); over the whole published ranges the same worlds spread from under a pascal to over 100 kPa (reports/gen1720-co2-attribution-2026-11-14.md)."""
    law = waterworld_carbon.SeafloorLaw(75.0e3, 0.27, 0.47)
    low_m, high_m = tuning.ABUNDANCE_CYCLE_OUTGASSING_RANGE[:2]
    low_l, high_l = tuning.ABUNDANCE_CYCLE_LAMBDA_RANGE[:2]
    for scale in (low_l, high_l):
        for exponent in (low_m, high_m):
            pascal = 1.0e5 * continental_weathering.pco2_bar(200.0, Q_EARTH * heat_over_earth, 0.29, 290.0, 238.0, (scale, exponent), law)
            assert 2.0 <= pascal <= 1.0e3, (scale, exponent, pascal)
