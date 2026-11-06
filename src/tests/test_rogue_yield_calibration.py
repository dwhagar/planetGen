# tests/test_rogue_yield_calibration.py

"""
TEST.499 (RP9): the rogue yield per star against the microlensing census and the terrestrial scale.

The planet architecture throws out bodies; those are the rogues of the stored donors. Two observed scales bound how many a star should leave behind:

* [R4254, R4255] (Sumi et al., the MOA-II 9-year survey): 21 free-floating planets per star over 0.33 Earth masses to Jupiter and beyond, with a lower 1-sigma of 8. The band is 8 to
  42 (the central value doubled; the survey's upper bound is looser). [R3191] limits Jupiter-mass planets, not the band's low end.
* [R0258] (Barclay et al.): 2.5 rocky free-floating bodies per star (embryos thrown out of a Mars- to Earth-mass range). The band is 0.8 to 7.5, a factor of three either way.

The observations are loose, so the bands are wide. The check counts what the architecture itself ejects, per star of the generator's own mass function, and what the Rogues levels let the
top-up bring a sector to. The donors alone fall far below both bands (about 0.1 body per star): the top-up, set by Boss's levels, makes up the rest. The test of the donor yield is a strict
expected failure, so it starts failing (and says so) the day the architecture's yield enters the bands; the finding is in
`/mnt/project-files/reports/rogue-yield-calibration-2026-11-06.md`.
"""

import random

import pytest

from planetgen import tuning
from planetgen.physics import constants, imf, planet_architecture as pa
from planetgen.util import draw

SUN = constants.SOLAR_MASS_TO_KG
EARTH = constants.EARTH_MASS_TO_KG

MICROLENSING_BAND = (8.0, 42.0)
"""tuple: Free-floating planets per star over 0.33 Earth masses and up: the lower 1-sigma of [R4254, R4255] to twice their central 21."""
MICROLENSING_FLOOR_EARTH = 0.33
TERRESTRIAL_BAND = (0.8, 7.5)
"""tuple: Rocky free-floating bodies per star, a factor of three either side of [R0258]'s 2.5."""
TERRESTRIAL_RANGE_EARTH = (0.1, 10.0)
STARS = 3000


def _donor_yields():
    """`(microlensing-range bodies, terrestrial bodies)` thrown out per star, over `STARS` stars of the generator's mass function and a disk metallicity spread."""
    rng = random.Random(11)
    micro = terrestrial = 0.0
    for i in range(STARS):
        mass = imf.sample(0.1, 2.0, rng)
        fe_h = rng.gauss(-0.1, 0.2)
        with draw.bound(20000 + i):
            arch = pa.build(mass * SUN, 2.7 * mass ** 2, 10.0 ** fe_h, fe_h, 100.0)
        for loss in pa.losses(arch):
            if not loss.ejected:
                continue
            each = loss.mass_kg / EARTH / loss.count
            ejected = loss.count * loss.ejected
            if each >= MICROLENSING_FLOOR_EARTH:
                micro += ejected
            if TERRESTRIAL_RANGE_EARTH[0] <= each < TERRESTRIAL_RANGE_EARTH[1]:
                terrestrial += ejected
    return micro / STARS, terrestrial / STARS


def test_the_donors_alone_leave_far_fewer_rogues_than_observed_and_are_not_zero():
    micro, terrestrial = _donor_yields()
    assert 0.0 < micro < 1.0 and 0.0 < terrestrial < 1.0
    assert micro < MICROLENSING_BAND[0] and terrestrial < TERRESTRIAL_BAND[0]


@pytest.mark.xfail(strict=True, reason="TEST.499 finding: the architecture ejects about 0.1 body per star, under both observed bands; the top-up makes up the difference")
def test_the_donor_yield_per_star_is_inside_the_observed_bands():
    micro, terrestrial = _donor_yields()
    assert MICROLENSING_BAND[0] <= micro <= MICROLENSING_BAND[1]
    assert TERRESTRIAL_BAND[0] <= terrestrial <= TERRESTRIAL_BAND[1]


def test_the_rogues_level_ceilings_reach_the_bands():
    ceilings = tuning.ROGUE_LEVEL_CEILINGS
    assert TERRESTRIAL_BAND[0] <= ceilings["Few"] <= TERRESTRIAL_BAND[1]           # Few is [R0258]'s 2.5
    assert TERRESTRIAL_BAND[0] <= ceilings["Normal"] <= TERRESTRIAL_BAND[1]
    assert MICROLENSING_BAND[0] <= ceilings["Many"] <= MICROLENSING_BAND[1]        # only Many reaches the microlensing census
    assert ceilings["None"] == 0.0
