# tests/test_radius_valley_sweep.py

"""
TEST.676: the radius valley of the generated sample, and the rogue planets that share its envelopes (GEN.1323).

`test_radius_valley.py` tunes the survey that GEN.822 fitted and `test_envelope_capture.py` tests the pieces. This file measures what the generator
makes: a sweep of G2V systems, the share of the planets of 1 to 4 Earth radii that are above 2.0 inside 100 days ([R1478]: about half), the
planets beyond 100 days, the largest a rocky planet comes out, and the rogue planets, which have no star to strip an envelope.
"""

import collections

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.generation.phenomena import starless_host
from planetgen.generation.star import Star
from planetgen.physics import constants, envelope_capture, planet_classes, planet_type, super_puff
from planetgen.util import draw

EARTH = constants.EARTH_MASS_TO_KG
SYSTEMS = 60
SEED = 20261107
INNER_DAYS = 100.0
"""float: [R1478] count the planets inside 100 days."""

Row = collections.namedtuple("Row", "days radius mass envelope density planet_type planet_class")


@pytest.fixture(scope="module")
def sweep():
    """One row per rocky planet of `SYSTEMS` G2V systems."""
    draw.set_run_seed(SEED)
    rows = []
    for _ in range(SYSTEMS):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        for planet in StarSystem(system_config=cfg).all_planets:
            if getattr(planet, "body_type", None) != "t" or planet.is_moon:
                continue
            rows.append(Row(planet.period * 365.25, planet.radius / constants.EARTH_RADIUS_KM, planet.mass / EARTH,
                            planet.envelope_mass_fraction or 0.0, planet.density, planet_type.planet_type(planet), planet.planet_class))
    return rows


def test_about_half_of_the_inner_planets_of_one_to_four_earth_radii_are_above_two(sweep):
    inner = [row for row in sweep if row.days < INNER_DAYS and 1.0 <= row.radius <= 4.0]
    assert len(inner) >= 25
    share = sum(row.radius > 2.0 for row in inner) / len(inner)
    assert 0.3 <= share <= 0.7, (len(inner), share)


def test_the_valley_is_populated_from_both_sides_and_the_peaks_are_not_empty(sweep):
    """The three bins of the valley (below 1.5, the gap 1.5 to 2.0, and above 2.0) each hold planets of the inner sample; a sample with an empty side has lost
    one of the two populations the valley divides."""
    inner = [row.radius for row in sweep if row.days < INNER_DAYS and 1.0 <= row.radius <= 4.0]
    bins = (sum(r < 1.5 for r in inner), sum(1.5 <= r <= 2.0 for r in inner), sum(r > 2.0 for r in inner))
    assert min(bins) >= 5, bins


def test_a_planet_beyond_a_hundred_days_of_four_earth_masses_or_more_is_mostly_a_sub_neptune(sweep):
    far = [row for row in sweep if row.days >= INNER_DAYS and 4.0 <= row.mass <= 20.0]
    assert len(far) >= 15
    envelope = sum(row.envelope >= planet_type.TERRESTRIAL_MAX_ENVELOPE for row in far) / len(far)
    assert envelope >= 0.6, (len(far), envelope)


def test_a_rocky_planet_is_never_larger_than_a_sub_neptune_can_be(sweep):
    """The start is held at 0.2 of the mass (EC5): the radius law of Lopez and Fortney is fitted to that, and no core of a few Earth masses is a puffed-up Saturn
    (9 Earth radii and more). A cold 4 Earth-mass core under the full 0.2 reaches 6.5 and a cold 3.5 Earth-mass one 7.05 (a lighter core puffs more; the super-puffs of Kepler-51 are the extreme at 6 to 9), so the bound is 7.5."""
    # GEN.1343: a super-puff is the extreme of Kepler-51, 6 to 9 Earth radii at under 0.1 g/cc. GEN.1406: the planet type names one only for 10 to 100 days and under 10 Earth masses
    # (`super_puff.is_super_puff`), while the puffed start (`draw_puff_start`) and the obliquity tide of GEN.1347 (which inflates a planet of 20 Earth fluxes and more by up to a factor of 2)
    # make a puff of a few Earth masses at 7 days or of 10.2 Earth masses at 26 days too: with the tide switched off none of the 60 systems has a planet over 6.5 Earth radii, with it
    # two planets of 8 and 9.9 Earth radii at 0.04 and 0.06 g/cc come out. A body under the definition's density edge is a puff whatever its type.
    puffed = [row for row in sweep if row.planet_type == planet_type.SUPER_PUFF or row.density < super_puff.DENSITY_RANGE_G_CM3[1]]
    ordinary = [row for row in sweep if row not in puffed]
    assert max(row.radius for row in ordinary) <= 7.5
    assert max((row.radius for row in puffed), default=0.0) <= 2.0 * 7.5       # a puff of 7.5 Earth radii under the tide's largest factor of 2
    assert max(row.density for row in ordinary if row.envelope >= planet_type.TERRESTRIAL_MAX_ENVELOPE) < 6.0


def test_no_puffed_planet_escapes_the_super_puff_type(sweep):
    """GEN.1470: a transitional planet of up to 30 Earth masses under 0.1 g/cc is typed a super-puff whatever its period (the type is mass and density)."""
    puffed = [row for row in sweep if row.envelope >= planet_type.TERRESTRIAL_MAX_ENVELOPE and row.density < 0.1 and row.mass <= 30.0]
    assert all(row.planet_type == planet_type.SUPER_PUFF for row in puffed), puffed


def test_the_super_puff_share_of_the_inner_sub_neptunes_is_counted_inside_a_hundred_days(sweep):
    """Row PF3 is a share of the 2 to 10 Earth-mass sub-Neptunes inside 100 days; with the type free of the period the count is made inside 100 days here, not by the type."""
    inner = [row for row in sweep if row.days <= INNER_DAYS and 2.0 <= row.mass <= 10.0 and row.envelope >= planet_type.TERRESTRIAL_MAX_ENVELOPE]
    assert len(inner) >= 10
    share = sum(row.planet_type == planet_type.SUPER_PUFF for row in inner) / len(inner)
    assert share <= 4.0 * super_puff.START_SHARE, (len(inner), share)        # one planet in a sweep this small is a tenth: a guard against the type sweeping up the sub-Neptunes, not a measurement of the row


def test_a_planet_with_an_envelope_is_a_gas_dwarf_or_an_ice_giant_not_a_rocky_letter(sweep):
    deep = [row for row in sweep if row.envelope >= planet_type.TERRESTRIAL_MAX_ENVELOPE]
    assert len(deep) >= 40
    letters = collections.Counter(row.planet_class for row in deep)
    assert sum(count for letter, count in letters.items() if letter in ("T", "I")) / len(deep) >= 0.9, letters


def test_the_light_planets_stay_rocky(sweep):
    light = [row for row in sweep if 1.0 <= row.mass < 1.5]
    assert len(light) >= 10
    assert sum(row.envelope < 0.01 for row in light) / len(light) >= 0.6


# --- rogue planets: the same constructor, no star, no loss ------------------------------------------------------------------------------------------

@pytest.fixture(scope="module")
def host_star():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    with draw.bound("valley-sweep:star"):
        star = Star(cfg)
    star.age = 4.5
    return star, cfg


@pytest.mark.parametrize("mass_earth", [5.0, 8.0, 12.0])
def test_a_rogue_planet_keeps_the_whole_envelope_it_started_with(host_star, mass_earth):
    star, cfg = host_star
    kept = []
    for seed in range(8):
        with draw.bound(f"valley-sweep:rogue:{mass_earth}:{seed}"):
            rogue = starless_host.build_ejected_planet(cfg, star, 4.5, mass_earth * EARTH, formation_au=0.3 + 0.08 * seed, moon_count=0)
        assert rogue.body_type == "t"
        if rogue._captured_envelope:
            kept.append(rogue)
            start = rogue._captured_envelope
            assert rogue.envelope_mass_fraction == pytest.approx(start)      # no light, no loss: what it started with
            assert 0.0 < start <= envelope_capture.MAX_START + 1e-9
            assert rogue.radius > rogue._bare_core["radius"] * 1.15
            assert planet_type.planet_type(rogue) in (planet_type.TRANSITIONAL, planet_type.TERRESTRIAL)
            assert planet_classes.classify(rogue) == rogue.planet_class
            if rogue.envelope_mass_fraction >= planet_type.TERRESTRIAL_MAX_ENVELOPE:
                assert rogue.planet_class in ("T", "I"), rogue.planet_class
    assert len(kept) == 8, len(kept)       # the start is keyed on the planet's mass, formation distance and star: each formation distance is another start


def test_a_rogue_planet_outside_the_mass_range_has_no_envelope(host_star):
    star, cfg = host_star
    for mass_earth in (0.4, 0.8):
        with draw.bound(f"valley-sweep:rogue-light:{mass_earth}"):
            rogue = starless_host.build_ejected_planet(cfg, star, 4.5, mass_earth * EARTH, formation_au=2.0, moon_count=0)
        assert not rogue._captured_envelope
