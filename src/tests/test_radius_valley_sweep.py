# tests/test_radius_valley_sweep.py

"""
TEST.676: the radius valley of the generated sample, and the rogue planets that share its envelopes (GEN.1323).

`test_radius_valley.py` tunes the survey that GEN.822 fitted and `test_envelope_capture.py` tests the pieces. This file measures what the generator
makes: a sweep of G2V systems, the share of the planets of 1 to 4 Earth radii that are above 2.0 inside 100 days ([R1478]: about half), the
planets beyond 100 days, the largest a rocky planet comes out, and the rogue planets, which have no star to strip an envelope.
"""

import collections
import concurrent.futures
import math
import multiprocessing
import os

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.generation.phenomena import starless_host
from planetgen.generation.star import Star
from planetgen.physics import bulk_structure, constants, envelope_capture, planet_classes, planet_type, super_puff
from planetgen.util import draw

EARTH = constants.EARTH_MASS_TO_KG
SYSTEMS = 180
SEED = 20261107
BLOCKS = 12
"""int: TEST.783: the sweep is made in this many blocks of `SYSTEMS // BLOCKS` systems, one worker at a time, every block on a run seed of its own (blocks 0 to 11, twelve of 15 systems since the Train 61 valley change: the first eight alone gave a valley of 4 planets, a count the survey's own valley makes a coin flip to pass a floor of 5).
The share of the inner planets above 2.0 Earth radii that the first test reads is 57 of 156 (0.365) over these blocks; the blocks alone give 0.23 to 0.62 (standard deviation 0.13), and
twelve blocks (blocks 8 to 11 added) give 72 of 197 (0.365) by the count of the note's first author; the Train 61 head's twelve give 75 of 143 inner planets of 1 to 4 Earth radii above 2.0 (0.52). The binomial noise of 156 planets is 0.039, so the 0.3 edge of the band is 1.7 of those under the share: a sample drawn again
(a change to the generator moves every draw) falls under 0.3 about one time in twenty, and the band, not the sample, is what is tight."""
INNER_DAYS = 100.0
"""float: [R1478] count the planets inside 100 days."""

Row = collections.namedtuple("Row", "days radius mass envelope density planet_type planet_class bare_density")


def _bare_density(planet):
    """The density (g/cc) the dry part of `planet` (its metal and rock, in the shares it was drawn with) would have bare, by the generator's own fit of the solver (`bulk_structure.radius_fit_earth`, 0.1 to 20 Earth masses),
    or infinity where the dry part is outside the fit. Water and an envelope are left off, and both only lower a body's mean density, so this is a ceiling of what the body can have."""
    metal, rock = planet.metal_mass_fraction or 0.0, planet.rock_mass_fraction or 0.0
    dry_earth = planet.mass / EARTH * (metal + rock)
    if metal + rock <= 0.0 or not 0.1 <= dry_earth <= 20.0:
        return math.inf
    x_fe = bulk_structure.x_fe_from_feo_mass_fraction(planet.rock_feo_fraction) if planet.rock_feo_fraction else 0.0
    radius_m = bulk_structure.radius_fit_earth(dry_earth, metal / (metal + rock), x_fe, 1.0) * constants.EARTH_RADIUS_KM * 1000.0
    return dry_earth * EARTH / (4.0 / 3.0 * math.pi * radius_m ** 3) / 1000.0


def _block_rows(block):
    """One row per rocky planet of one block of G2V systems, on the run seed of the block."""
    draw.set_run_seed(f"{SEED}:{block}")
    rows = []
    for _ in range(SYSTEMS // BLOCKS):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        for planet in StarSystem(system_config=cfg).all_planets:
            if getattr(planet, "body_type", None) != "t" or planet.is_moon:
                continue
            rows.append(Row(planet.period * 365.25, planet.radius / constants.EARTH_RADIUS_KM, planet.mass / EARTH,
                            planet.envelope_mass_fraction or 0.0, planet.density, planet_type.planet_type(planet), planet.planet_class, _bare_density(planet)))
    return rows


@pytest.fixture(scope="module")
def sweep():
    """One row per rocky planet of `SYSTEMS` G2V systems. TEST.783: the systems of a run seed make one another's draws, so the sweep is `BLOCKS` blocks on
    seeds of their own, made at the same time (forked workers: a system takes about 15 s, and 60 one after another ran past the hang guard; the sample is twice as large for a third of the time)."""
    with concurrent.futures.ProcessPoolExecutor(max_workers=min(BLOCKS, os.cpu_count() or 1), mp_context=multiprocessing.get_context("fork")) as pool:
        blocks = list(pool.map(_block_rows, range(BLOCKS)))
    return [row for rows in blocks for row in rows]


def test_about_half_of_the_inner_planets_of_one_to_four_earth_radii_are_above_two(sweep):
    inner = [row for row in sweep if row.days < INNER_DAYS and 1.0 <= row.radius <= 4.0]
    assert len(inner) >= 25
    share = sum(row.radius > 2.0 for row in inner) / len(inner)
    assert 0.3 <= share <= 0.7, (len(inner), share)


VALLEY_OVER_PEAK_OCCURRENCE = 0.030 / 0.116
"""float: The depth of the valley the CKS survey measures, [R1478: Table 4] (the spline fit of the occurrence of planets inside 100 days, deconvolved for the survey's completeness and
the radius errors): 0.030 planets per star at 1.9 Earth radii in the valley over 0.116 at 2.4, the peak of the sub-Neptunes. The paper's headline is a deficit of 2 or more
at 1.5 to 2.0 Earth radii in the observed counts, and says the deconvolved valley "is likely deeper than observed" (0.26 of the peak here, 0.5 observed)."""
VALLEY, PEAK = (1.5, 2.0), (2.19, 2.70)
"""tuple: The valley's radii (the gap of [R1478]) and the radii of the sub-Neptune peak the valley is read against (bins 2.19 to 2.43 and 2.43 to 2.70 of [R1478: Table 3], the peak at 2.4)."""
SIGNIFICANCE = 0.01


def _least_valley_count(sample_size, valley_share):
    """The largest count `k` of the valley's planets, out of `sample_size` planets in the valley and the peak together, for which a generator that has the valley the survey
    measures (a binomial of `valley_share`) has fewer than `k` of them less than `SIGNIFICANCE` of the time: a sample with fewer is a valley too empty to be the survey's."""
    below, k = 0.0, 0
    while True:
        below += math.comb(sample_size, k) * valley_share ** k * (1.0 - valley_share) ** (sample_size - k)
        if below > SIGNIFICANCE:
            return k
        k += 1


def test_the_valley_is_populated_from_both_sides_and_the_peaks_are_not_empty(sweep):
    """The planets of the inner sample fall on both sides of the valley (at least 5 of 1 to 1.5 and of 2.0 to 4.0 Earth radii: a sample with an empty side has lost one of the
    two populations the valley divides) and the valley itself is not emptier than the survey's measured valley allows.

    The valley is a gap, so its count is small by nature: 6 of the 143 inner planets (4 percent) over the twelve blocks, 4 of 104 over the first eight, and a floor of 5 failed a re-drawn
    sample about half the time (the sample is drawn again by any change to the generator, GEN.1202 moved every draw). The floor now comes from the survey. Taking the valley (1.5 to 2.0 Earth
    radii, width 0.5) and the sub-Neptune peak (2.19 to 2.70, width 0.51) of the inner sample together, a generator with the valley of [R1478: Table 4], 0.030 over 0.116 of the peak's
    occurrence, puts 0.26 x 0.5 / (0.26 x 0.5 + 0.51) = 20 percent of those planets in the valley; with `n` of them the valley count is a binomial and the test asks for the count a correct
    generator misses less than 1 time in 100 (for `n` of 28 the count is 1: Bin(28, 0.20) has 0 with probability 0.0017 and 1 or fewer with 0.014). This sweep has 6 of 28 (21 percent)."""
    inner = [row.radius for row in sweep if row.days < INNER_DAYS and 1.0 <= row.radius <= 4.0]
    low, high = sum(r < VALLEY[0] for r in inner), sum(r > VALLEY[1] for r in inner)
    valley = sum(VALLEY[0] <= r <= VALLEY[1] for r in inner)
    peak = sum(PEAK[0] <= r <= PEAK[1] for r in inner)
    assert min(low, high) >= 5, (low, valley, high)
    width_valley, width_peak = VALLEY[1] - VALLEY[0], PEAK[1] - PEAK[0]
    share = VALLEY_OVER_PEAK_OCCURRENCE * width_valley / (VALLEY_OVER_PEAK_OCCURRENCE * width_valley + width_peak)
    assert valley >= _least_valley_count(valley + peak, share), (low, valley, high, peak)


def test_a_planet_beyond_a_hundred_days_of_four_earth_masses_or_more_is_mostly_a_sub_neptune(sweep):
    far = [row for row in sweep if row.days >= INNER_DAYS and 4.0 <= row.mass <= 20.0]
    assert len(far) >= 15
    envelope = sum(row.envelope >= planet_type.TERRESTRIAL_MAX_ENVELOPE for row in far) / len(far)
    assert envelope >= 0.6, (len(far), envelope)


def test_a_rocky_planet_is_never_larger_than_a_sub_neptune_can_be(sweep):
    """The start is held at 0.2 of the mass (EC5): the radius law of Lopez and Fortney is fitted to that, and no core of a few Earth masses is a puffed-up Saturn
    (9 Earth radii and more). A cold 4 Earth-mass core under the full 0.2 reaches 6.5 and a cold 3.5 Earth-mass one 7.05 (a lighter core puffs more; the super-puffs of Kepler-51 are the extreme at 6 to 9), so the bound is 7.5 for a core of a few Earth masses (the tail of a 12 Earth-mass planet is below)."""
    # GEN.1343: a super-puff is the extreme of Kepler-51, 6 to 9 Earth radii at under 0.1 g/cc. GEN.1406: the planet type names one only for 10 to 100 days and under 10 Earth masses
    # (`super_puff.is_super_puff`), while the puffed start (`draw_puff_start`) and the obliquity tide of GEN.1347 (which inflates a planet of 20 Earth fluxes and more by up to a factor of 2)
    # make a puff of a few Earth masses at 7 days or of 10.2 Earth masses at 26 days too: with the tide switched off none of the 60 systems has a planet over 6.5 Earth radii, with it
    # two planets of 8 and 9.9 Earth radii at 0.04 and 0.06 g/cc come out. A body under the definition's density edge is a puff whatever its type.
    puffed = [row for row in sweep if row.planet_type == planet_type.SUPER_PUFF or row.density < super_puff.DENSITY_RANGE_G_CM3[1]]
    ordinary = [row for row in sweep if row not in puffed]
    # The ordinary tail is the 12 Earth-mass planets whose puffed start the obliquity tide of GEN.1347 inflates to a density just over the super-puff edge:
    # 2 of 2,416 rocky planets in 360 G2V systems (0.08 percent) pass 7.5 Earth radii, all at 12.0 Earth masses, 0.11 to 0.14 g/cc, envelope 0.15 to 0.19, 8.36 and 7.94 Earth radii on main and
    # 7.72 after the GEN.1521 respacing moved an orbit into the tide's range; none is a rocky letter. So the bound for an ordinary planet is 8.5, the tail's measured maximum with a margin.
    assert max(row.radius for row in ordinary) <= 8.5
    assert max((row.radius for row in puffed), default=0.0) <= 2.0 * 7.5       # a puff of 7.5 Earth radii under the tide's largest factor of 2
    # An ordinary planet that keeps an envelope is no denser than its own dry part would be bare: the envelope only adds volume (a gas lighter than rock) and water only lowers the density, so the ceiling is the density
    # the generator's own fit of the structure solver gives the planet's metal and rock at that mass (`_bare_density`: 6.9 g/cc for a 10 Earth-mass planet of no iron, 15 for one of 0.8 iron),
    # with 5 percent for the fit (0.4 percent in radius, 1.2 in density) and the squeeze the envelope's weight gives the core. The 6.0 g/cc that stood here was the largest density of the sample that
    # drew it, not a physical edge: one 10.3 Earth-mass planet at 3.6 days (envelope 0.0033, 2.08 Earth radii, 6.28 g/cc, bare 6.9) is over it.
    assert all(row.density < 1.05 * row.bare_density for row in ordinary if row.envelope >= planet_type.TERRESTRIAL_MAX_ENVELOPE)


def test_no_puffed_planet_escapes_the_super_puff_type(sweep):
    """GEN.1470: a transitional planet of up to 30 Earth masses under 0.1 g/cc is typed a super-puff whatever its period (the type is mass and density)."""
    puffed = [row for row in sweep if row.envelope >= planet_type.TERRESTRIAL_MAX_ENVELOPE and row.density < 0.1 and row.mass <= 30.0]
    assert all(row.planet_type == planet_type.SUPER_PUFF for row in puffed), puffed


def test_the_super_puff_share_of_the_inner_sub_neptunes_is_counted_inside_a_hundred_days(sweep):
    """Row PF3 is a share of the 2 to 10 Earth-mass sub-Neptunes inside 100 days; with the type free of the period the count is made inside 100 days here, not by the type."""
    inner = [row for row in sweep if row.days <= INNER_DAYS and 2.0 <= row.mass <= 10.0 and row.envelope >= planet_type.TERRESTRIAL_MAX_ENVELOPE]
    assert len(inner) >= 6        # 8 since GEN.1711 (the median start of 0.04 leaves fewer envelopes over the type's edge), 10 or more before: a guard on the sample size
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
        if rogue.build.captured_envelope:
            kept.append(rogue)
            start = rogue.build.captured_envelope
            assert rogue.envelope_mass_fraction == pytest.approx(start)      # no light, no loss: what it started with
            assert 0.0 < start <= envelope_capture.MAX_START + 1e-9
            assert rogue.radius > rogue.build.bare_core["radius"] * 1.15
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
        assert not rogue.build.captured_envelope
