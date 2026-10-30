# tests/test_irregular_population.py

"""
GEN.1077: how many irregular moons a planet holds and which are stored (`planetgen.physics.irregular_population`): the capture windows, the grinding, the order-statistic draw into tiers, and the head counts of
the four giants (GEN.1075 targets: about 8 above 5 km at Jupiter, about 1 above 40 km, 8, 5, 6 and 8 above the storage floor, rocky planets none).
"""

import math
import random

import pytest

from planetgen.physics import belt_ledger, irregular_population as pop

SUN = 1.989e30
JUPITER = dict(planet_mass_kg=1.898e27, envelope_fraction=0.8, distance_au=5.203)
SATURN = dict(planet_mass_kg=5.683e26, envelope_fraction=0.8, distance_au=9.537)
URANUS = dict(planet_mass_kg=8.681e25, envelope_fraction=0.15, distance_au=19.19)
NEPTUNE = dict(planet_mass_kg=1.024e26, envelope_fraction=0.15, distance_au=30.07)
AGE = 4.5e9
MEAN_ENCOUNTERS = 129
FLOOR_RADIUS_KM = 10.6              # 1e16 kg at 2000 kg/m3, the radius table of CM6


class Stream:
    """A stand-in for `draw.Stream`: `random()` from a seeded generator."""

    def __init__(self, seed):
        self._rng = random.Random(seed)

    def random(self):
        return self._rng.random()


def population(planet, encounters=MEAN_ENCOUNTERS, seed=1, **more):
    return pop.irregular_population(star_mass_kg=SUN, encounters=encounters, age_yr=AGE, stream=Stream(seed), **planet, **more)


def mass_of(radius_km):
    return pop._radius_mass_kg(radius_km, 2000.0)


# --- who captures ---------------------------------------------------------------------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("mass_kg, envelope, distance", [(5.972e24, 0.0, 1.0), (1.0e25, 0.005, 1.0), (1.3e22, 0.0, 40.0), (3.0e25, 0.0, 0.1)])
def test_rocky_planets_and_dwarf_planets_capture_nothing(mass_kg, envelope, distance):
    p = population(dict(planet_mass_kg=mass_kg, envelope_fraction=envelope, distance_au=distance))
    assert p.windows == pop.NO_WINDOWS and p.survivors is None and p.ledger is None
    assert p.count_above(1.0) == 0.0 and p.stored_above(1.0) == 0


def test_a_hot_giant_with_its_regular_moons_beyond_its_stable_limit_captures_nothing():
    p = population(dict(JUPITER, distance_au=0.05), outermost_regular_moon_au=0.002)
    assert p.windows == pop.NO_WINDOWS and p.ledger is None


def test_a_warm_giant_captures_in_proportion_to_its_small_hill_volume_and_stores_nothing():
    warm = population(dict(JUPITER, distance_au=0.5))
    assert 0.0 < warm.windows.total < 1.0e-3 * population(JUPITER).windows.total
    assert warm.count_above(5.0) < 0.01 and warm.ledger.masses_kg == ()


def test_a_bigger_hill_sphere_than_jupiters_does_not_capture_more_per_encounter():
    ten = population(dict(JUPITER, planet_mass_kg=1.898e28))
    assert ten.windows.encounter == pytest.approx(population(JUPITER).windows.encounter)


# --- the windows ----------------------------------------------------------------------------------------------------------------------------------------------------------------

def test_the_encounter_channel_is_the_yield_per_encounter_and_gives_the_per_planet_probabilities_of_the_papers():
    for encounters in (59, 280):
        windows = pop.capture_windows(**JUPITER, star_mass_kg=SUN, encounters=encounters)
        assert windows.encounter == pytest.approx(encounters * pop.ENCOUNTER_YIELD)
    # CM2: 1.3e-8 to 3.6e-8 per particle at Jupiter in the histories of 59, 280 and 49 encounters (the 280 one is the high one)
    assert 1.0e-8 < pop.capture_windows(**JUPITER, star_mass_kg=SUN, encounters=59).total < 2.0e-8
    assert 3.0e-8 < pop.capture_windows(**JUPITER, star_mass_kg=SUN, encounters=280).total < 7.0e-8


def test_the_gas_channel_is_a_fifth_of_the_mean_encounter_channel_and_only_for_a_massive_envelope():
    no_instability = pop.capture_windows(**JUPITER, star_mass_kg=SUN, encounters=0)
    assert no_instability.encounter == 0.0
    assert no_instability.gas == pytest.approx(0.2 * pop.capture_windows(**JUPITER, star_mass_kg=SUN, encounters=129).encounter)
    assert pop.capture_windows(**URANUS, star_mass_kg=SUN, encounters=0).total == 0.0           # an ice giant with no instability holds none
    assert pop.capture_windows(**URANUS, star_mass_kg=SUN, encounters=129).gas == 0.0
    assert pop.capture_windows(**SATURN, star_mass_kg=SUN, encounters=0).gas > 0.0


def test_exchange_is_an_ice_giant_event_with_an_instability():
    assert pop.capture_windows(**NEPTUNE, star_mass_kg=SUN, encounters=129).triton_class == pytest.approx(pop.TRITON_EXCHANGE_PROBABILITY)
    assert pop.capture_windows(**NEPTUNE, star_mass_kg=SUN, encounters=0).triton_class == 0.0
    assert pop.capture_windows(**JUPITER, star_mass_kg=SUN, encounters=129).triton_class == 0.0


@pytest.mark.parametrize("kwargs", [dict(encounters=-1), dict(envelope_fraction=1.5)])
def test_nonsense_arguments_raise(kwargs):
    arguments = dict(JUPITER, star_mass_kg=SUN, encounters=10)
    arguments.update(kwargs)
    with pytest.raises(ValueError):
        pop.capture_windows(**arguments)


# --- the reservoir and the grinding ---------------------------------------------------------------------------------------------------------------------------------------------

def test_the_outer_disc_has_the_counted_bodies_and_a_mass_inside_the_papers_range():
    disc = pop.outer_disc()
    assert belt_ledger.count_above(disc.law, disc.mass_kg, mass_of(40.0)) == pytest.approx(3.5e7, rel=1e-9)
    assert 20.0 < disc.mass_kg / pop.EARTH_MASS_KG < 35.0
    assert pop.outer_disc(55.0).mass_kg == pytest.approx(2.0 * disc.mass_kg)


def test_the_grinding_reaches_32_km_at_the_age_of_the_solar_system_and_scales_with_the_collision_time():
    assert pop.full_radius_km(AGE) == pytest.approx(32.8, abs=0.5)
    assert pop.full_radius_km(0.0) == 0.0
    assert pop.full_radius_km(AGE, 2.0) == pytest.approx(pop.full_radius_km(AGE) * 2.0 ** -1.25)
    assert pop.full_radius_km(2.0e9) < pop.full_radius_km(4.0e9)


def test_a_young_planet_is_not_ground_and_has_the_reservoirs_own_count():
    disc = pop.outer_disc()
    young = pop.survivors(disc, 2.0e-8, 0.0)
    assert young.law == disc.law and young.count_above(40.0) == pytest.approx(2.0e-8 * 3.5e7)


def test_the_ground_population_has_the_observed_size_law_and_the_reservoirs_above_the_reach():
    disc = pop.outer_disc()
    probability = 2.0e-8
    left = pop.survivors(disc, probability, AGE)
    index = lambda lo, hi: -math.log(left.count_above(hi) / left.count_above(lo)) / math.log(hi / lo)
    assert index(1.0, 10.0) == pytest.approx(1.0)                                    # CM5 and CM17: cumulative 1.0 below the reach
    assert index(100.0, 400.0) == pytest.approx(5.0)                                 # CM3: the reservoir's above its break
    ungrounded = pop.survivors(disc, probability, 0.0)
    for radius in (left.full_radius_km, 70.0, 200.0):
        assert left.count_above(radius) == pytest.approx(ungrounded.count_above(radius), rel=1e-9)
    assert left.count_above(3.0) < ungrounded.count_above(3.0)
    assert 0.8 < left.mass_kg / ungrounded.mass_kg < 1.0                             # the mass is in the biggest bodies


def test_a_very_old_planet_is_ground_beyond_the_reservoirs_break_and_the_law_stays_valid():
    left = pop.survivors(pop.outer_disc(), 2.0e-8, 1.0e10)
    assert left.full_radius_km > 70.0 and len(left.law.exponents) == len(left.law.breaks_km) + 1
    assert left.law.exponents == pytest.approx((1.0, 4.2, 5.0))


# --- the head counts of the giants ----------------------------------------------------------------------------------------------------------------------------------------------

def test_a_jupiter_with_the_mean_history_holds_about_eight_above_5_km_and_one_above_40_km():
    p = population(JUPITER)
    assert p.count_above(5.0) == pytest.approx(8.0, rel=0.25)                       # GEN.1075: 8, not the design's 10 to 20
    assert p.count_above(40.0) == pytest.approx(1.0, rel=0.3)


def test_the_head_counts_of_the_four_giants_lie_inside_what_the_histories_span():
    # CM6: above 5 km 8, 11, 9, 8; above the storage floor 8, 5, 6, 8. The histories read have 49 to 280 encounters; a real giant's count is one Poisson draw about the mean of its history.
    low = [population(giant, 49).count_above(FLOOR_RADIUS_KM) for giant in (JUPITER, SATURN, URANUS, NEPTUNE)]
    high = [population(giant, 280).count_above(FLOOR_RADIUS_KM) for giant in (JUPITER, SATURN, URANUS, NEPTUNE)]
    for observed, lo, hi in zip((8, 5, 6, 8), low, high):
        assert lo < observed < hi + 2.0 * math.sqrt(hi)
    for giant, observed in ((JUPITER, 8), (SATURN, 11), (URANUS, 9), (NEPTUNE, 8)):
        high_mean = population(giant, 280).count_above(5.0)
        assert population(giant, 49).count_above(5.0) < observed < high_mean + 2.0 * math.sqrt(high_mean)


def test_the_four_giants_come_out_within_one_order_of_magnitude_of_each_other():
    counts = [population(giant).count_above(5.0) for giant in (JUPITER, SATURN, URANUS, NEPTUNE)]
    assert max(counts) / min(counts) < 10.0
    masses = [population(giant).survivors.mass_kg for giant in (JUPITER, SATURN, URANUS, NEPTUNE)]
    assert max(masses) / min(masses) < 10.0


def test_the_mass_of_a_jupiters_irregulars_is_the_observed_one_to_a_factor_two():
    assert 3.9e18 / 2.0 < population(JUPITER).survivors.mass_kg < 3.9e18 * 2.0         # CM6


# --- the draw into tiers --------------------------------------------------------------------------------------------------------------------------------------------------------

def test_the_stored_bodies_above_the_floor_average_the_expected_count():
    expected = population(JUPITER).count_above(FLOOR_RADIUS_KM)
    counts = [population(JUPITER, seed=seed).ledger.above_cut for seed in range(300)]
    assert sum(counts) / len(counts) == pytest.approx(expected, rel=0.12)
    assert 3.0 < expected < 6.0


def test_a_jupiter_has_a_body_above_40_km_in_the_share_of_histories_a_poisson_draw_gives():
    expected = population(JUPITER).count_above(40.0)
    share = sum(1 for seed in range(400) if population(JUPITER, seed=seed).stored_above(mass_of(40.0)) > 0) / 400.0
    assert share == pytest.approx(1.0 - math.exp(-expected), abs=0.08)


def test_the_stored_list_follows_the_storage_rule():
    for seed in range(60):
        ledger = population(JUPITER, seed=seed).ledger
        masses = ledger.masses_kg
        assert list(masses) == sorted(masses, reverse=True) and len(masses) <= pop.STORAGE_MAX_BODIES
        assert all(m >= pop.STORAGE_MIN_MASS_KG for m in masses) and not ledger.scaled
        assert sum(1 for m in masses if m >= pop.STORAGE_FLOOR_KG) == ledger.above_cut
        assert ledger.covered >= pop.STORAGE_COVERAGE or len(masses) == pop.STORAGE_MAX_BODIES
        assert ledger.belt_mass_kg == pytest.approx(sum(masses) + ledger.remainder_kg)


def test_the_draw_is_deterministic_takes_the_same_stream_whatever_the_planet_and_keeps_its_ranks():
    first, again = population(JUPITER, seed=7), population(JUPITER, seed=7)
    assert first.ledger.masses_kg == again.ledger.masses_kg
    after_giant, after_rock = Stream(9), Stream(9)
    pop.irregular_population(**JUPITER, star_mass_kg=SUN, encounters=129, age_yr=AGE, stream=after_giant)
    pop.irregular_population(5.972e24, 0.0, 1.0, SUN, 129, AGE, after_rock)
    assert after_giant.random() == after_rock.random()
    heavier = population(JUPITER, seed=7, reservoir=pop.outer_disc(60.0))
    assert heavier.ledger.masses_kg[0] > first.ledger.masses_kg[0]


# --- the shared ledger ----------------------------------------------------------------------------------------------------------------------------------------------------------

def test_the_mass_below_a_mass_adds_to_the_mass_above_it():
    law, belt = belt_ledger.ROCKY_LAW, 4.0e21
    assert belt_ledger.mass_below(law, belt, 1.0e-9) == pytest.approx(0.0, abs=1e-6 * belt)
    assert belt_ledger.mass_below(law, belt, 1.0e30) == pytest.approx(belt, rel=1e-3)
    previous = 0.0
    for mass in (1.0e12, 1.0e15, 1.0e18, 1.0e19, 1.0e21):
        below = belt_ledger.mass_below(law, belt, mass)
        assert below > previous
        previous = below
    # the mass above, by integrating the law, against the closed form
    low, steps = 1.0e18, 40000
    ratio = (1.0e30 / low) ** (1.0 / steps)
    total, m = 0.0, low
    for _ in range(steps):
        upper = m * ratio
        total += 0.5 * (m + upper) * (belt_ledger.count_above(law, belt, m) - belt_ledger.count_above(law, belt, upper))
        m = upper
    assert belt - belt_ledger.mass_below(law, belt, low) == pytest.approx(total, rel=2e-3)


def test_a_belt_ledger_with_the_default_rule_is_unchanged_by_the_new_arguments():
    exponentials = belt_ledger.draw_exponentials(Stream(3))
    law, belt = belt_ledger.ROCKY_LAW, 4.0e21
    plain = belt_ledger.ledger_from(law, belt, exponentials)
    named = belt_ledger.ledger_from(law, belt, exponentials, belt_ledger.L1_MASS_KG, belt_ledger.LEDGER_COVERAGE, belt_ledger.LEDGER_MAX_BODIES, 0.0, False)
    assert plain == named and plain.belt_mass_kg == belt
