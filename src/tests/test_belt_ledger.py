# tests/test_belt_ledger.py

"""GEN.1061: the ledger of a belt, the stored bodies of the mass cut and the mass left over (`physics/belt_ledger.py`; belt-evolution.md sections 3.1 and 11)."""

import math
import random

import pytest

from planetgen.physics import belt_ledger as bl, belt_mass, constants, formation

EARTH_KG = constants.EARTH_MASS_TO_KG
MAIN_BELT_KG = 4.008e-4 * EARTH_KG        # Pitjeva and Pitjev 2018; a test point only
KUIPER_BELT_KG = 1.97e-2 * EARTH_KG       # the same, with the 31 largest objects; a test point only
KENYON_BROMLEY_MAX_KG = 4.0 / 3.0 * math.pi * (1.75e6) ** 3 * 1000.0     # the largest icy body of a young disc: 1 750 km radius, 1 g/cm3


class Stream:
    """A stand-in for `draw.Stream`: `random()` from a seeded generator."""

    def __init__(self, seed):
        self._rng = random.Random(seed)

    def random(self):
        return self._rng.random()


def ledgers(law, mass_kg, count=300):
    return [bl.draw_ledger(law, mass_kg, Stream(seed)) for seed in range(count)]


def test_the_law_adds_up_to_the_belts_mass():
    """Integrating `m dN` over the stretches of the law gives the belt's mass: numerically, over a fine grid of masses."""
    for law, mass in ((bl.ROCKY_LAW, MAIN_BELT_KG), (bl.ICY_LAW, KUIPER_BELT_KG)):
        low, high = 1.0e5, 1.0e26
        steps = 60000
        ratio = (high / low) ** (1.0 / steps)
        total, m = 0.0, low
        for _ in range(steps):
            upper = m * ratio
            total += 0.5 * (m + upper) * (bl.count_above(law, mass, m) - bl.count_above(law, mass, upper))
            m = upper
        assert total == pytest.approx(mass, rel=0.02)


def test_the_inverse_is_the_inverse_and_the_law_is_continuous():
    for law in (bl.ROCKY_LAW, bl.ICY_LAW):
        for count in (0.5, 1.0, 3.0, 40.0, 3.0e3, 1.0e6):
            mass = bl.mass_at_count(law, MAIN_BELT_KG, count)
            assert bl.count_above(law, MAIN_BELT_KG, mass) == pytest.approx(count, rel=1e-9)
        for d_km in law.breaks_km:
            m = math.pi / 6.0 * (d_km * 1000.0) ** 3 * law.density_kg_m3
            assert bl.count_above(law, MAIN_BELT_KG, m * (1 + 1e-9)) == pytest.approx(bl.count_above(law, MAIN_BELT_KG, m * (1 - 1e-9)), rel=1e-6)
        masses = [bl.mass_at_count(law, MAIN_BELT_KG, c) for c in (1.0, 2.0, 10.0, 100.0, 1.0e4)]
        assert masses == sorted(masses, reverse=True)


def test_a_law_must_have_finite_mass_at_both_ends():
    with pytest.raises(ValueError):
        bl.SizeLaw((100.0,), (3.5, 4.0), 1000.0)
    with pytest.raises(ValueError):
        bl.SizeLaw((100.0,), (1.0, 2.5), 1000.0)
    with pytest.raises(ValueError):
        bl.SizeLaw((100.0, 50.0), (1.0, 3.0, 4.0), 1000.0)


def test_the_ledger_conserves_the_belts_mass():
    for law, mass in ((bl.ROCKY_LAW, MAIN_BELT_KG), (bl.ICY_LAW, KUIPER_BELT_KG), (bl.ROCKY_LAW, 1.0e16), (bl.ICY_LAW, 1.0e25)):
        for ledger in ledgers(law, mass, 60):
            assert sum(ledger.masses_kg) + ledger.remainder_kg == pytest.approx(mass, rel=1e-12)
            assert ledger.remainder_kg >= (1.0 - bl.LEDGER_MAX_SHARE) * mass * (1.0 - 1e-12)
            assert list(ledger.masses_kg) == sorted(ledger.masses_kg, reverse=True)
            assert 1 <= len(ledger.masses_kg) <= bl.LEDGER_MAX_BODIES


def test_every_body_above_the_cut_is_stored_and_the_rule_stops_at_the_coverage_or_the_cap():
    for ledger in ledgers(bl.ROCKY_LAW, MAIN_BELT_KG, 100):
        stored_above = sum(1 for m in ledger.masses_kg if m >= bl.L1_MASS_KG)
        assert stored_above == min(ledger.above_cut, bl.LEDGER_MAX_BODIES) or ledger.scaled
        if len(ledger.masses_kg) > ledger.above_cut and not ledger.scaled:
            # bodies were added below the cut: only because the ones above it held under the coverage, and it stopped at the coverage or the cap
            assert sum(ledger.masses_kg[:ledger.above_cut]) < bl.LEDGER_COVERAGE * MAIN_BELT_KG
            assert ledger.covered >= bl.LEDGER_COVERAGE or len(ledger.masses_kg) == bl.LEDGER_MAX_BODIES
            assert len(ledger.masses_kg) == bl.LEDGER_MAX_BODIES or sum(ledger.masses_kg[:-1]) < bl.LEDGER_COVERAGE * MAIN_BELT_KG


def test_a_belt_whose_mass_sits_in_thousands_of_bodies_is_held_to_the_cap_with_most_of_its_mass_in_the_rest():
    """The main belt's mass is mostly in bodies of 50 to 200 km, thousands of them: 300 bodies cannot cover 90 percent, and the rest is the L2 and L3 tiers."""
    ledger = bl.draw_ledger(bl.ROCKY_LAW, MAIN_BELT_KG, Stream(4))
    assert len(ledger.masses_kg) == bl.LEDGER_MAX_BODIES
    assert ledger.covered < bl.LEDGER_COVERAGE
    assert ledger.remainder_kg > 0.0


def test_a_light_belt_with_nothing_above_the_cut_adds_bodies_until_ninety_percent_is_covered_or_the_cap():
    ledger = bl.draw_ledger(bl.ROCKY_LAW, 1.0e17, Stream(1))
    assert ledger.above_cut == 0
    assert len(ledger.masses_kg) > 1
    assert ledger.covered >= bl.LEDGER_COVERAGE or len(ledger.masses_kg) == bl.LEDGER_MAX_BODIES


def test_a_rich_disc_is_held_to_the_cap():
    ledger = bl.draw_ledger(bl.ICY_LAW, KUIPER_BELT_KG, Stream(3))
    assert len(ledger.masses_kg) == bl.LEDGER_MAX_BODIES
    assert ledger.above_cut >= bl.LEDGER_MAX_BODIES


def test_ranks_are_stable_when_the_belt_loses_mass():
    """The same stream with a lighter belt: the same bodies in the same order, each lighter by the same factor above the top break."""
    heavy = bl.draw_ledger(bl.ROCKY_LAW, MAIN_BELT_KG, Stream(11))
    light = bl.draw_ledger(bl.ROCKY_LAW, 0.5 * MAIN_BELT_KG, Stream(11))
    factor = 0.5 ** (3.0 / bl.ROCKY_LAW.exponents[-1])
    top_break_kg = math.pi / 6.0 * (bl.ROCKY_LAW.breaks_km[-1] * 1000.0) ** 3 * bl.ROCKY_LAW.density_kg_m3
    checked = 0
    for a, b in zip(heavy.masses_kg, light.masses_kg):
        if b > top_break_kg:
            assert b / a == pytest.approx(factor, rel=1e-9)
            checked += 1
    assert checked >= 1
    assert all(b <= a for a, b in zip(heavy.masses_kg, light.masses_kg))


def test_the_same_stream_gives_the_same_ledger_and_takes_the_same_draws_whatever_the_mass():
    first, second = Stream(5), Stream(5)
    assert bl.draw_ledger(bl.ICY_LAW, KUIPER_BELT_KG, first).masses_kg == bl.draw_ledger(bl.ICY_LAW, KUIPER_BELT_KG, second).masses_kg
    other = Stream(5)
    bl.draw_ledger(bl.ICY_LAW, 1.0e20, other)
    assert other.random() == first.random()


def test_the_number_above_a_mass_follows_the_law():
    """The draw is a Poisson process in the law: the mean number of drawn bodies above a mass is `N(>m)`."""
    samples = []
    for seed in range(400):
        exps = bl.draw_exponentials(Stream(seed))
        samples.append(sum(1 for m in bl.ledger_from(bl.ROCKY_LAW, MAIN_BELT_KG, exps).masses_kg if m >= 3.0e19))
    mean = sum(samples) / len(samples)
    assert mean == pytest.approx(bl.count_above(bl.ROCKY_LAW, MAIN_BELT_KG, 3.0e19), rel=0.15)


def test_the_solar_system_main_belt_is_a_test_target():
    """At the main belt's mass, tens of bodies are above 1e19 kg (the design's 10 to 30) and a Ceres-size body (9.4e20 kg) is a rare but possible top."""
    counts = sorted(l.above_cut for l in ledgers(bl.ROCKY_LAW, MAIN_BELT_KG))
    assert 8 <= counts[len(counts) // 2] <= 30
    tops = sorted(l.masses_kg[0] for l in ledgers(bl.ROCKY_LAW, MAIN_BELT_KG))
    assert 3.0e19 < tops[len(tops) // 2] < 1.0e21
    assert tops[-1] > 3.0e20


def test_the_kuiper_belt_is_a_test_target():
    """At the Kuiper belt's mass: hundreds of bodies above the cut (so the cap holds), a largest body of 1e21 to 1e22 kg, none beyond a young disc's 1 750 km."""
    tops = sorted(l.masses_kg[0] for l in ledgers(bl.ICY_LAW, KUIPER_BELT_KG))
    assert 3.0e20 < tops[len(tops) // 2] < 1.0e22
    assert sum(1 for t in tops if t > KENYON_BROMLEY_MAX_KG) <= 0.1 * len(tops)


def test_a_stored_body_cannot_outweigh_the_belt():
    for seed in range(200):
        ledger = bl.draw_ledger(bl.ICY_LAW, 1.0e21, Stream(seed))
        assert ledger.masses_kg[0] < ledger.belt_mass_kg


def test_the_largest_body_and_the_belts_mass_are_settled_together():
    sun_w = 3.828e26
    snow = formation.snow_line_au(sun_w)
    primordial = belt_mass.primordial_mass_earth(2.2, 3.4, constants.SOLAR_MASS_TO_KG, snow)
    for seed in range(5):
        result = bl.solve_belt(primordial, 2.2, 3.4, 1.0, 4500.0, bl.ROCKY_LAW, Stream(seed))
        assert result.ledger.belt_mass_kg == pytest.approx(result.mass_earth * EARTH_KG)
        assert result.largest_km == pytest.approx(result.ledger.largest_km)
        assert bl.diameter_km(result.ledger.masses_kg[0], bl.ROCKY_LAW.density_kg_m3) == pytest.approx(result.largest_km)
        again = belt_mass.belt_mass_earth(primordial, 2.2, 3.4, 1.0, 4500.0, result.largest_km)
        assert result.mass_earth == pytest.approx(again, rel=2 * bl.SOLVE_TOLERANCE)
        assert result.mass_earth < primordial
        assert sum(result.ledger.masses_kg) + result.ledger.remainder_kg == pytest.approx(result.ledger.belt_mass_kg, rel=1e-12)


def test_the_solar_systems_belt_comes_out_within_a_few_times_the_inferred_mass():
    """A test and not an input: 2.2 to 3.4 AU round a Sun-like star at 4.5 Gyr, drawn belts have a few tenths of a thousandth to a few thousandths of an Earth mass."""
    snow = formation.snow_line_au(3.828e26)
    primordial = belt_mass.primordial_mass_earth(2.2, 3.4, constants.SOLAR_MASS_TO_KG, snow)
    masses = sorted(bl.solve_belt(primordial, 2.2, 3.4, 1.0, 4500.0, bl.ROCKY_LAW, Stream(s)).mass_earth for s in range(40))
    assert 1.0e-4 < masses[len(masses) // 2] < 1.0e-2


def test_a_young_belt_keeps_its_primordial_mass_without_a_search():
    snow = formation.snow_line_au(3.828e26)
    primordial = belt_mass.primordial_mass_earth(2.2, 3.4, constants.SOLAR_MASS_TO_KG, snow)
    result = bl.solve_belt(primordial, 2.2, 3.4, 1.0, 0.0, bl.ROCKY_LAW, Stream(2))
    assert result.passes == 2      # the primordial ledger and the check that nothing was ground
    assert result.mass_earth == pytest.approx(primordial)


def test_the_diameter_of_a_body():
    assert bl.diameter_km(4.0 / 3.0 * math.pi * 500.0 ** 3 * 1000.0 * 1e9, 1000.0) == pytest.approx(1000.0)
    with pytest.raises(ValueError):
        bl.diameter_km(-1.0, 1000.0)
