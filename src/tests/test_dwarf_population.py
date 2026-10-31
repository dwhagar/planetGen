# tests/test_dwarf_population.py

"""GEN.1090: the dwarf planets of a belt or outer disc, from its mass and size law (`physics/dwarf_population.py`; belt-origins.md section 5)."""

import math
import random

import pytest

from planetgen.physics import belt_ledger as bl, body_strength as bs, constants, dwarf_population as dp, moon_composition

EARTH_KG = constants.EARTH_MASS_TO_KG
MAIN_BELT_KG = 4.008e-4 * EARTH_KG        # Pitjeva and Pitjev 2018; a test point only
KUIPER_BELT_KG = 1.97e-2 * EARTH_KG       # the same; a test point only


class Stream:
    def __init__(self, seed):
        self._rng = random.Random(seed)

    def random(self):
        return self._rng.random()


def test_the_count_above_a_diameter_is_the_size_laws():
    law = bl.ROCKY_LAW
    mass = math.pi / 6.0 * 150.0e3 ** 3 * law.density_kg_m3
    assert dp.count_above_km(law, MAIN_BELT_KG, 150.0) == pytest.approx(bl.count_above(law, MAIN_BELT_KG, mass))
    assert dp.count_above_km(law, MAIN_BELT_KG, 100.0) > dp.count_above_km(law, MAIN_BELT_KG, 200.0) > dp.count_above_km(law, 0.5 * MAIN_BELT_KG, 200.0)
    with pytest.raises(ValueError):
        dp.count_above_km(law, MAIN_BELT_KG, 0.0)


def test_the_typical_largest_body_has_one_body_above_it():
    for law, mass in ((bl.ROCKY_LAW, MAIN_BELT_KG), (bl.ICY_LAW, KUIPER_BELT_KG)):
        top_kg, top_km = dp.typical_largest(law, mass)
        assert bl.count_above(law, mass, top_kg) == pytest.approx(1.0)
        assert dp.count_above_km(law, mass, top_km) == pytest.approx(1.0)
        assert dp.chance_of_body_above(law, mass, top_kg) == pytest.approx(1.0 - math.exp(-1.0))


def test_the_chance_of_a_ceres_rises_with_the_belts_mass_and_is_small_at_the_main_belts():
    ceres = 9.4e20
    low = dp.chance_of_body_above(bl.ROCKY_LAW, MAIN_BELT_KG, ceres)
    high = dp.chance_of_body_above(bl.ROCKY_LAW, 10.0 * MAIN_BELT_KG, ceres)
    assert 0.0 < low < high < 1.0
    assert low < 0.2


def test_the_round_mass_is_the_engines_and_icy_bodies_turn_round_at_a_smaller_mass_than_rock():
    assert 1.0e19 <= dp.round_mass_kg(1.0) <= 1.0e20          # Mimas, 3.8e19 kg, is the smallest round ice body of the Solar System
    assert 1.0e20 <= dp.round_mass_kg(0.0) <= 1.0e21
    assert dp.round_mass_kg(0.0) > 3.0 * dp.round_mass_kg(0.5)
    assert bs.ratio_for(dp.round_mass_kg(0.3), 0.3) == pytest.approx(1.0, rel=1e-6)
    with pytest.raises(ValueError):
        dp.round_mass_kg(1.5)


def test_margots_discriminant_is_the_papers_for_jupiter_and_mars():
    """Margot 2015 prints Pi = 4.0e4 for Jupiter and 54 for Mars; the dwarf planets are under 1."""
    sun = 1.0
    assert dp.margot_pi(317.8 * EARTH_KG, 5.2, sun) == pytest.approx(4.0e4, rel=0.1)
    assert dp.margot_pi(0.107 * EARTH_KG, 1.524, sun) == pytest.approx(54.0, rel=0.1)
    assert dp.margot_pi(9.4e20, 2.77, sun) < 1.0          # Ceres
    assert dp.margot_clearing_mass_earth(1.0, 1.0) == pytest.approx(1.9e-4 * (2.0 * math.sqrt(3.0)) ** 1.5)
    assert dp.margot_clearing_mass_earth(5.0, 1.0) > dp.margot_clearing_mass_earth(1.0, 1.0)
    assert dp.margot_clearing_mass_earth(1.0, 0.5) < dp.margot_clearing_mass_earth(1.0, 1.0)


def test_soters_mu_puts_ceres_among_the_dwarfs_and_a_planet_in_a_clear_zone():
    assert dp.soter_mu(9.4e20, MAIN_BELT_KG - 9.4e20) == pytest.approx(0.6, abs=0.1)
    assert dp.status(9.4e20, 0.3, MAIN_BELT_KG - 9.4e20) == dp.DWARF      # Ceres, about a quarter to a third water by mass
    assert dp.status(6.0e24, 0.1, 1.0e22) == dp.PLANET
    assert dp.status(1.0e19, 0.3, 1.0e21) == dp.MINOR
    assert dp.status(5.0e20, 0.5, 1.0e22) == dp.DWARF
    with pytest.raises(ValueError):
        dp.soter_mu(1.0, 0.0)


def test_class_d_is_the_conditions_rule_on_a_stored_airless_round_body():
    assert dp.class_letter(1.0e21, 0.3) == "D"
    assert dp.class_letter(1.0e21, 0.0) == "D"             # round rock is D too: the water rule is gone (GEN.1160)
    assert dp.class_letter(2.0e19, 0.8) is None            # under a ratio of 1: an Asteroid, no letter
    assert dp.class_letter(1.0e20, 0.0) is None            # a rock of this mass is still not round
    assert dp.class_letter(5.0e18, 0.5) is None            # not stored: under the ledger's cut
    assert dp.class_letter(dp.AIRLESS_MAX_MASS_KG * 1.01, 0.5) is None
    assert dp.AIRLESS_MAX_MASS_KG == moon_composition.DWARF_PLANET_MAX_MASS_KG


def test_the_members_of_a_ledger_are_said_what_they_are():
    ledger = bl.draw_ledger(bl.ICY_LAW, KUIPER_BELT_KG, Stream(2))
    listed = dp.members(ledger, 0.4, 40.0, 1.0)
    assert len(listed) == len(ledger.masses_kg)
    assert [m.rank for m in listed] == list(range(len(listed)))
    assert [m.mass_kg for m in listed] == list(ledger.masses_kg)
    for m in listed:
        assert m.diameter_km == pytest.approx(bl.diameter_km(m.mass_kg, bl.ICY_LAW.density_kg_m3))
        assert (m.status == dp.MINOR) == (m.mass_kg < dp.round_mass_kg(0.4))
        assert m.class_letter == ("D" if m.mass_kg >= dp.round_mass_kg(0.4) else None)
        assert m.margot_pi < 1.0 and m.soter_mu < dp.SOTER_MU_PLANET
    dwarfs = dp.dwarf_planets(listed)
    assert all(m.status == dp.DWARF for m in dwarfs)
    assert [m.rank for m in dwarfs] == sorted(m.rank for m in dwarfs)


def test_the_ice_fraction_can_be_stated_for_each_body():
    ledger = bl.draw_ledger(bl.ROCKY_LAW, MAIN_BELT_KG, Stream(4))
    fractions = [0.0 if rank % 2 else 0.3 for rank in range(len(ledger.masses_kg))]
    listed = dp.members(ledger, fractions, 2.8, 1.0)
    assert [m.class_letter for m in listed] == ["D" if m.mass_kg >= bs.round_mass_for(f) else None for f, m in zip(fractions, listed)]


def test_the_main_belt_is_a_test_target():
    """Targets and never inputs: the main belt has 200 to 250 bodies of 100 km and up; the model's count at the belt's mass is within a factor of three, and a Ceres is a rare draw."""
    assert 225.0 / 3.0 < dp.count_above_km(bl.ROCKY_LAW, MAIN_BELT_KG, 100.0) < 225.0 * 3.0
    flagged = sum(1 for seed in range(300) if bl.draw_ledger(bl.ROCKY_LAW, MAIN_BELT_KG, Stream(seed)).masses_kg[0] > 3.2e20)
    assert 0 < flagged < 60


def test_the_kuiper_belt_is_a_test_target():
    """Targets and never inputs. The design note's 15 to 20 dwarf candidates of a Kuiper-type disc are the strongest icy bodies (the top of the strength range turns round at about 1e20 kg); at the median
    strength the engine turns bodies of a few 1e19 kg round (Mimas is round at 3.8e19 kg), and the count is a few hundred."""
    strong = bs.round_mass_for(1.0, quantile=1.0)
    assert 17.5 / 3.0 < bl.count_above(bl.ICY_LAW, KUIPER_BELT_KG, strong) < 17.5 * 3.0
    median = bl.count_above(bl.ICY_LAW, KUIPER_BELT_KG, dp.round_mass_kg(0.4))
    assert 100.0 < median < 500.0
    counts = sorted(len(dp.dwarf_planets(dp.members(bl.draw_ledger(bl.ICY_LAW, KUIPER_BELT_KG, Stream(s)), 0.4, 40.0, 1.0))) for s in range(200))
    assert 100 <= counts[len(counts) // 2] <= 500


def test_a_thin_belt_lists_none_and_mass_is_conserved():
    ledger = bl.draw_ledger(bl.ROCKY_LAW, 1.0e18, Stream(1))
    assert dp.dwarf_planets(dp.members(ledger, 0.0, 2.8, 1.0)) == ()
    assert sum(ledger.masses_kg) + ledger.remainder_kg == pytest.approx(1.0e18, rel=1e-12)
