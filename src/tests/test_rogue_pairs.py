# tests/test_rogue_pairs.py

"""GEN.332 slice D: a rogue planet may leave its donor bound to a companion, at the observed young share that fades with the age (`physics/rogue_pairs.py`), and the donor's
rogues stay a function of the donor alone (`rogue_yield.bind_pairs`)."""

import math

import pytest

from planetgen.generation import ejection_events, rogue_yield
from planetgen.physics import constants, rogue_pairs

JUPITER_KG = constants.JUPITER_MASS_TO_KG


def test_the_share_is_the_observed_young_share_and_fades_by_the_same_factor_every_10_myr():
    assert rogue_pairs.bound_share(0.0) == pytest.approx(0.09)
    assert rogue_pairs.survival(0.010) == pytest.approx(0.02)
    assert rogue_pairs.bound_share(0.020) == pytest.approx(0.09 * 0.02 ** 2)
    assert rogue_pairs.bound_share(1.0) < 1.0e-100      # a Gyr after the ejection nothing is left
    with pytest.raises(rogue_pairs.RoguePairError):
        rogue_pairs.survival(-1.0)


def test_the_separation_is_log_uniform_inside_the_observed_range():
    low, high = rogue_pairs.SEPARATION_RANGE_AU
    two_jupiters = 2.0 * JUPITER_KG
    assert rogue_pairs.separation_au(0.0, two_jupiters) == pytest.approx(low)
    assert rogue_pairs.separation_au(0.5, two_jupiters) == pytest.approx(math.sqrt(low * high))
    assert rogue_pairs.separation_au(0.999999, two_jupiters) == pytest.approx(high, rel=1e-4)
    with pytest.raises(rogue_pairs.RoguePairError):
        rogue_pairs.separation_au(1.0, two_jupiters)


def test_the_galaxys_tide_holds_two_jupiters_far_outside_the_observed_range_and_unbinds_a_featherweight_pair():
    assert 2.0e4 < rogue_pairs.tidal_limit_au(2.0 * JUPITER_KG) < 6.0e4          # about 35,000 AU
    assert rogue_pairs.separation_au(0.5, 1.0e12) is None                         # a pair of 1e12 kg cannot stay bound even at 28 AU
    assert rogue_pairs.tidal_limit_au(1.0e12) < rogue_pairs.SEPARATION_RANGE_AU[0]


def test_the_period_is_keplers_third_law():
    sun, earth = 1.988409870698051e30, 5.972e24
    assert rogue_pairs.period_yr(1.0, sun, earth) == pytest.approx(1.0, rel=1e-3)
    p = rogue_pairs.period_yr(100.0, 2.0 * JUPITER_KG, 3.0 * JUPITER_KG)
    assert p == pytest.approx(rogue_pairs.period_yr(100.0, 5.0 * JUPITER_KG, 1.0e-3), rel=1e-6)     # only the total mass counts
    assert rogue_pairs.period_yr(200.0, JUPITER_KG, JUPITER_KG) == pytest.approx(rogue_pairs.period_yr(100.0, JUPITER_KG, JUPITER_KG) * 2 ** 1.5)


STAR_A, STAR_B = object(), object()


def _event(index, since_gy, cause="thrown", mass_jupiters=2.0, kick=(5.0, 1.0, 0.0)):
    return ejection_events.EjectionEvent("uid", index, mass_jupiters * JUPITER_KG, 1.0e24, 5.0, 5.0, 3.0, 1.0e27, True, 0.001 * index, since_gy + 0.001, since_gy, kick, (1.0, 2.0, 3.0), 1.0e6, cause)


def _young(count, star=STAR_A):
    return [(star, _event(i, 0.0, kick=(float(i), 1.0, 0.0))) for i in range(count)]


def test_a_young_donors_rogues_pair_at_about_the_observed_share_and_the_two_are_one_system():
    found = _young(4000)
    events, pairs = rogue_yield.bind_pairs("donor", found)
    assert len(events) == len(found)
    share = len(pairs) / len(found)                  # 9 percent of the systems are pairs, a pair is two of the bodies: 2 * 0.09 / 1.09
    assert 0.12 < share < 0.21
    for index, link in pairs.items():
        other = pairs[link.partner_index]
        assert other.partner_index == index and other.separation_au == link.separation_au and other.period_yr == link.period_yr
        assert 28.0 <= link.separation_au <= 384.0 and link.period_yr > 0.0
    for first, second in ((a, b) for a, b in zip(events, events[1:]) if a.index in pairs and pairs[a.index].partner_index == b.index):
        assert (first.since_ejection_gy, first.kick_km_s, first.orbital_velocity_km_s, first.delay_gy) == (second.since_ejection_gy, second.kick_km_s, second.orbital_velocity_km_s, second.delay_gy)


def test_rogues_that_are_old_never_pair_and_their_events_are_untouched():
    found = [(STAR_A, _event(i, 1.0)) for i in range(2000)]
    events, pairs = rogue_yield.bind_pairs("donor", found)
    assert pairs == {} and events == [event for _star, event in found]


def test_the_pairs_of_a_donor_are_the_same_every_time_and_do_not_depend_on_the_uid_of_another():
    found = _young(500)
    assert rogue_yield.bind_pairs("donor", found) == rogue_yield.bind_pairs("donor", found)
    assert rogue_yield.bind_pairs("donor", found)[1] != rogue_yield.bind_pairs("another donor", found)[1]


def test_two_bodies_of_different_stars_or_different_causes_are_never_a_pair():
    split = [(STAR_A if i % 2 == 0 else STAR_B, _event(i, 0.0)) for i in range(2000)]
    assert rogue_yield.bind_pairs("donor", split)[1] == {}
    mixed = [(STAR_A, _event(i, 0.0, cause="thrown" if i % 2 == 0 else "cluster")) for i in range(2000)]
    assert rogue_yield.bind_pairs("donor", mixed)[1] == {}


def test_a_body_belongs_to_one_pair_only():
    _events, pairs = rogue_yield.bind_pairs("donor", _young(3000))
    partners = [link.partner_index for link in pairs.values()]
    assert len(partners) == len(set(partners)) == len(pairs)


def test_a_rogue_with_no_companion_has_no_pair_and_the_note_says_so_only_when_it_has_one():
    link = rogue_yield.Pair(7, 3.0 * JUPITER_KG, 120.0, 4200.0)
    assert "companion of 3.0 Jupiter masses, 120 AU away" in rogue_yield._pair_note(link) and "4,200 years" in rogue_yield._pair_note(link)
    assert rogue_yield.Rogue(None, None, None, "note").pair is None
