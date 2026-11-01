# tests/test_belt_history.py

"""GEN.1064: the events that act on a belt over a system's life (`physics/belt_history.py`; belt-evolution.md sections 2.5 and 4)."""

import random

import pytest

from planetgen.physics import belt_history as bh, belt_mass, constants, mass_loss as stellar_mass_loss, orbits

EARTH_KG = constants.EARTH_MASS_TO_KG


class Stream:
    def __init__(self, seed):
        self._rng = random.Random(seed)

    def random(self):
        return self._rng.random()

    def uniform(self, a, b):
        return self._rng.uniform(a, b)


def test_a_formed_belt_has_all_its_mass_and_no_history():
    state = bh.BeltState.formed(2.2, 3.4)
    assert (state.mass_fraction, state.scattered_fraction, state.history, state.lost) == (1.0, 0.0, (), False)
    assert bh.history_text(state) == ""
    with pytest.raises(ValueError):
        bh.BeltState.formed(3.0, 2.0)


def test_a_grand_tack_leaves_a_few_thousandths_and_an_instability_a_tenth_to_a_half():
    for seed in range(60):
        tack = bh.giant_event(bh.BeltState.formed(2.2, 3.4), bh.GRAND_TACK, Stream(seed))
        unstable = bh.giant_event(bh.BeltState.formed(2.2, 3.4), bh.INSTABILITY, Stream(seed))
        assert bh.GRAND_TACK_SURVIVAL[0] <= tack.mass_fraction <= bh.GRAND_TACK_SURVIVAL[1]
        assert bh.INSTABILITY_SURVIVAL[0] <= unstable.mass_fraction <= bh.INSTABILITY_SURVIVAL[1]
        assert len(tack.history) == 1 and "percent" in tack.history[0]
    with pytest.raises(ValueError):
        bh.giant_event(bh.BeltState.formed(2.2, 3.4), "flyby", Stream(1))


def test_the_survival_is_log_uniform_over_its_range():
    shares = sorted(bh.giant_event(bh.BeltState.formed(2.2, 3.4), bh.INSTABILITY, Stream(s)).mass_fraction for s in range(2000))
    assert shares[1000] == pytest.approx((0.1 * 0.5) ** 0.5, rel=0.1)


def test_the_stream_takes_the_same_draw_whatever_the_state():
    lost = bh.engulf(bh.BeltState.formed(1.0, 2.0), 5.0)
    assert lost.lost
    a, b = Stream(3), Stream(3)
    bh.giant_event(lost, bh.INSTABILITY, a)
    bh.giant_event(bh.BeltState.formed(2.0, 3.0), bh.INSTABILITY, b)
    assert a.random() == b.random()
    assert bh.giant_event(lost, bh.INSTABILITY, Stream(1)) is lost


def test_a_red_giant_swallows_the_inner_part_of_a_belt_and_the_outer_part_survives():
    """A belt of 1.5 to 4 au, an envelope that reaches 2.5 au: the part inside is lost, the rest keeps its share of the mass."""
    state = bh.engulf(bh.BeltState.formed(1.5, 4.0), 2.5)
    assert (state.lower_au, state.upper_au) == (2.5, 4.0)
    k = 2.0 - bh.SURFACE_DENSITY_SLOPE
    assert state.mass_fraction == pytest.approx((4.0 ** k - 2.5 ** k) / (4.0 ** k - 1.5 ** k))
    assert 0.0 < state.mass_fraction < 1.0
    assert "swelled out to 2.50 au" in state.history[0]
    assert not state.lost


def test_a_belt_wholly_inside_the_envelope_is_lost_and_one_wholly_outside_is_untouched():
    swallowed = bh.engulf(bh.BeltState.formed(1.0, 2.0), 3.0)
    assert swallowed.lost and swallowed.mass_fraction == 0.0 and "swallowed the whole belt" in swallowed.history[0]
    untouched = bh.BeltState.formed(5.0, 8.0)
    assert bh.engulf(untouched, 3.0) is untouched
    with pytest.raises(ValueError):
        bh.engulf(untouched, 0.0)


def test_the_engulfment_radius_comes_from_the_stars_own_module():
    assert hasattr(stellar_mass_loss, "engulfment_radius_au")


def test_a_white_dwarfs_belt_widens_by_the_mass_ratio_and_throws_a_share_inward_if_giants_remain():
    """A 5 solar mass star ends as a 0.9 solar mass white dwarf: a belt at 20 to 40 au goes to 111 to 222 au."""
    formed = bh.BeltState.formed(20.0, 40.0)
    with_giants = bh.mass_loss(formed, 5.0, 0.9, True)
    assert (with_giants.lower_au, with_giants.upper_au) == pytest.approx((20.0 * 5.0 / 0.9, 40.0 * 5.0 / 0.9))
    thrown = min(1.0, bh.MASS_LOSS_SCATTER_SHARE * (1.0 - 0.9 / 5.0))
    assert with_giants.scattered_fraction == pytest.approx(thrown)
    assert with_giants.mass_fraction == pytest.approx(1.0 - thrown)
    assert with_giants.mass_fraction + with_giants.scattered_fraction == pytest.approx(1.0)
    alone = bh.mass_loss(formed, 5.0, 0.9, False)
    assert alone.scattered_fraction == 0.0 and alone.mass_fraction == 1.0
    assert (alone.lower_au, alone.upper_au) == (with_giants.lower_au, with_giants.upper_au)
    assert "lost 82 percent of its mass" in with_giants.history[0]


def test_no_mass_loss_changes_nothing_and_a_gain_is_refused():
    state = bh.BeltState.formed(3.0, 5.0)
    assert bh.mass_loss(state, 1.0, 1.0, True) is state
    with pytest.raises(ValueError):
        bh.mass_loss(state, 1.0, 1.5, True)


def test_a_companion_truncates_the_outer_belt_at_the_holman_wiegert_limit():
    """A companion of half the pair's mass on a circular orbit at 60 au: stable out to 0.464 - 0.38 x 0.5 = 0.274 of the separation, 16.4 au."""
    limit_au = orbits.holman_wiegert_critical_semimajor_axis(60.0, 0.5, 0.0)
    assert limit_au == pytest.approx(60.0 * (0.464 - 0.38 * 0.5))
    state = bh.limit(bh.BeltState.formed(10.0, 40.0), 0.0, limit_au, "a companion at 60 au")
    assert state.upper_au == pytest.approx(limit_au) and state.lower_au == 10.0
    assert 0.0 < state.mass_fraction < 1.0
    assert state.history[0].startswith("A companion at 60 au limits the belt")
    far = bh.limit(bh.BeltState.formed(1.0, 5.0), 0.0, limit_au, "a companion at 60 au")
    assert far.mass_fraction == 1.0 and far.history == ()


def test_a_circumbinary_belt_starts_beyond_the_inner_limit_and_one_wholly_inside_it_is_lost():
    inner = orbits.holman_wiegert_circumbinary_a_crit_au(1.0, 0.5, 0.0)
    state = bh.limit(bh.BeltState.formed(1.0, 6.0), inner, 1.0e9, "a close companion")
    assert state.lower_au == pytest.approx(inner) and state.upper_au == 6.0
    gone = bh.limit(bh.BeltState.formed(1.0, 1.5), inner, 1.0e9, "a close companion")
    assert gone.lost and "was lost" in gone.history[0]
    with pytest.raises(ValueError):
        bh.limit(state, 5.0, 2.0, "x")


def test_the_events_in_order_keep_the_history_and_conserve_mass():
    state = bh.BeltState.formed(2.0, 6.0)
    state = bh.giant_event(state, bh.INSTABILITY, Stream(4))
    state = bh.engulf(state, 3.0)
    state = bh.mass_loss(state, 3.0, 0.7, True)
    assert len(state.history) == 3
    assert bh.history_text(state).count(".") >= 3
    assert state.mass_fraction + state.scattered_fraction <= 1.0
    assert bh.mass_kg(state, 1.0e22) == pytest.approx(state.mass_fraction * 1.0e22)
    assert bh.scattered_mass_kg(state, 1.0e22) == pytest.approx(state.scattered_fraction * 1.0e22)


def test_scattered_rock_makes_the_observed_pollution_rates():
    """[R2415]: a few 1e5 to a few 1e8 g/s. A belt of 4e17 to 3e20 kg scattered over 100 Myr spans it."""
    assert bh.pollution_rate_g_s(1.0e18) == pytest.approx(1.0e21 / (100.0 * 3.15576e13), rel=1e-9)
    assert bh.is_polluted(bh.pollution_rate_g_s(1.0e18))
    assert not bh.is_polluted(bh.pollution_rate_g_s(1.0e14))
    assert 1.0e5 < bh.pollution_rate_g_s(4.0e17) < 1.0e6
    assert bh.pollution_rate_g_s(3.0e20) < 3.0e9
    assert bh.pollution_rate_g_s(2.0e18, 50.0) == pytest.approx(2.0 * bh.pollution_rate_g_s(2.0e18, 100.0))
    with pytest.raises(ValueError):
        bh.pollution_rate_g_s(-1.0)


def test_a_main_belt_mass_belt_at_a_white_dwarf_is_a_test_target():
    """A test and not an input: the Solar System's belt (4e-4 Earth masses), widened and scattered at its Sun's death, makes a rate near the top of the observed range."""
    formed_kg = 4.0e-4 * EARTH_KG
    state = bh.mass_loss(bh.BeltState.formed(2.2, 3.4), 1.0, 0.54, True)
    rate = bh.pollution_rate_g_s(bh.scattered_mass_kg(state, formed_kg))
    assert 1.0e7 < rate < 1.0e9
    assert bh.is_polluted(rate)


def test_the_mass_of_a_belt_that_survived_goes_into_the_grinding_unchanged():
    """The history gives the mass left of what formed; `belt_mass` grinds that with the age."""
    state = bh.giant_event(bh.BeltState.formed(2.2, 3.4), bh.INSTABILITY, Stream(7))
    left = belt_mass.belt_mass_earth(state.mass_fraction * 2.0, 2.2, 3.4, 1.0, 4500.0, 500.0)
    assert 0.0 < left < state.mass_fraction * 2.0
