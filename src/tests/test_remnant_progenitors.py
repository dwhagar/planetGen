# tests/test_remnant_progenitors.py

"""GEN.554: which massive stars leave a black hole and which a neutron star. The failure windows of initial mass hold the
failed share of the 8 to 150 solar mass stars, nest as the share grows, keep the verified islands, and every progenitor the
generator draws lies in its kind's windows (docs/design/systems/explosion-islands.md)."""

from types import SimpleNamespace

import pytest

from planetgen import tuning
from planetgen.galaxy import structure
from planetgen.generation import multistar
from planetgen.physics import black_hole_formation, imf, remnant_progenitors as progenitors
from planetgen.util import draw

BOTTOM, TOP = tuning.CORE_COLLAPSE_MASS_RANGE_SOL
SHARES = (0.0, 0.01, 0.03, 0.05, 0.07, 0.09, 0.15, 0.2, 0.32, 0.45, 1.0)


def _share(windows):
    return sum(imf.number(low, high) for low, high in windows) / imf.number(BOTTOM, TOP)


def _inside(windows, mass):
    return any(low <= mass <= high for low, high in windows)


@pytest.mark.parametrize("share", SHARES)
def test_the_failure_windows_hold_the_failed_share_of_the_massive_stars(share):
    assert _share(progenitors.failure_windows(share)) == pytest.approx(share, abs=1e-5)


@pytest.mark.parametrize("share", SHARES)
def test_the_windows_are_ordered_apart_and_inside_the_core_collapse_range(share):
    windows = progenitors.failure_windows(share)
    assert all(BOTTOM <= low < high <= TOP for low, high in windows)
    assert all(a[1] < b[0] for a, b in zip(windows, windows[1:]))


@pytest.mark.parametrize("kind", progenitors.KINDS)
def test_the_two_kinds_share_the_range_between_them(kind):
    share = 0.2
    other = progenitors.KINDS[1 - progenitors.KINDS.index(kind)]
    mine, theirs = progenitors.kind_windows(kind, share), progenitors.kind_windows(other, share)
    assert _share(mine) + _share(theirs) == pytest.approx(1.0, abs=1e-5)
    for mass in (8.5, 14.0, 23.5, 24.9, 30.0, 40.0, 60.0, 120.0):
        assert _inside(mine, mass) != _inside(theirs, mass)


def test_a_larger_share_keeps_every_mass_a_smaller_one_fails():
    previous = progenitors.failure_windows(SHARES[0])
    for share in SHARES[1:]:
        windows = progenitors.failure_windows(share)
        for low, high in previous:
            for fraction in (0.0, 0.5, 1.0):
                assert _inside(windows, low + fraction * (high - low))
        previous = windows


def test_the_verified_islands_are_the_windows_when_the_share_is_theirs():
    islands = tuning.EXPLOSION_FAILURE_ISLANDS_SOL
    share = _share(islands)
    assert 0.0 < share <= 0.07 + 1e-3, "[R3327]: 0 to 7% of the massive stars at solar metallicity"
    windows = progenitors.failure_windows(share)
    assert len(windows) == len(islands)
    for (low, high), (want_low, want_high) in zip(windows, islands):
        assert (low, high) == pytest.approx((want_low, want_high), rel=1e-3)


def test_the_widest_verified_share_reaches_down_to_the_lowest_island():
    """45% [R4251] fails every star above about 14.5 solar masses; [R1286] put the lowest island at 15."""
    low = progenitors.failure_windows(0.45)[0][0]
    assert 14.0 < low < 15.5
    assert progenitors.failure_windows(0.45)[-1][1] == TOP


def test_a_galaxy_with_no_failures_leaves_neutron_stars_only():
    assert progenitors.failure_windows(0.0) == ()
    assert progenitors.kind_windows("neutron-star", 0.0) == ((BOTTOM, TOP),)
    assert progenitors.sample_mass_sol("black-hole", failed_share=0.0) is None


def test_the_reference_share_is_the_middle_draws_solar_share_inside_the_verified_range():
    assert progenitors.reference_failed_share() == black_hole_formation.reference_failed_shares(structure.REFERENCE_STRUCTURE.neutron_star_unbound_share)[0]
    low, high = tuning.BLACK_HOLE_FAILED_SHARE_SOLAR_RANGE
    assert low <= progenitors.reference_failed_share() <= high
    assert progenitors.failure_windows() == progenitors.failure_windows(progenitors.reference_failed_share())


def test_a_progenitor_lies_in_its_kinds_windows_and_no_lighter_than_asked():
    share = 0.15
    with draw.bound(554):
        for kind in progenitors.KINDS:
            windows = progenitors.kind_windows(kind, share)
            for lowest in (0.0, 12.0, 30.0):
                masses = [progenitors.sample_mass_sol(kind, lowest, draw, share) for _ in range(200)]
                assert all(mass is not None and mass >= lowest - 1e-9 and _inside(windows, mass) for mass in masses)


def test_progenitors_follow_the_imf_over_the_windows():
    """The draw's share in each window is that window's share of the IMF."""
    share = 0.2
    windows = progenitors.kind_windows("black-hole", share)
    wide = max(windows, key=lambda w: w[1] - w[0])
    weight = imf.number(*wide) / sum(imf.number(*w) for w in windows)
    with draw.bound(5541):
        masses = [progenitors.sample_mass_sol("black-hole", 0.0, draw, share) for _ in range(4000)]
    assert sum(wide[0] <= mass <= wide[1] for mass in masses) / len(masses) == pytest.approx(weight, abs=0.03)


def test_no_progenitor_when_no_window_reaches_the_lowest_mass():
    assert progenitors.sample_mass_sol("neutron-star", TOP + 1.0, draw) is None
    assert progenitors.sample_mass_sol("black-hole", 100.0, draw, 0.01) is None


def test_a_seed_gives_the_same_progenitor():
    first = progenitors.sample_mass_sol("neutron-star", 0.0, draw.Stream(7))
    assert first == progenitors.sample_mass_sol("neutron-star", 0.0, draw.Stream(7))


def test_the_number_of_a_kind_is_its_share_of_the_massive_stars():
    total = imf.number(BOTTOM, TOP)
    assert progenitors.number("black-hole", 0.0, 0.25) / total == pytest.approx(0.25, abs=1e-5)
    assert progenitors.number("neutron-star", 0.0, 0.25) / total == pytest.approx(0.75, abs=1e-5)
    assert progenitors.number("neutron-star", TOP, 0.25) == 0.0


def test_a_black_hole_collapses_directly_when_its_progenitor_failed():
    inside = progenitors.failure_windows()[0]
    assert progenitors.fails(0.5 * (inside[0] + inside[1]))
    assert not progenitors.fails(BOTTOM)
    placed = SimpleNamespace(initial_mass_solar=0.5 * (inside[0] + inside[1]), mass_solar=5.0)
    assert multistar._placed_direct_collapse(placed)
    placed = SimpleNamespace(initial_mass_solar=None, mass_solar=5.0)
    assert not multistar._placed_direct_collapse(placed)
    assert multistar._placed_direct_collapse(SimpleNamespace(initial_mass_solar=None, mass_solar=25.0))


def test_the_remnants_the_multi_star_builder_draws_come_from_the_windows():
    with draw.bound(5542):
        for kind in progenitors.KINDS:
            windows = progenitors.kind_windows(kind)
            leaves = [multistar._remnant(kind, 1.0e10, 0.0, draw) for _ in range(100)]
            assert all(leaf is not None and _inside(windows, leaf.initial) for leaf in leaves)
        assert all(leaf.direct_collapse for leaf in (multistar._remnant("black-hole", 1.0e10, 0.0, draw) for _ in range(20)))
