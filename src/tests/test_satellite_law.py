# tests/test_satellite_law.py

"""
GEN.369 and GEN.624: a galaxy's satellites are drawn from its seed and its halo by stated laws, and scale with the host's
halo (none beyond its virial radius, none at all for a galaxy with no halo). The Milky Way's own rows
(`tests/reference_satellites.py`) are the plausibility target, not the source.
"""

import math
import statistics

import pytest

from planetgen.galaxy import mass_chain, satellite_history
from planetgen.generation import plan_outline, satellites
from tests import reference_satellites


def _seeds(count):
    return [bytes([i % 256, i // 256]) + bytes(14) for i in range(count)]


def _population(stellar_mass, count=40):
    drawn = []
    for seed in _seeds(count):
        chain = mass_chain.draw_mass_chain(seed, stellar_mass)
        drawn.append((chain, satellites.draw_satellites(seed, chain)))
    return drawn


def test_the_subhalo_count_is_the_laws():
    assert satellites.subhalo_count_above(1.0e-5) == pytest.approx(640.0, rel=0.02)     # [R0508]: about 640
    assert satellites.subhalo_count_above(0.01) == pytest.approx(0.95, rel=0.02)
    assert satellites.subhalo_count_above(0.001) > satellites.subhalo_count_above(0.01)


def test_occupation_follows_the_erf_about_its_centre():
    assert satellites.occupation_probability(5.0e5) == 0.0
    assert satellites.occupation_probability(satellites.OCCUPATION_M50_SOLAR) == pytest.approx(0.5)
    assert satellites.occupation_probability(3.0e8) == pytest.approx(1.0, abs=1e-6)


def test_the_two_size_laws_meet_where_they_hand_over():
    join = 10.0 ** satellites.SIZE_JOIN_LOG_MASS
    below = 10.0 ** (satellites.CARLSTEN_INTERCEPT + satellites.CARLSTEN_SLOPE * satellites.SIZE_JOIN_LOG_MASS)
    assert satellites.half_light_radius_pc(join * 1.0001, 0.0) == pytest.approx(below, rel=0.05)
    assert 100.0 < satellites.half_light_radius_pc(1.0e4, 0.0) < 300.0
    assert satellites.mean_fe_h(1.0e6) == pytest.approx(-1.69) and satellites.mean_fe_h(1.0e8) == pytest.approx(-1.09)


def test_a_galaxy_with_no_halo_has_no_satellites():
    chain = mass_chain.MassChain(1.0e6, None, None, None)
    assert satellites.draw_satellites(bytes(16), chain) == []


def test_the_count_scales_with_the_hosts_halo_at_three_masses():
    small = [len(found) for _chain, found in _population(3.0e7)]
    middle = [len(found) for _chain, found in _population(3.0e9)]
    milky_way = [len(found) for _chain, found in _population(5.0e10)]
    assert statistics.median(small) <= 3 and max(small) <= 8
    assert statistics.median(small) < statistics.median(middle) < statistics.median(milky_way)
    assert 20 <= statistics.median(milky_way) <= 150


def test_a_small_galaxy_pays_for_no_satellite_systems():
    """A 3e7 solar mass galaxy no longer holds the Milky Way's 6.4e9 satellite systems (PERF.146)."""
    for _chain, found in _population(3.0e7):
        assert sum(c.mass_solar for c in found) < 0.05 * 3.0e7


def test_no_satellite_lies_beyond_the_virial_radius_and_the_radius_is_the_halos():
    for stellar_mass in (3.0e9, 5.0e10):
        for chain, found in _population(stellar_mass, count=15):
            r200 = chain.virial_radius_kpc * 1000.0
            for component in found:
                distance = math.sqrt(sum(c * c for c in component.center_pc))
                assert satellites.INNER_OVER_R200 * r200 * 0.999 <= distance <= r200 * 1.001


def test_the_milky_way_mass_population_is_plausible_against_the_milky_ways_own_rows():
    """The Milky Way has 12 satellites above 1e5 solar masses [R2972]; the drawn median of those above 1e5 is of that order."""
    reference = sum(1 for row in reference_satellites.MILKY_WAY_SATELLITES if row.stellar_mass_sol > 1.0e5)
    drawn = [sum(1 for c in found if c.mass_solar > 1.0e5) for _chain, found in _population(5.0e10, count=60)]
    assert reference / 3.0 <= statistics.median(drawn) <= reference * 3.0
    # a Large Magellanic Cloud class satellite (above 1e9) is a minority outcome, not a given
    heavy = sum(1 for _chain, found in _population(5.0e10, count=60) if any(c.mass_solar > 1.0e9 for c in found))
    assert 0 < heavy < 30


def test_the_same_seed_draws_the_same_satellites_and_two_seeds_two():
    chain = mass_chain.draw_mass_chain(bytes(16), 5.0e10)
    first, again = satellites.draw_satellites(bytes(16), chain), satellites.draw_satellites(bytes(16), chain)
    assert first == again and first
    other_seed = bytes([1]) + bytes(15)
    other = satellites.draw_satellites(other_seed, mass_chain.draw_mass_chain(other_seed, 5.0e10))
    assert [c.center_pc for c in first] != [c.center_pc for c in other]
    assert len({c.seed for c in first}) == len(first)
    assert [c.mass_solar for c in first] == sorted((c.mass_solar for c in first), reverse=True)


def test_a_drawn_satellite_is_a_generated_row_with_a_profile_and_a_metallicity():
    _chain, found = _population(5.0e10, count=1)[0]
    assert found and {c.kind for c in found} == {"satellite"} and {c.source for c in found} == {"generated"}
    assert all(c.catalogue_id is None and c.velocity_kms == (0.0, 0.0, 0.0) for c in found)
    assert all(c.profile == ("exponential" if c.mass_solar >= satellites.EXPONENTIAL_MIN_STELLAR_MASS_SOLAR else "plummer")
               for c in found)
    assert all(-3.5 <= c.metallicity_fe_h <= 0.6 and 2.0 <= c.age_gy <= 13.0 and c.tidal_radius_pc > c.half_mass_radius_pc
               for c in found)
    assert all(c.mass_solar >= satellites.MIN_STELLAR_MASS_SOLAR for c in found)


def test_the_plan_draws_its_satellites_from_the_halo():
    built = plan_outline.build(plan_outline._parser().parse_args(["--mass", "3e7", "--disk-scale-length", "700pc"]),
                               bytes(16), 4.0)
    assert sum(c.mass_solar for c in built.satellites) < 0.05 * 3.0e7


def test_no_satellite_out_weighs_a_quarter_of_its_host():
    for stellar_mass in (3.0e9, 5.0e10):
        for chain, found in _population(stellar_mass, count=100):
            assert all(c.mass_solar <= satellites.SATELLITE_MAX_STELLAR_FRACTION * chain.stellar_mass_solar for c in found)


def test_the_stripping_radius_is_270_kpc_for_the_milky_ways_halo_and_scales_with_the_virial_radius():
    """GEN.1193: [R1671] put the cutoff of the dwarfs' HI at 250 to 280 kpc."""
    milky_way = mass_chain.virial_radius_kpc(mass_chain.median_halo_mass(5.0e10))
    assert satellites.stripping_radius_kpc(milky_way) == pytest.approx(270.0)
    assert satellites.stripping_radius_kpc(2.0 * milky_way) == pytest.approx(540.0)


def _forming_row(mass=1.0e8, radius=800.0):
    """A satellite that forms stars at the young end of the age range and not at the old end."""
    limit = satellite_history.oldest_forming_age_gy(mass, radius, *satellites.AGE_RANGE_GY)
    assert limit is not None and limit < satellites.AGE_RANGE_GY[1]
    return mass, radius, limit


def test_a_light_satellite_inside_the_stripping_radius_forms_no_stars_whatever_age_it_was_drawn():
    mass, radius, limit = _forming_row()
    inside = satellites.age_with_gas_by_distance(mass, radius, limit - 1.0, 100.0e3, 270.0e3, 0.0)
    assert inside > limit
    assert not satellite_history.history(mass, radius, inside).forming
    beyond = satellites.age_with_gas_by_distance(mass, radius, limit - 1.0, 400.0e3, 270.0e3, 0.0)
    assert beyond == limit - 1.0 and satellite_history.history(mass, radius, beyond).forming


def test_a_satellite_as_heavy_as_the_small_cloud_keeps_its_gas_inside_the_stripping_radius():
    mass = satellites.STRIPPING_RETAINED_MIN_STELLAR_MASS_SOLAR
    radius = satellites.half_light_radius_pc(mass, 0.0)
    youngest = satellites.AGE_RANGE_GY[0]
    assert satellite_history.history(mass, radius, youngest).forming
    assert satellites.age_with_gas_by_distance(mass, radius, youngest, 50.0e3, 270.0e3, 0.0) == youngest


def test_a_satellite_beyond_the_stripping_radius_keeps_its_gas_for_the_share_that_does():
    mass, radius, limit = _forming_row()
    old = limit + 2.0
    kept = satellites.age_with_gas_by_distance(mass, radius, old, 400.0e3, 270.0e3, 0.5)
    assert kept == limit and satellite_history.history(mass, radius, kept).forming
    lost = satellites.age_with_gas_by_distance(
        mass, radius, old, 400.0e3, 270.0e3, satellites.GAS_DETECTED_BEYOND_STRIPPING + 0.01)
    assert lost == old
    assert satellites.age_with_gas_by_distance(1.0e4, 100.0, 12.0, 400.0e3, 270.0e3, 0.0) == 12.0       # forms at no age


def test_no_drawn_satellite_inside_the_stripping_radius_forms_stars_unless_it_is_heavy():
    for chain, found in _population(5.0e10, count=12):
        radius_pc = satellites.stripping_radius_kpc(chain.virial_radius_kpc) * 1000.0
        for component in found:
            if math.sqrt(sum(c * c for c in component.center_pc)) < radius_pc and component.mass_solar < \
                    satellites.STRIPPING_RETAINED_MIN_STELLAR_MASS_SOLAR:
                assert not satellite_history.history(
                    component.mass_solar, component.half_mass_radius_pc, component.age_gy).forming, component.name
