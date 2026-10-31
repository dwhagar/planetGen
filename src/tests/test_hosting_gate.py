# tests/test_hosting_gate.py

"""
GEN.579: the hosting gate (`planetgen.physics.hosting_gate`): the three parts in order (a satellite must fit, a stable
orbit must exist, a chance that one exists at all), the size cut, the share table with its sources, a moon that has no row,
and the tidal-survival radius of a submoon (Kollmeier and Raymond 2019, equation 1) against the moons the paper names.
"""

import math

import pytest

from planetgen import tuning
from planetgen.generation import possibility_table as pt
from planetgen.physics import constants, hosting_gate as hg, relations

SUN = 1.989e30
AU = constants.AU_M


def body(diameter_m, population="belt-small-body", density=2000.0, **kw):
    mass = density * math.pi / 6.0 * diameter_m ** 3
    kw.setdefault("parent_mass_kg", SUN)
    kw.setdefault("orbit_radius_m", 2.5 * AU)
    return hg.Host(population, diameter_m, mass, **kw)


def test_a_two_kilometre_belt_body_passes_all_three_parts():
    verdict = hg.can_host(body(2000.0))
    assert verdict.allowed and verdict.part == 0 and "0.160" in verdict.reason


def test_a_body_below_the_size_cut_has_no_satellites_even_with_room():
    verdict = hg.can_host(body(500.0))
    assert not verdict.allowed and verdict.part == 3 and "size cut" in verdict.reason
    assert hg.satellite_share(body(500.0)) == 0.0


def test_a_host_no_bigger_than_the_floor_fails_part_one():
    verdict = hg.can_host(body(10.0))
    assert not verdict.allowed and verdict.part == 1
    verdict = hg.can_host(body(12.0))            # a 10 m satellite would outweigh the host past a tenth
    assert not verdict.allowed and verdict.part == 1


def test_the_floor_satellite_fits_a_host_whose_mass_is_ten_times_it():
    floor_mass = hg.smallest_satellite_mass_kg()
    d = hg.FLOOR_DIAMETER_M * (1.0 / hg.MASS_RATIO_LIMIT) ** (1.0 / 3.0)
    assert hg.satellite_fits(body(d * 1.001)) and not hg.satellite_fits(body(d * 0.99))
    assert floor_mass == pytest.approx(2000.0 * math.pi / 6.0 * 1000.0)


def test_a_body_close_to_a_star_has_no_stable_zone():
    # 1 km body at 0.001 au: its Hill radius is under its own size
    verdict = hg.can_host(body(2000.0, orbit_radius_m=0.001 * AU))
    assert not verdict.allowed and verdict.part == 2


def test_the_stable_zone_is_the_hill_fraction_and_the_radius_plus_floor():
    host = body(20_000.0)
    inner, outer = hg.stable_zone_m(host)
    assert inner == pytest.approx((20_000.0 + 10.0) / 2.0)
    hill = relations.hill_radius(host.orbit_radius_m, host.mass_kg, SUN)
    assert outer == pytest.approx(0.4895 * hill)
    assert hg.hill_radius_m(body(20_000.0, eccentricity=0.5)) == pytest.approx(0.5 * hill)


def test_a_free_body_has_no_outer_edge():
    host = hg.Host("belt-small-body", 2000.0, 1e13)
    assert hg.stable_zone_m(host)[1] is None and hg.stable_orbit_exists(host)


def test_a_fluid_host_has_the_roche_limit_as_its_inner_edge():
    host = body(20_000.0, fluid=True)
    inner, _ = hg.stable_zone_m(host)
    assert inner == pytest.approx(relations.ROCHE_FLUID_COEFFICIENT * 10_000.0)   # same density as the satellite
    assert inner > hg.stable_zone_m(body(20_000.0))[0]


def test_host_arguments_are_checked():
    with pytest.raises(hg.HostingError):
        hg.Host("belt-small-body", 0.0, 1.0)
    with pytest.raises(hg.HostingError):
        hg.Host("belt-small-body", 1.0, 1.0, parent_mass_kg=1.0)
    with pytest.raises(hg.HostingError):
        hg.Host("belt-small-body", 1.0, 1.0, parent_mass_kg=1.0, orbit_radius_m=1.0, eccentricity=1.0)


@pytest.mark.parametrize("population, diameter, share", [
    ("inner-small-body", 2000.0, 0.15), ("belt-small-body", 2000.0, 0.16),
    ("inner-small-body", 50_000.0, 0.05), ("belt-small-body", 200_000.0, 0.05),
    ("trans-neptunian-cold", 200_000.0, 0.22), ("trans-neptunian-hot", 200_000.0, 0.055),
    ("trans-neptunian-cold", 50_000.0, 0.0), ("trans-neptunian-hot", 5000.0, 0.0)])
def test_the_share_rows(population, diameter, share):
    host = body(diameter, population, orbit_radius_m=40 * AU)
    assert hg.satellite_share(host) == pytest.approx(share)


def test_a_population_the_table_does_not_cover_defers_to_its_host_class():
    host = body(2.0e6, "planet-like")
    assert hg.satellite_share(host) is None
    verdict = hg.can_host(host)
    assert verdict.allowed and "host-class table" in verdict.reason


def test_a_moon_hosts_nothing_until_a_row_exists():
    moon = hg.Host("moon", 3.0e6, 1e22, parent_mass_kg=1.9e27, orbit_radius_m=1.9e9)
    verdict = hg.can_host(moon)
    assert not verdict.allowed and verdict.part == 3 and "no" in verdict.reason
    assert hg.satellite_share(moon) == 0.0


def test_the_table_is_registered_with_a_source_and_state_on_every_row():
    table = pt.REGISTRY.get("satellite-share")
    assert [row.name for row in table.rows][-1] == "no satellites"
    assert all(row.source and row.state in (pt.PARTLY_VERIFIED, pt.STAND_IN) for row in table.rows)
    assert all(row.replacement for row in table.rows if row.state == pt.STAND_IN)
    for population in hg.SATELLITE_SHARE_CLASSES:
        for diameter in (1500.0, 20_000.0, 300_000.0):
            subject = pt.Subject(population, {"diameter_m": diameter})
            assert sum(p for _, p in table.probabilities(subject)) == pytest.approx(1.0, abs=1e-12)


def test_the_tidal_radius_matches_the_equation_and_the_paper_s_named_moons():
    jupiter, saturn, earth = 1.898e27, 5.683e26, 5.972e24
    cases = {"moon": (1737e3, 3344.0, earth, 384.4e6, True), "callisto": (2410e3, 2500.0, jupiter, 1.883e9, True),
             "titan": (2575e3, 2500.0, saturn, 1.222e9, True), "iapetus": (735e3, 2500.0, saturn, 3.561e9, True),
             "io": (1822e3, 2500.0, jupiter, 4.217e8, False), "europa": (1561e3, 2500.0, jupiter, 6.71e8, False),
             "rhea": (764e3, 2500.0, saturn, 5.27e8, False)}
    for name, (radius, rho, planet, a, passes) in cases.items():
        assert hg.moon_survives_tides(radius, rho, planet, a) is passes, name
    assert hg.tidal_survival_radius_m(2500.0, saturn, 1.222e9) / 1e3 == pytest.approx(353.6, rel=0.01)


def test_a_smaller_submoon_or_a_wider_orbit_needs_a_smaller_moon():
    base = hg.tidal_survival_radius_m(2500.0, 5.683e26, 1.222e9)
    assert hg.tidal_survival_radius_m(2500.0, 5.683e26, 1.222e9, submoon_mass_kg=1e15) < base
    assert hg.tidal_survival_radius_m(2500.0, 5.683e26, 2.4e9) < base
    assert hg.tidal_survival_radius_m(2500.0, 5.683e26, 1.222e9, tidal_q=1000.0) < base
    with pytest.raises(hg.HostingError):
        hg.tidal_survival_radius_m(0.0, 1.0, 1.0)
