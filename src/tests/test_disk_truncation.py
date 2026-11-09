# tests/test_disk_truncation.py

"""GEN.1451: a companion star cuts a star's disk, and so its belts, before their mass is worked out (`planetgen.physics.disk_truncation`, `belt_origin.build(disk_limit_au=...)`)."""

import math
from types import SimpleNamespace

import pytest

from planetgen.physics import belt_origin as bo, disk_truncation as dt, radiation_field as rf
from planetgen.physics.constants import SOLAR_MASS_TO_KG as SUN
from tests.test_belt_origin import SOL_LIKE, build, giant, rocky


def star(mass_kg=SUN):
    return SimpleNamespace(mass=mass_kg, luminosity=3.8e26)


def with_companion(host, axis_au, mass_kg=SUN, eccentricity=0.0):
    rf.attach(host, [rf.Companion(star(mass_kg), 1, axis_au, eccentricity)])
    return host


def test_the_roche_lobe_of_an_equal_pair_is_the_fit_of_eggleton():
    assert dt.roche_lobe_share(1.0) == pytest.approx(0.379, abs=0.002)
    assert dt.roche_lobe_share(10.0) > dt.roche_lobe_share(1.0) > dt.roche_lobe_share(0.1), "the heavier star holds the bigger lobe"


def test_a_disk_round_an_equal_pair_ends_near_a_quarter_of_the_periastron():
    radius = dt.truncation_radius_au(SUN, SUN, 20.0, 0.0)
    assert radius == pytest.approx(dt.TRUNCATION_OF_ROCHE_LOBE * dt.roche_lobe_share(1.0) * 20.0)
    assert 0.25 < radius / 20.0 < 0.31
    assert dt.truncation_radius_au(SUN, SUN, 20.0, 0.5) == pytest.approx(0.5 * radius), "the periastron sets it"


def test_a_lone_star_has_no_truncation():
    lone = star()
    rf.attach(lone, ())
    assert dt.star_truncation_au(lone) == math.inf
    assert dt.truncation_radius_au(SUN, 0.0, 20.0, 0.0) == math.inf


def test_the_tightest_companion_sets_the_radius():
    host = star()
    rf.attach(host, [rf.Companion(star(), 1, 80.0, 0.0), rf.Companion(star(0.2 * SUN), 2, 12.0, 0.1)])
    assert dt.star_truncation_au(host) == pytest.approx(dt.truncation_radius_au(SUN, 0.2 * SUN, 12.0, 0.1))


def test_no_limit_gives_the_belts_of_before_exactly():
    plain = build(SOL_LIKE)
    for limit in (math.inf, 1.0e9):
        again = build(SOL_LIKE, disk_limit_au=limit)
        assert [(b.origin, b.inner_au, b.outer_au, b.primordial_mass_earth, b.mass_earth) for b in again] == \
               [(b.origin, b.inner_au, b.outer_au, b.primordial_mass_earth, b.mass_earth) for b in plain]


def test_a_limit_cuts_the_outer_disk_and_its_mass_to_the_stretch_that_is_left():
    chain = [rocky(1.0), giant(5.2)]
    whole = [b for b in build(chain) if b.origin == "outer-disk"][0]
    cut = [b for b in build(chain, disk_limit_au=40.0) if b.origin == "outer-disk"][0]
    assert whole.outer_au == pytest.approx(100.0) and cut.outer_au == pytest.approx(40.0)
    assert cut.inner_au == whole.inner_au
    assert cut.primordial_mass_earth < whole.primordial_mass_earth
    assert bo.conservation_error(cut) < 1e-9


def test_a_limit_inside_a_stretch_leaves_no_belt_when_nothing_of_it_is_left():
    chain = [rocky(1.0), giant(5.2)]
    start = [b for b in build(chain) if b.origin == "outer-disk"][0].inner_au
    belts = build(chain, disk_limit_au=0.97 * start)
    assert all(b.outer_au <= 0.97 * start for b in belts)
    assert not [b for b in belts if b.origin == "outer-disk"], "the outer disk starts past the giant's chaotic zone, beyond the limit"


def test_a_planets_chaotic_zone_that_is_not_a_giants_lies_inside_its_keep_out():
    """The register row of GEN.1451: only giants clear spans from a belt, because for any smaller planet the zone is already inside the stretch no belt may take (five Hill radii)."""
    for mass_earth in (0.1, 1.0, 5.0, 20.0, 50.0):
        body = rocky(2.0, mass_earth * 5.97e24)
        low, high = bo.keep_out(body, SUN)
        half = bo.chaotic_half_width_au(body, SUN)
        assert body.distance_au + half <= high and body.distance_au - half >= low, mass_earth


def test_the_disk_limit_is_the_nearer_of_the_cut_and_the_stability_limit():
    host = with_companion(star(), 20.0)
    host.a_crit_au = 3.0
    assert dt.star_disk_limit_au(host) == 3.0
    host.a_crit_au = 50.0
    assert dt.star_disk_limit_au(host) == dt.star_truncation_au(host)
    lone = star()
    rf.attach(lone, ())
    assert dt.star_disk_limit_au(lone) == math.inf
    lone.a_crit_au = None
    assert dt.star_disk_limit_au(lone) == math.inf
