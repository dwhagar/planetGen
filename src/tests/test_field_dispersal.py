# tests/test_field_dispersal.py

"""How an interstellar asteroid field spreads and ends (GEN.1100, `physics/field_dispersal.py`): radius and age are one number and reversible, the design's age table, the mass is conserved, freeze-out
radii and times of design section 5.3, the shear, the unbound belt's speed, and bodies that are pure functions of (seed, rank, age)."""

import math

import pytest

from planetgen.physics import constants, field_dispersal as fd

KM = 1000.0
YEAR = constants.SECONDS_PER_YEAR
EARTH = constants.EARTH_MASS_TO_KG
AU = constants.AU_M
LY = constants.LIGHTYEAR_M


@pytest.mark.parametrize("kms,to_63au_yr,to_1ly_yr", [(0.1, 3e3, 3e6), (1, 300, 3e5), (5, 60, 6e4), (10, 30, 3e4), (15, 20, 2e4)])
def test_the_age_at_a_radius_follows_the_designs_table(kms, to_63au_yr, to_1ly_yr):
    assert fd.age_yr_at_radius(63 * AU, kms * KM) == pytest.approx(to_63au_yr, rel=0.1)
    assert fd.age_yr_at_radius(LY, kms * KM) == pytest.approx(to_1ly_yr, rel=0.1)
    assert fd.retirement_age_yr(kms * KM) == pytest.approx(to_1ly_yr, rel=0.1)


def test_radius_and_age_are_one_number_and_reversible():
    for speed in (0.1 * KM, 3 * KM, 12 * KM):
        for r0 in (0.0, 1.4e7, 5e10):
            for age in (0.0, 1.0, 1e3, 2.5e5):
                radius = fd.rms_radius_m(age, speed, r0)
                assert radius >= r0
                assert fd.age_yr_at_radius(radius, speed, r0) == pytest.approx(age, rel=1e-9, abs=1e-9) or radius == r0
    assert fd.age_yr_at_radius(1.0, 5 * KM, 1e9) == 0.0                  # a radius under the initial one is the birth
    radii = [fd.rms_radius_m(a, 5 * KM) for a in (0, 10, 100, 1e4)]
    assert radii == sorted(radii)


def test_a_field_is_retired_past_a_light_year_and_ages_are_uniform_up_to_it():
    speed = 5 * KM
    limit = fd.retirement_age_yr(speed)
    assert not fd.is_retired(limit * 0.99, speed) and fd.is_retired(limit * 1.01, speed)
    assert fd.rms_radius_m(limit, speed) == pytest.approx(LY)
    assert fd.age_from_unit_yr(0.0, speed) == 0.0 and fd.age_from_unit_yr(0.5, speed) == pytest.approx(limit / 2)
    assert all(not fd.is_retired(fd.age_from_unit_yr(u / 100, speed), speed) for u in range(100))
    with pytest.raises(fd.FieldDispersalError):
        fd.age_from_unit_yr(1.0, speed)


def test_the_mass_is_conserved_as_the_field_spreads():
    mass = 0.1 * EARTH
    speed = 8 * KM
    for age in (10.0, 1e3, 1e4):
        radius = fd.rms_radius_m(age, speed)
        assert fd.column_density_kg_m2(mass, radius) * math.pi * radius ** 2 == pytest.approx(mass)
    ratio = fd.column_density_kg_m2(mass, fd.rms_radius_m(1e3, speed)) / fd.column_density_kg_m2(mass, fd.rms_radius_m(1e4, speed))
    assert ratio == pytest.approx(100.0, rel=0.01)                                  # falls as 1 / t^2


@pytest.mark.parametrize("diameter_m,radius_au,time", [(1e5, 6.5e-4, 5 * 3600), (1e3, 6.5e-3, 2 * 86400), (1.0, 0.2, 60 * 86400), (1e-3, 6.5, 6 * YEAR)])
def test_freeze_out_follows_the_designs_table(diameter_m, radius_au, time):
    assert fd.freeze_out_radius_m(EARTH, diameter_m) / AU == pytest.approx(radius_au, rel=0.1)
    assert fd.freeze_out_time_s(EARTH, diameter_m, 5 * KM) == pytest.approx(time, rel=0.2)
    assert fd.freeze_out_radius_m(EARTH, diameter_m / 2) > fd.freeze_out_radius_m(EARTH, diameter_m)       # smaller pieces collide out to a larger radius


def test_shear_is_small_beside_the_spreading_until_late():
    pc = constants.PARSEC_M
    assert fd.shear_displacement_m(pc, 1e7) / pc == pytest.approx(0.3, rel=0.05)
    assert fd.shear_displacement_m(0.0, 1e7) == 0.0
    with pytest.raises(fd.FieldDispersalError):
        fd.shear_displacement_m(pc, -1.0)


def test_an_unbound_belt_leaves_at_the_speed_of_the_design():
    msun = constants.SOLAR_MASS_TO_KG
    v = fd.unbound_belt_speed_ms(10 * msun, 1.4 * msun, 50 * AU)
    assert v / KM == pytest.approx(11.3, abs=0.1)                                 # 13.3 km/s orbital, G (M0 - 2 M_rem) / a at infinity
    assert fd.unbound_belt_speed_ms(10 * msun, 5.1 * msun, 50 * AU) == 0.0       # over half the mass stays: the belt stays bound
    assert fd.retirement_age_yr(v) == pytest.approx(2.7e4, rel=0.1)              # the design's 25 to 37 kyr at 8 to 12 km/s


def test_a_body_is_a_pure_function_of_seed_rank_and_age():
    a = fd.member_position_m("F1", 7, 5000.0, 6 * KM)
    assert a == fd.member_position_m("F1", 7, 5000.0, 6 * KM)
    assert a != fd.member_position_m("F1", 8, 5000.0, 6 * KM) and a != fd.member_position_m("F2", 7, 5000.0, 6 * KM)
    assert fd.member_velocity_ms("F1", 7, 6 * KM) == fd.member_velocity_ms("F1", 7, 6 * KM)
    assert fd.member_position_m("F1", 7, 0.0, 6 * KM) == (0.0, 0.0, 0.0)
    v = fd.member_velocity_ms("F1", 7, 6 * KM)
    young, old = fd.member_position_m("F1", 7, 100.0, 6 * KM), fd.member_position_m("F1", 7, 200.0, 6 * KM)
    assert old[0] == pytest.approx(2 * young[0]) and old[2] == pytest.approx(2 * young[2])        # linear in the age (x and z carry no shear)
    assert young[0] == pytest.approx(v[0] * 100.0 * YEAR)


def test_the_rms_radius_of_the_bodies_is_the_radius_of_the_field():
    speed, age = 6 * KM, 2e4
    radii = [math.sqrt(sum(c * c for c in fd.member_position_m("F1", k, age, speed))) for k in range(3000)]
    rms = math.sqrt(sum(r * r for r in radii) / len(radii))
    assert rms == pytest.approx(fd.rms_radius_m(age, speed), rel=0.06)           # shear at 6 km/s for 20 kyr is tiny


def test_an_annulus_is_a_thin_ring_that_spreads_on_straight_tangents():
    speed, a, age = 11 * KM, 50 * AU, 1000.0
    for rank in range(20):
        x, y, z = fd.member_position_m("SN1", rank, age, speed, shape="annulus", belt_radius_m=a)
        assert z == 0.0
        assert math.hypot(x, y) == pytest.approx(math.hypot(a, speed * age * YEAR), rel=1e-4)       # shear is far below this at 1 kyr
    assert fd.member_position_m("SN1", 3, 0.0, speed, shape="annulus", belt_radius_m=a) == pytest.approx(
        (a * math.cos(fd.member_azimuth_rad("SN1", 3)), a * math.sin(fd.member_azimuth_rad("SN1", 3)), 0.0))


def test_bad_arguments_are_refused():
    with pytest.raises(fd.FieldDispersalError):
        fd.rms_radius_m(-1.0, 1000.0)
    with pytest.raises(fd.FieldDispersalError):
        fd.rms_radius_m(1.0, 0.0)
    with pytest.raises(fd.FieldDispersalError):
        fd.member_velocity_ms("F", -1, 1000.0)
    with pytest.raises(fd.FieldDispersalError):
        fd.member_velocity_ms("F", 1, 1000.0, shape="torus")
    with pytest.raises(fd.FieldDispersalError):
        fd.member_position_m("F", 1, 1.0, 1000.0, shape="annulus")
