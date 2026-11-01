"""
The core's own thermal history (GEN.838, `planetgen.physics.core_heat`) beside the mantle box of `interior_heat`: the core's temperature, the
heat across the core-mantle boundary, the inner core and the dynamo of the Earth and of the worlds around it.
"""

import pytest

from planetgen.physics import constants, core_heat, interior_heat as ih, radiogenic

EARTH = constants.EARTH_MASS_TO_KG
EARTH_AGE_GY = 4.568
MID = ih.Interior(1.0e12, 250.0)


def _content():
    k, th, u = radiogenic.BULK_SILICATE_EARTH
    return radiogenic.Content(k * radiogenic.formation_mix("K")[1], th * radiogenic.formation_mix("Th")[1],
                              u * radiogenic.formation_mix("U")[1], radiogenic.SOLAR_AL26_RATIO)


def _heat(mass_earth=1.0, radius_km=6371.0, metal=0.325, age=EARTH_AGE_GY, lid=ih.STAGNANT, interior=MID):
    rock = radiogenic.Rock(1.0 - metal, metal, 0.0, radiogenic.BULK_SILICATE_EARTH_AL_MASS_FRACTION)
    return ih.heat_flow(_content(), interior, mass_earth * EARTH, radius_km, age, rock, None, lid)


def _stepper(mass_earth=1.0, radius_km=6371.0, metal=0.325):
    dry_m = radius_km * 1000.0
    return core_heat.Stepper(MID, mass_earth * EARTH * metal, mass_earth * EARTH * (1.0 - metal), dry_m * metal ** 0.5, dry_m)


def test_the_pressures_of_the_shells_give_the_earths_boundary_and_centre():
    earth = _stepper()
    assert earth.cmb_gpa == pytest.approx(136.0, rel=0.03)
    assert earth.centre_gpa == pytest.approx(364.0, rel=0.1)
    assert earth.liquidus_k(core_heat.LIQUIDUS_ANCHOR[0]) == pytest.approx(core_heat.LIQUIDUS_ANCHOR[1])
    assert earth.liquidus_k(0.0) == pytest.approx(core_heat.LIQUIDUS_ZERO_K)


def test_the_earths_core_has_an_inner_core_a_boundary_flow_and_a_running_dynamo():
    core = _heat().core
    assert 3800.0 < core.temperature_k < 4300.0           # 4,000 K or so at the top of the core
    assert 0.25 < core.inner_fraction < 0.5               # the inner core's 1,220 km of the core's 3,480 km
    assert 4.0e12 < core.cmb_w < 1.5e13                   # 5 to 15 TW across the boundary ([R2589], second-hand)
    assert core.cmb_w > core.adiabatic_w > 0.0
    assert core.running and core.dynamo_ended_gy is None
    assert core.dipole_scale == pytest.approx(1.0, abs=0.15)


def test_the_core_cools_and_freezes_as_the_body_ages():
    cores = [_heat(age=age).core for age in (0.5, 1.5, 3.0, 4.5)]
    temperatures = [core.temperature_k for core in cores]
    assert temperatures == sorted(temperatures, reverse=True)
    fractions = [core.inner_fraction for core in cores]
    assert fractions == sorted(fractions) and fractions[0] == 0.0 < fractions[-1]
    assert cores[0].running


def test_the_mantle_box_does_not_read_the_core(monkeypatch):
    with_core = _heat()
    monkeypatch.setattr(core_heat.Stepper, "advance", lambda self, *args: None)
    without = _heat()
    assert (with_core.total_w, with_core.mantle_temperature_k, with_core.released_j, with_core.urey_ratio) == \
        (without.total_w, without.mantle_temperature_k, without.released_j, without.urey_ratio)


def test_a_body_with_no_metal_has_no_core():
    rock = radiogenic.Rock(1.0, 0.0, 0.0, radiogenic.BULK_SILICATE_EARTH_AL_MASS_FRACTION)
    found = ih.heat_flow(_content(), MID, 0.1 * EARTH, 3000.0, 4.0, rock)
    assert found.core is None and found.total_w > 0.0
    assert ih.heat(radiogenic.heat_terms(_content(), 0.0), MID, EARTH, 6371.0, 1.0, radiogenic.Rock(0.0, 0.0, 1.0, 0.0)).core is None


def test_the_inner_core_follows_the_temperature_between_all_liquid_and_all_solid():
    earth = _stepper()
    assert earth.inner_radius_m(6000.0) == 0.0
    assert earth.inner_radius_m(2000.0) == earth.core_m
    radii = [earth.inner_radius_m(kelvin) for kelvin in (4400.0, 4200.0, 4000.0, 3800.0)]
    assert radii == sorted(radii) and 0.0 == radii[0] < radii[-1] <= earth.core_m


def test_the_core_loses_heat_to_a_cooler_mantle_and_gains_it_from_a_hotter_one():
    earth = _stepper()
    cold, hot = earth._flux(4000.0, 1500.0), earth._flux(4000.0, 3000.0)
    assert cold > 0.0 > hot
    assert earth._flux(4000.0, 1500.0 + 0.0) > earth._flux(4000.0, 2500.0) > 0.0        # a smaller drop carries less
    assert earth._flux(earth.boundary_k(2000.0), 2000.0) == 0.0


def test_the_freezing_iron_slows_the_cooling(monkeypatch):
    earth = _stepper()
    temperature_k = 4000.0
    assert earth.inner_radius_m(temperature_k) > 0.0
    with_latent = abs(earth._cooling(temperature_k, 1500.0)[0])
    monkeypatch.setattr(core_heat, "LATENT_HEAT_J_KG", 0.0)
    monkeypatch.setattr(core_heat, "GRAVITATIONAL_ENERGY_J_KG", 0.0)
    assert abs(earth._cooling(temperature_k, 1500.0)[0]) > with_latent


def test_a_mars_sized_core_has_stopped_its_dynamo_and_had_one_when_young():
    old = _heat(mass_earth=0.107, radius_km=3390.0, metal=0.22).core
    assert not old.running and old.dipole_scale == 0.0
    assert old.dynamo_ended_gy is not None and 0.0 < old.dynamo_ended_gy < EARTH_AGE_GY
    young = _heat(mass_earth=0.107, radius_km=3390.0, metal=0.22, age=0.3).core
    assert young.running and young.dynamo_ended_gy is None


def test_a_thin_liquid_shell_stops_the_dynamo(monkeypatch):
    earth = _stepper()
    monkeypatch.setattr(core_heat, "MIN_LIQUID_SHELL", 0.9)
    assert earth._power(4000.0, 2000.0) == 0.0


def test_a_super_earths_core_stays_liquid_and_its_dynamo_runs():
    heavy = _heat(mass_earth=5.0, radius_km=10000.0).core
    assert heavy.inner_fraction == 0.0 and heavy.running
    assert heavy.radius_m > _heat().core.radius_m
    assert heavy.dipole_scale > 1.0


def test_the_lid_decides_the_mantle_the_core_sits_under_and_not_whether_earth_has_a_field():
    stagnant, mobile = _heat(lid=ih.STAGNANT).core, _heat(lid=ih.MOBILE).core
    assert stagnant.boundary_k > mobile.boundary_k
    assert stagnant.running and mobile.running


def test_a_tide_that_melts_the_mantle_still_gives_a_core():
    """Io on Jupiter's orbit: the melt-limited branch of the box hands back the core of the run that warmed it."""
    from planetgen.physics import tidal_heating
    forcing = tidal_heating.Forcing(1.898e27, 4.217e8, 0.0041)
    rock = radiogenic.Rock(0.675, 0.325, 0.0, radiogenic.BULK_SILICATE_EARTH_AL_MASS_FRACTION)
    found = ih.heat_flow(_content(), MID, 0.015 * EARTH, 1822.0, 4.5, rock, forcing)
    assert found.tidal_w > 0.0 and found.core is not None
