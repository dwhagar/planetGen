# src/tests/test_dark_light.py

"""GEN.1501: the dark members of an orbit tree give light: a brown dwarf by its cooling, a neutron star by its age, a black hole by the winds that feed it (`physics/dark_light.py`)."""

import math

import pytest

from planetgen.generation.stored_trees import DarkMember
from planetgen.physics import constants, dark_light as dl, radiation_field, stellar_wind

SOLAR_MASS = 1.989e30


def _member(kind, mass, radius_au=None, number=2):
    return DarkMember(number, 0, 1, kind, mass, radius_au)


class _Star:
    """A star as the wind law reads it."""
    def __init__(self, mass_sol=1.0, luminosity_sol=1.0, radius_km=695700.0, type="G2V Yellow Main Sequence Star", yerkes_class="V"):
        self.mass = mass_sol * SOLAR_MASS
        self.luminosity = luminosity_sol * constants.SOLAR_LUMINOSITY
        self.radius = radius_km
        self.type = type
        self.yerkes_class = yerkes_class


def test_the_wind_law_is_the_one_the_heliosphere_uses():
    from planetgen.generation.star import Star
    sun = _Star()
    mdot, speed = stellar_wind.wind(sun.mass, sun.luminosity, sun.radius, sun.type, sun.yerkes_class)
    assert 1e8 < mdot < 1e11 and 2e5 < speed < 1e6, "about the Sun's wind: 1e9 kg/s at a few hundred km/s"
    assert Star._calculate_heliosphere_radius_static(sun.mass, sun.luminosity, sun.radius, sun.type, sun.yerkes_class) > 0


def test_a_neutron_star_cools_with_age_and_the_curve_is_continuous_at_the_end_of_the_neutrino_era():
    ages = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 5.0]
    temperatures = [dl.surface_temperature_k(age) for age in ages]
    assert temperatures == sorted(temperatures, reverse=True)
    assert dl.surface_temperature_k(1e-4) == pytest.approx(dl.NEUTRON_STAR_PHOTON_ERA_T_K)
    before, after = dl.surface_temperature_k(1e-4 * (1 - 1e-9)), dl.surface_temperature_k(1e-4 * (1 + 1e-9))
    assert before == pytest.approx(after, rel=1e-6)
    assert 3e3 < dl.surface_temperature_k(5.0) < 1e4, "an old neutron star is a few thousand kelvin"
    assert dl.neutron_star_luminosity_w(12.0, 1.0) == pytest.approx(
        constants.STEFAN_BOLTZMANN_CONSTANT * 4 * math.pi * 12000.0 ** 2 * dl.surface_temperature_k(1.0) ** 4)


def test_a_neutron_star_is_younger_than_its_system_by_its_progenitors_life():
    star = dl.source_of(_member("neutron-star", 1.4, radius_au=12.0 / constants.AU_TO_KM), 5.0)
    assert star.age == pytest.approx(5.0 - dl.NEUTRON_STAR_PROGENITOR_LIFE_GY)
    assert star.state == "neutron-star" and star.luminosity > 0 and star.radius == pytest.approx(12.0)
    young = dl.source_of(_member("neutron-star", 1.4, radius_au=12.0 / constants.AU_TO_KM), 0.01)
    assert young.age == dl.NEUTRON_STAR_MIN_AGE_GY, "never under the youngest age the curve is read at"
    assert young.luminosity > star.luminosity * 1e3


def test_a_brown_dwarf_gives_the_light_of_its_cooling_at_the_age_of_the_system():
    from planetgen.physics import brown_dwarf
    source = dl.source_of(_member("brown-dwarf", 0.05), 5.0)
    state = brown_dwarf.state(0.05, 5.0)
    assert source.luminosity == pytest.approx(state["luminosity_sol"] * constants.SOLAR_LUMINOSITY)
    assert source.type.endswith("dwarf") and source.state == "brown-dwarf"
    assert dl.source_of(_member("brown-dwarf", 0.05), 0.5).luminosity > source.luminosity, "a younger one is brighter"


def test_a_black_hole_shines_only_when_a_wind_feeds_it_and_more_the_nearer_the_star_is():
    hole = _member("black-hole", 10.0, radius_au=1e-6)
    assert dl.source_of(hole, 5.0) is None, "nothing feeds it: it gives nothing"
    blue = _Star(mass_sol=20.0, luminosity_sol=5e4, radius_km=8e6, type="O5V", yerkes_class="V")
    near = dl.source_of(hole, 5.0, [(blue, 1.0)])
    far = dl.source_of(hole, 5.0, [(blue, 10.0)])
    assert near.luminosity > far.luminosity > 0
    assert near.luminosity == pytest.approx(far.luminosity * 100.0, rel=0.1), "the capture falls as the square of the separation while the wind outruns the orbit"
    assert near.state == "black-hole-disk" and near.temperature > 1e4


def test_what_a_hole_takes_is_never_more_than_the_wind_and_its_light_never_more_than_eddington():
    mass = 10 * SOLAR_MASS
    assert dl.wind_capture_kg_s(mass, 1e-6, 1e12, 3e6, SOLAR_MASS) <= 1e12
    assert dl.disk_luminosity_w(mass, 1e30) == pytest.approx(dl.EDDINGTON_W_PER_SOLAR_MASS * 10, rel=1e-3)
    assert dl.disk_luminosity_w(mass, 0.0) == 0.0
    assert dl.wind_capture_kg_s(mass, 0.0, 1e9, 3e5, SOLAR_MASS) == 0.0


def test_the_disk_is_taken_at_the_innermost_stable_orbit():
    mass = 10 * SOLAR_MASS
    assert dl.disk_radius_km(mass) == pytest.approx(3.0 * 2.0 * constants.G * mass / constants.SPEED_OF_LIGHT_M_S ** 2 / 1e3)
    assert dl.disk_radius_km(mass) == pytest.approx(88.6, rel=0.02), "88 km for ten solar masses"


def test_an_unlit_kind_or_a_member_with_no_mass_or_age_gives_none():
    assert dl.source_of(_member("stripped-helium", 1.0), 5.0) is None
    assert dl.source_of(_member("neutron-star", 1.4), 5.0) is None, "no radius"
    assert dl.source_of(_member("brown-dwarf", 0.05), None) is None
    assert dl.source_of(_member("brown-dwarf", 0.0), 5.0) is None


def test_a_dark_source_is_a_source_of_the_field_with_the_bands_it_has():
    source = dl.source_of(_member("neutron-star", 1.4, radius_au=12.0 / constants.AU_TO_KM), 0.05)
    per = radiation_field._per_au(source)
    assert per["bolometric_w_m2"] == pytest.approx(source.luminosity / (4 * math.pi * constants.AU_M ** 2))
    assert per["xuv_flux_earth"] is None and per["exposure_index"] is None and per["flare_irradiation_index"] is None


def _system_with_a_brown_dwarf():
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw
    for seed in (2, 74, 90, 128, 172, 276):         # seeds that drew one (a bounded search; any of them will do)
        config = SystemConfig()
        config.PLANETS = False
        with draw.bound(seed):
            system = StarSystem(system_config=config)
        if system.dark_members and system.dark_members[0].kind == "brown-dwarf":
            return system
    pytest.fail("no system with a brown dwarf in its tree among the seeds tried")


def _dark_companions(system):
    return [c for star in system.stars for c in radiation_field.companions_of(star) if c.dark]


def test_a_brown_dwarf_of_the_tree_is_a_dark_source_of_the_stars_it_hangs_beside():
    system = _system_with_a_brown_dwarf()
    found = _dark_companions(system)
    assert found and all(c.star.state == "brown-dwarf" and c.member_no >= len(system.stars) for c in found)
    host = system.star
    with_dark = radiation_field.field_at(host, 1.0).bolometric_w_m2
    assert any(c.member_no in with_dark.by_source for c in found if c in radiation_field.companions_of(host))


def test_a_neutron_star_in_the_tree_lights_the_planets_of_the_star_beside_it_and_a_lone_star_is_unchanged():
    import dataclasses
    system = _system_with_a_brown_dwarf()
    host = system.star
    before = radiation_field.starlight_w_m2(host, 1.0)
    system.dark_members = (dataclasses.replace(system.dark_members[0], kind="neutron-star", mass_solar=1.4, radius_au=12.0 / constants.AU_TO_KM),) + tuple(system.dark_members[1:])
    system.primary_star.age = 0.05          # the system is young: the neutron star is hot
    system.build_ledger()
    neutron = [c for c in radiation_field.companions_of(host) if c.dark]
    assert neutron and neutron[0].star.state == "neutron-star"
    assert radiation_field.starlight_w_m2(host, 1.0) > before
    assert radiation_field.starlight_w_m2(host, 1.0) == pytest.approx(before + sum(
        radiation_field.field_at(host, 1.0).bolometric_w_m2.by_source.get(c.member_no, 0.0) for c in neutron) - 0.0, rel=1.0)       # the total is the sum of its parts
    assert not any(c.dark for c in radiation_field.companions_of(_lone_star()))


def _lone_star():
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw
    config = SystemConfig()
    config.PLANETS = False
    config.BINARY_SYSTEM = False
    with draw.bound(5):
        return StarSystem(system_config=config).star


def test_the_truncation_of_a_disk_does_not_count_a_glowing_dark_member():
    from planetgen.physics import disk_truncation
    system = _system_with_a_brown_dwarf()
    star = system.star
    dark_only = [c for c in radiation_field.companions_of(star) if c.dark]
    stars_only = [c for c in radiation_field.companions_of(star) if not c.dark]
    assert dark_only
    radiation_field.attach(star, stars_only)
    without = disk_truncation.star_truncation_au(star)
    radiation_field.attach(star, stars_only + dark_only)
    assert disk_truncation.star_truncation_au(star) == without
