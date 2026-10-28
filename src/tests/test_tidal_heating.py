# tests/test_tidal_heating.py

"""
GEN.737: the tide of an eccentric orbit as a heat source of the body, for rocky and icy bodies (`planetgen.physics.tidal_heating`,
`interior_heat`): Io, Europa, Enceladus and the Moon, and the accounting of the heat a generated body stores.
"""

import math
import types

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import interior_heat, radiogenic, tidal_heating
from planetgen.physics import constants as pc
from planetgen.physics.mathcheck import _bulk_silicate_earth
from planetgen.util import draw
from tests.fuzz_support import deterministic_entropy

JUPITER = pc.JUPITER_MASS_TO_KG
SATURN = 5.683e26
EARTH = pc.EARTH_MASS_TO_KG
AGE_GY = 4.5
STIFF = interior_heat.Interior(1.0e12, 250.0)
"""Interior: A mantle rheology in the middle of the drawn ranges."""


def _rock(silicate, metal):
    return radiogenic.Rock(silicate, metal, 0.0, radiogenic.BULK_SILICATE_EARTH_AL_MASS_FRACTION)


def _heat(mass_kg, radius_km, rock, forcing, interior=STIFF):
    return interior_heat.heat_flow(_bulk_silicate_earth(), interior, mass_kg, radius_km, AGE_GY, rock, forcing)


def test_peale_power_has_the_published_scaling():
    forcing = tidal_heating.Forcing(JUPITER, 6.711e8, 0.0094)
    base = tidal_heating.peale_power_w(3.0e-4, forcing, 1560.8e3)
    assert tidal_heating.peale_power_w(6.0e-4, forcing, 1560.8e3) == pytest.approx(2.0 * base)
    assert tidal_heating.peale_power_w(3.0e-4, tidal_heating.Forcing(JUPITER, 6.711e8, 0.0188), 1560.8e3) == pytest.approx(4.0 * base)
    assert tidal_heating.peale_power_w(3.0e-4, forcing, 2 * 1560.8e3) == pytest.approx(32.0 * base)
    farther = tidal_heating.Forcing(JUPITER, 2 * 6.711e8, 0.0094)
    assert tidal_heating.peale_power_w(3.0e-4, farther, 1560.8e3) == pytest.approx(base / 2 ** 7.5)


def test_europa_gets_its_published_tide_from_its_real_eccentricity():
    """Europa: 1.4e11 W at e = 0.0094 (Peale, Cassen and Reynolds 1979 with the ice shell's k2/Q); the heat flow over its surface is 4.6e-3 W/m2."""
    forcing = tidal_heating.Forcing(JUPITER, 6.711e8, 0.0094)
    assert tidal_heating.peale_power_w(tidal_heating.ICY_K2_OVER_Q, forcing, 1560.8e3) == pytest.approx(1.4e11, rel=0.1)
    assert tidal_heating.icy_flux_w_m2(forcing, 1560.8) == pytest.approx(4.6e-3, rel=0.1)


def test_enceladus_gets_the_equilibrium_tide_and_not_its_unexplained_remainder():
    """Enceladus: the equilibrium tide of its eccentricity of 0.0047 is about 1 GW (Meyer and Wisdom 2007: 1.1 GW), its radiated 16 GW (Howett et
    al. 2011) is more than that tide makes, an open problem the old 1 to 50 resonance boost papered over. This is the tide alone."""
    forcing = tidal_heating.Forcing(SATURN, 2.3802e8, 0.0047)
    power = tidal_heating.icy_flux_w_m2(forcing, 252.1) * 4 * math.pi * (252.1e3) ** 2
    assert power == pytest.approx(1.1e9, rel=0.7)
    assert power < 1.5e10


def test_the_moon_is_flexed_solid_and_heated_by_a_couple_of_hundred_megawatts():
    """The Moon: k2 = 0.021 from its rigidity (observed 0.0242, Williams and Boggs 2015) and Q = 100 give about 2e8 W at e = 0.0549 (the lunar Q
    of 38 would give 6e8 W); it stays solid."""
    mass, radius_km = 7.346e22, 1737.4
    assert tidal_heating.rigid_k2(mass, radius_km * 1000.0) == pytest.approx(0.0242, rel=0.2)
    forcing = tidal_heating.Forcing(EARTH, 3.84399e8, 0.0549)
    found = _heat(mass, radius_km, _rock(0.98, 0.02), forcing)
    assert 1.0e8 < found.tidal_w < 1.0e9
    assert found.mantle_temperature_k < tidal_heating.HEAT_PIPE_TEMPERATURE_K
    assert found.tidal_w < 1.0e-2 * found.own_w
    assert found.own_flux_w_m2 == pytest.approx(found.flux_w_m2, rel=1.0e-2)


def test_io_is_melt_limited_and_radiates_about_two_watts_a_square_metre():
    """Io at e = 0.0041: its tide in the solid state (1.6e12 W) is more than its own heat (4e11 W) and warms its mantle to the melt that erupts, so
    it holds the heat-pipe temperature and radiates 9e13 W, 2.2 W/m2 (Veeder et al. 2004: 2.24 W/m2, 0.9 to 1.5e14 W)."""
    forcing = tidal_heating.Forcing(JUPITER, 4.218e8, 0.0041)
    found = _heat(8.93e22, 1821.6, _rock(0.8, 0.2), forcing)
    assert found.tidal_w == pytest.approx(9.0e13, rel=0.5)
    assert found.flux_w_m2 == pytest.approx(2.2, rel=0.5)
    assert found.mantle_temperature_k == pytest.approx(tidal_heating.HEAT_PIPE_TEMPERATURE_K)
    assert found.own_w < 0.05 * found.tidal_w


def test_a_soft_mantle_can_lose_the_tide_by_convection_before_it_melts():
    """The melt state depends on the mantle's drawn rheology (a soft one convects the 1.6e12 W away at a cooler temperature): the same Io then
    stays solid and its tide is the solid-state one."""
    forcing = tidal_heating.Forcing(JUPITER, 4.218e8, 0.0041)
    soft = _heat(8.93e22, 1821.6, _rock(0.8, 0.2), forcing, interior_heat.Interior(2.6e10, 200.0))
    assert soft.mantle_temperature_k < tidal_heating.HEAT_PIPE_TEMPERATURE_K
    assert soft.tidal_w == pytest.approx(tidal_heating.solid_power_w(forcing, 8.93e22, 1821.6e3))


def test_a_circle_and_no_host_give_no_tide_and_leave_the_heat_as_it_was():
    assert tidal_heating.forcing_of(types.SimpleNamespace(primary_mass_kg=EARTH, eccentricity=0.0, distance=0.00257)) is None
    assert tidal_heating.forcing_of(types.SimpleNamespace(primary_mass_kg=None, eccentricity=0.1, distance=1.0)) is None
    plain = _heat(EARTH, 6371.0, _rock(0.675, 0.325), None)
    assert plain.tidal_w == 0.0 and plain.own_w == plain.total_w
    earth = _heat(EARTH, 6371.0, _rock(0.675, 0.325), tidal_heating.Forcing(1.989e30, pc.AU_TO_M, 0.0167))
    assert 0.0 < earth.tidal_w < 1.0e-6 * earth.total_w
    assert earth.total_w == pytest.approx(plain.total_w, rel=1.0e-6)
    assert earth.urey_ratio == pytest.approx(plain.urey_ratio, rel=1.0e-6)


def test_the_heat_adds_up_and_the_urey_ratio_leaves_the_tide_out():
    forcing = tidal_heating.Forcing(JUPITER, 4.218e8, 0.0041)
    found = _heat(8.93e22, 1821.6, _rock(0.8, 0.2), forcing)
    assert found.total_w == pytest.approx(found.radiogenic_w + found.primordial_w + found.tidal_w)
    assert found.own_w == pytest.approx(found.total_w - found.tidal_w)
    assert found.urey_ratio == pytest.approx(found.radiogenic_w / found.own_w)


def test_more_eccentric_means_more_heat():
    last = 0.0
    for eccentricity in (0.001, 0.003, 0.01, 0.03):
        found = _heat(7.346e22, 1737.4, _rock(0.98, 0.02), tidal_heating.Forcing(EARTH, 3.84399e8, eccentricity))
        assert found.tidal_w > last
        last = found.tidal_w


def _bodies(system):
    for planet in system.planets + system.secondary_planets:
        if not hasattr(planet, "moons"):
            continue
        yield planet
        yield from planet.moons


def test_every_generated_body_stores_the_tide_of_its_orbit_and_the_habitability_energy_counts_it():
    from planetgen.physics import hydrosphere
    seen = icy = rocky = 0
    for seed in range(10):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(SystemConfig(), galactic_center_dist_ly=26000.0)
        for body in _bodies(system):
            if body.body_type in ("g", "a"):
                assert body.tidal_flux_w_m2 is None
                continue
            seen += 1
            forcing = tidal_heating.forcing_of(body)
            if forcing is None:  # a planet or captured moon the tide has circularised has none left
                assert body.eccentricity == 0.0 and not body.tidal_flux_w_m2
                continue
            if (body.water_mass_fraction or 0.0) >= hydrosphere.ICY_MIN_WATER_FRACTION:
                icy += 1
                assert body.tidal_flux_w_m2 == pytest.approx(tidal_heating.icy_flux_w_m2(forcing, body.radius))
            elif body.tidal_flux_w_m2 is not None:
                rocky += 1
                assert body.tidal_flux_w_m2 > 0.0
                area = 4 * math.pi * (body.radius * 1000.0) ** 2
                # never more than the melt-limited state makes (the box works on the body's dry radius, a few percent from the stored one)
                assert body.tidal_flux_w_m2 * area <= tidal_heating.molten_power_w(forcing, body.radius * 1000.0) * 1.25
            # the own heat flow does not contain the tide: the habitability energy term adds it
            assert body.internal_heat_flux_w_m2 is not None and body.internal_heat_flux_w_m2 > 0.0
    assert seen > 10 and rocky > 0
