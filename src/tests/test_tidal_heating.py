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
from planetgen.physics import icy_tide, interior_heat, radiogenic, spin_torques, tidal_heating
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
COOL = interior_heat.Interior(1.0e13, 200.0)
"""Interior: One whose mantle runs cooler (a stiffer prefactor with a lower activation energy) than `STIFF` at the age of these bodies."""
CORNERS = [interior_heat.Interior(eta0, activation) for eta0 in (9.3e7, 1.9e11) for activation in (260.0, 390.0)]
"""list: The four corners of the drawn viscosity ranges (`interior_heat.draw_interior`)."""


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


def test_europa_gets_a_tide_of_the_published_size_from_its_ice_shell_with_nothing_fitted():
    """Europa: 1e11 to 1e12 W of tide at e = 0.0094 in a shell of 10 to 30 km ([R2060]; the old 1.4e11 W was a fitted k2/Q): the shell's own rheology
    (`icy_tide`, GEN.960) gives it from the temperature profile and the strain, with the radiogenic heat of its rock under it."""
    from planetgen.physics import hydrosphere, icy_tide
    flux, power = icy_tide.flux_from_tide(4.8e22, 1560.8e3, JUPITER, 6.711e8, 0.0094, 7.0e-4, 100.0)
    assert 1.0e11 < power < 1.0e12
    thickness_km = hydrosphere.ice_shell_thickness_km(flux, 100.0)
    assert 10.0 < thickness_km < 40.0
    assert flux == pytest.approx(7.0e-4 + power / (4.0 * math.pi * 1560.8e3 ** 2), rel=1e-9)


def test_enceladus_gets_a_tide_within_a_factor_of_ten_of_what_it_radiates():
    """Enceladus: 16 GW radiated at its eccentricity of 0.0047; the equilibrium shell of the model makes 1e11 W, six times more (the shell
    is 6 km against an observed 20 or so), the one body the stand-ins of `icy_tide` overshoot, within the factor of ten they claim."""
    from planetgen.physics import icy_tide
    _, power = icy_tide.flux_from_tide(1.08e20, 252.1e3, SATURN, 2.3802e8, 0.0047, 3.0e-4, 75.0)
    assert 1.6e9 < power < 1.6e11


def test_the_moon_is_flexed_solid_and_heated_by_a_fraction_of_a_percent_of_its_own_heat():
    """The Moon: k2 = 0.021 from its rigidity (observed 0.0242, [R4711]); on a cool mantle the Andrade Q of tens to hundreds gives the
    few times 1e8 to 1e10 W of the tide at e = 0.0549 (the lunar Q of 38 gives 6e8 W), and it stays solid."""
    mass, radius_km = 7.346e22, 1737.4
    assert tidal_heating.rigid_k2(mass, radius_km * 1000.0) == pytest.approx(0.0242, rel=0.2)
    forcing = tidal_heating.Forcing(EARTH, 3.84399e8, 0.0549)
    found = _heat(mass, radius_km, _rock(0.98, 0.02), forcing, COOL)
    assert 1.0e8 < found.tidal_w < 1.0e11
    assert found.mantle_temperature_k < interior_heat.mantle_of(mass, radius_km, _rock(0.98, 0.02)).heat_pipe_k
    assert found.tidal_w < 0.1 * found.own_w
    assert found.own_flux_w_m2 == pytest.approx(found.flux_w_m2, rel=0.1)


def test_io_is_melt_limited_and_radiates_a_watt_or_more_a_square_metre_with_no_number_fitted_to_it():
    """Io at e = 0.0041: the tide in the box (the Andrade rheology at the mantle's temperature) is far more than its own heat (4e11 W), so its melt
    erupts, it holds the temperature of the critical melt fraction and radiates the tide there. The power depends on the mantle's drawn viscosity
    (9e12 to 5e13 W over the corners of the drawn range) and Io's observed 0.9 to 1.5e14 W (2.24 W/m2) lies a factor of 2 to 3 above the top of it."""
    forcing = tidal_heating.Forcing(JUPITER, 4.218e8, 0.0041)
    rock = _rock(0.8, 0.2)
    pipe_k = interior_heat.mantle_of(8.93e22, 1821.6, rock).heat_pipe_k
    found = _heat(8.93e22, 1821.6, rock, forcing)
    assert found.mantle_temperature_k == pytest.approx(pipe_k)
    assert found.own_w < 0.05 * found.tidal_w
    assert found.flux_w_m2 == pytest.approx(found.total_w / (4 * math.pi * (1821.6e3) ** 2))
    powers = [_heat(8.93e22, 1821.6, rock, forcing, corner).tidal_w for corner in CORNERS]
    assert min(powers) < 1.0e13 < max(powers) and 0.9e14 / 4.0 < max(powers)
    assert all(1.0e12 < power < 1.0e15 for power in powers)
    assert all(_heat(8.93e22, 1821.6, rock, forcing, corner).mantle_temperature_k == pytest.approx(pipe_k, rel=0.01) for corner in CORNERS)   # GEN.1185: the softest corner sits 0.7 percent under the pipe


def test_a_weaker_tide_warms_the_mantle_without_melting_it_to_the_pipe():
    """Io at e = 3e-4 (a seventh of its tide): the tide is more than its own heat but the mantle conducts it away under the heat-pipe temperature,
    warmer than it would be with no tide, and the tide there is the rheology's at that temperature (the box is run with it)."""
    forcing = tidal_heating.Forcing(JUPITER, 4.218e8, 3.0e-4)
    rock = _rock(0.8, 0.2)
    pipe_k = interior_heat.mantle_of(8.93e22, 1821.6, rock).heat_pipe_k
    found = _heat(8.93e22, 1821.6, rock, forcing, COOL)
    alone = _heat(8.93e22, 1821.6, rock, None, COOL)
    assert alone.mantle_temperature_k < found.mantle_temperature_k < pipe_k
    assert found.tidal_w > alone.total_w * 0.5
    mantle = interior_heat.mantle_of(8.93e22, 1821.6, rock)
    assert found.tidal_w == pytest.approx(interior_heat._tide_w(forcing, mantle, found.mantle_temperature_k, COOL))


def test_the_tide_is_continuous_in_the_eccentricity_through_the_heat_pipe():
    """From a hundred times less tide to the full one the heat the mantle gives out never falls with the eccentricity by more than a third."""
    rock = _rock(0.8, 0.2)
    last = 0.0
    for eccentricity in (1.0e-4, 3.0e-4, 1.0e-3, 2.0e-3, 4.1e-3):
        found = _heat(8.93e22, 1821.6, rock, tidal_heating.Forcing(JUPITER, 4.218e8, eccentricity))
        assert found.total_w > 0.66 * last
        last = found.total_w


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
            if forcing is None:  # a planet or captured moon the tide has circularised has none left, and a stored 0.0 says so (GEN.778)
                assert body.eccentricity == 0.0 and body.tidal_flux_w_m2 == 0.0
                continue
            if (body.water_mass_fraction or 0.0) >= hydrosphere.ICY_MIN_WATER_FRACTION:
                icy += 1
                _, power = icy_tide.flux_from_tide(body.mass, body.radius * 1000.0, forcing.host_mass_kg, forcing.semi_major_axis_m, forcing.eccentricity,
                                                   body.internal_heat_flux_w_m2, interior_heat.icy_shell_top_k(body))
                assert body.tidal_flux_w_m2 == pytest.approx(power / (4.0 * math.pi * (body.radius * 1000.0) ** 2))
            else:
                rocky += 1
                assert body.tidal_flux_w_m2 >= 0.0
                area = 4 * math.pi * (body.radius * 1000.0) ** 2
                # never more than a fluid Love number over a Q of 1 would make (the box works on the dry radius, a few percent from the stored one)
                bound = tidal_heating.peale_power_w(1.5, forcing, body.radius * 1000.0) if forcing.eccentricity > 0.0 else 0.0
                if forcing.moons:      # and the moons' tides at the same ceiling (GEN.959): N |spin - n| with a k2/Q of 1.5
                    planet = spin_torques.Planet(body.mass, body.radius * 1000.0)
                    torques = spin_torques.Torques(planet, (), spin_torques.constant_dissipation(1.5, 1.0), None, spin_torques.GRAVITATIONAL_TORQUE_FACTOR)
                    for mass_kg, distance_m, retrograde, rate in forcing.moons:
                        bound += abs(torques.gravitational(forcing.spin_rate_rad_s, spin_torques.Raiser(mass_kg, distance_m, rate))) * abs(forcing.spin_rate_rad_s - rate)
                assert body.tidal_flux_w_m2 * area <= bound * 1.25
            # the own heat flow does not contain the tide: the habitability energy term adds it
            assert body.internal_heat_flux_w_m2 is not None and body.internal_heat_flux_w_m2 > 0.0
    assert seen > 10 and rocky > 0


def test_an_icy_body_on_an_orbit_so_round_the_tide_underflows_stores_zero():
    """An eccentricity of 1e-203 squares to 0.0, and the tide of 0.0 W/m2 is stored as that, as every rocky body's is (GEN.778)."""
    from planetgen.physics import hydrosphere
    for seed in range(12):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(SystemConfig(), galactic_center_dist_ly=26000.0)
        for body in _bodies(system):
            if body.body_type in ("g", "a") or (body.water_mass_fraction or 0.0) < hydrosphere.ICY_MIN_WATER_FRACTION:
                continue
            if tidal_heating.forcing_of(body) is None:
                continue
            body.eccentricity = 1.0e-203
            interior_heat.update_heat(body)
            assert body.tidal_flux_w_m2 == 0.0
            return
    pytest.skip("no icy body on an orbit in these seeds")


def test_the_solidus_runs_through_the_three_points_it_was_made_from():
    """1393 K at the surface, 2200 K at 10 GPa and about 4180 K at the core-mantle boundary, and Hirschmann's 1750 K at 3 GPa without being fitted to it."""
    assert tidal_heating.solidus_k(0.0) == pytest.approx(1393.0)
    assert tidal_heating.solidus_k(3.0) == pytest.approx(1750.0, rel=0.02)
    assert tidal_heating.solidus_k(10.0) == pytest.approx(2200.0, rel=0.03)
    assert tidal_heating.solidus_k(135.0) == pytest.approx(4180.0, rel=0.03)
    assert tidal_heating.liquidus_k(3.0) == tidal_heating.solidus_k(3.0) + tidal_heating.MELT_RANGE_K


def test_the_heat_pipe_temperature_is_the_critical_melt_fraction_of_the_range():
    from planetgen.physics import tidal_rheology
    solidus, liquidus = tidal_heating.solidus_k(4.0), tidal_heating.liquidus_k(4.0)
    pipe = tidal_heating.heat_pipe_temperature_k(4.0)
    assert tidal_rheology.melt_fraction(pipe, solidus, liquidus) == pytest.approx(tidal_rheology.CRITICAL_MELT_FRACTION)


def test_the_tide_temperature_depends_on_the_content_and_the_draw_alone():
    from tests.tide_support import rocky
    body = rocky(7.342e22, 1737.4)
    assert interior_heat.tide_temperature_k(body) == interior_heat.tide_temperature_k(rocky(7.342e22, 1737.4))
    assert interior_heat.tide_temperature_k(body) != interior_heat.tide_temperature_k(rocky(7.342e22, 1737.4, eta0_pa_s=9.3e7, activation_kj_mol=260.0))


def test_a_rocky_body_without_its_draws_has_no_tidal_response():
    from tests.tide_support import rocky
    bare = rocky(7.342e22, 1737.4, mantle_eta0_pa_s=None)
    with pytest.raises(interior_heat.InteriorHeatError):
        interior_heat.tidal_response(bare, 1.0e-6)


def test_the_moon_stays_solid_at_the_pipe_temperature_and_io_does_not():
    """The heat-pipe temperature lies above a Moon's cool mantle and the mantle of an Io is held at it (one body's draw is not the other's)."""
    moon_rock = _rock(0.98, 0.02)
    pipe = interior_heat.mantle_of(7.346e22, 1737.4, moon_rock).heat_pipe_k
    assert _heat(7.346e22, 1737.4, moon_rock, tidal_heating.Forcing(EARTH, 3.84399e8, 0.0549), COOL).mantle_temperature_k < pipe
