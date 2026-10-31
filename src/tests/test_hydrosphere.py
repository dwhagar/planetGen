"""
Hydrosphere and ocean chemistry (GEN.88, `planetgen.physics.hydrosphere`):
the design doc's section 5 numbers, the ocean-class rules and what every
generated body stores.
"""

import types

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants, hydrosphere as hy, tidal_heating, water_layer as wl
from planetgen.util import draw
from planetgen.util.random import log_uniform

from tests.air_support import air, fold_air
from tests.fuzz_support import deterministic_entropy


def _generate_hydrosphere(body, internal_flux_w_m2):
    """`generate_hydrosphere` after the surface water is drawn, as `planets.calculate_atmospheric_conditions` does."""
    return hy.generate_hydrosphere(body, internal_flux_w_m2, hy.draw_surface_water(body) if body.body_type != "g" else 0.0)


EARTH = constants.EARTH_MASS_TO_KG


@pytest.mark.parametrize("bottom_k, pressure_mpa", [
    (251.2, 210.0), (256.2, 350.0), (272.1, 620.0), (286.8, 800.0), (300.4, 1000.0), (327.0, 1500.0), (353.8, 2200.0),
])
def test_the_liquidus_follows_the_recomputed_melting_curves(bottom_k, pressure_mpa):
    assert wl.liquidus_pressure_pa(bottom_k) == pytest.approx(pressure_mpa * 1e6, rel=1e-9)


def test_the_deepest_liquid_at_earth_gravity_is_about_seventy_kilometres():
    """The design doc 5.2 had 73 km at a fixed 1050 kg/m3; the water is lighter than that over the first 700 MPa."""
    deepest = wl.column(EARTH, 6371.0, 0.02 * EARTH, bottom_temperature_k=280.0)
    assert deepest.liquid_km == pytest.approx(65.5, rel=0.03)
    assert deepest.hp_ice_km > 0 and deepest.base_pressure_pa > wl.liquidus_pressure_pa(280.0)


def test_earths_water_covers_seven_tenths():
    layer = hy.water_layer_depth_km(EARTH, 6371.0, hy.EARTH_OCEAN_FRACTION)
    assert layer == pytest.approx(2.74, rel=0.03)
    assert hy.ocean_fraction(layer, constants.EARTH_GRAVITY) == pytest.approx(0.72, abs=0.02)
    assert hy.ocean_fraction(hy.water_layer_depth_km(EARTH, 6371.0, 2e-3), constants.EARTH_GRAVITY) == 1.0


def test_europas_lid_and_tides():
    assert hy.ice_shell_thickness_km(0.05, 100.0) == pytest.approx(11.4, abs=0.2)
    assert hy.pressure_melted_shell_km(0.05, 100.0, 1.3) < hy.ice_shell_thickness_km(0.05, 100.0)
    europa = tidal_heating.peale_power_w(3e-4, tidal_heating.Forcing(constants.JUPITER_MASS_TO_KG, 6.709e8, 0.009), 1560.8e3)
    assert europa == pytest.approx(1.4e11, rel=0.15)


def test_the_lid_counts_by_mass_and_high_pressure_ice_seals_the_ocean():
    # A lid of 10 km of ice is 9.17e6 kg/m2: a column lighter than that freezes solid.
    area = 4 * 3.141592653589793 * (6371.0e3) ** 2
    light = wl.column(EARTH, 6371.0, 9.0e6 * area, shell_km=10.0, bottom_temperature_k=273.15)
    assert light.liquid_km == 0.0 and light.hp_ice_km == 0.0 and light.ice_shell_km == pytest.approx(9.0e6 / 917.0 / 1000.0)
    open_ocean = wl.column(EARTH, 6371.0, 2.0e7 * area, shell_km=10.0, bottom_temperature_k=273.15)
    assert open_ocean.ice_shell_km == 10.0 and open_ocean.liquid_km == pytest.approx(10.34, rel=0.02) and open_ocean.hp_ice_km == 0.0
    deep = wl.column(EARTH, 6371.0, 0.05 * EARTH, bottom_temperature_k=300.0)
    assert deep.liquid_km == pytest.approx(88.1, rel=0.03) and deep.hp_ice_km > 0


def _ocean(**gases):
    return types.SimpleNamespace(mantle_redox=gases.pop("redox", "intermediate"), **air(**gases))


def test_ocean_classes_follow_the_rules():
    assert hy.ocean_class(_ocean(), 5.0, 1e-3, 1.0) == "ice-sealed"
    assert hy.ocean_class(_ocean(), 0.0, 1e-5, 0.1) == "chloride brine"
    acid = _ocean(redox="oxidized", p_so2_kpa=1.0, p_co2_kpa=10.0)
    assert hy.ocean_class(acid, 0.0, 1e-3, 0.5) == "acid sulfate"
    soda = _ocean(redox="reduced", p_so2_kpa=0.0, p_co2_kpa=50.0)
    assert hy.ocean_class(soda, 0.0, 1e-3, 0.95) == "soda"
    assert hy.ocean_class(soda, 0.0, 1e-3, 0.5) == "neutral"


# The water mass fraction range a test body of each kind of world draws from (log-uniform): an Earth-like world, a dry rock, a snowball,
# an ocean world and an icy moon.
TEST_WATER_MASS_FRACTION = {"M": (1e-4, 5e-4), "C": (1e-6, 1e-4), "P": (1e-4, 1e-3), "O": (1e-2, 0.1), "D": (0.05, 0.5)}


def _body(planet_class, temperature, pressure_pa, mass_earth=1.0, radius_km=6371.0, **extra):
    fields = dict(body_type="t", planet_class=planet_class, mass=mass_earth * EARTH, radius=radius_km,
                  water_mass_fraction=log_uniform(*TEST_WATER_MASS_FRACTION[planet_class]),
                  gravity=mass_earth / (radius_km / 6371.0) ** 2, surface_temperature=temperature,
                  atmospheric_pressure=pressure_pa, is_moon=False, mantle_redox="intermediate",
                  p_co2_kpa=0.04, p_so2_kpa=0.0, p_h2_kpa=0.0)
    fields.update(extra)
    return types.SimpleNamespace(**fold_air(fields))


def test_states_by_temperature_and_pressure():
    draw.set_run_seed(88)
    earth = _body("M", 288.0, 101325.0)
    _generate_hydrosphere(earth, 0.087)
    assert earth.hydrosphere == "surface ocean" and earth.ocean_class in hy.OCEAN_CLASSES
    assert earth.ocean_fraction + earth.land_fraction == pytest.approx(1.0)
    assert 6.5 <= earth.ocean_ph <= 8.5 or earth.ocean_class != "neutral"
    venus = _body("M", 740.0, 9.2e6)
    _generate_hydrosphere(venus, 0.087)
    assert venus.hydrosphere == "vapour" and venus.ocean_depth_km is None
    moon = _body("C", 250.0, 0.0)
    _generate_hydrosphere(moon, 0.02)
    assert moon.hydrosphere == "dry" and moon.land_fraction == 1.0
    snowball = _body("P", 200.0, 50000.0)
    _generate_hydrosphere(snowball, 0.087)
    assert snowball.hydrosphere in ("ice", "ice-covered ocean") and snowball.ice_shell_km > 0
    hycean = _body("O", 300.0, 1e6, p_h2_kpa=900.0)
    _generate_hydrosphere(hycean, 0.087)
    assert hycean.hydrosphere == "hycean"
    giant = types.SimpleNamespace(body_type="g", water_mass_fraction=0.3)
    _generate_hydrosphere(giant, 1.0)
    assert giant.water_mass_fraction == 0.3
    assert all(getattr(giant, field) is None for field in hy.HYDROSPHERE_FIELDS if field != "water_mass_fraction")


def test_a_liquid_surface_of_chloride_brine_is_a_brine_seep_not_a_sea():
    """GEN.908: surface water under `BRINE_MAX_WATER_FRACTION` of the mass makes a concentrated brine, a puddle and not an ocean."""
    draw.set_run_seed(908)
    for temperature, pressure in ((288.0, 101325.0), (250.0, 30000.0), (265.0, 600.0)):
        seep = _body("M", temperature, pressure, water_mass_fraction=3e-6)
        _generate_hydrosphere(seep, 0.087)
        if seep.ocean_depth_km:
            assert seep.ocean_class == "chloride brine" and seep.hydrosphere == hy.BRINE_SEEP == "brine seep"
    wet = _body("M", 288.0, 101325.0, water_mass_fraction=3e-4)
    _generate_hydrosphere(wet, 0.087)
    assert wet.hydrosphere == "surface ocean" and wet.ocean_class != "chloride brine"
    warm_seep = _body("M", 288.0, 101325.0, water_mass_fraction=3e-6)
    _generate_hydrosphere(warm_seep, 0.087)
    assert warm_seep.hydrosphere == "brine seep" and warm_seep.ocean_depth_km


def test_a_brine_seep_is_liquid_water_but_not_a_sunlit_sea_nor_a_wet_world_nor_an_induction_ocean():
    from planetgen.physics import biogenic_oxygen, habitability_world, magnetism, planet_classes
    assert "brine seep" not in biogenic_oxygen.LIVING_HYDROSPHERES
    assert "brine seep" not in magnetism.INDUCTION_HYDROSPHERES
    seep = types.SimpleNamespace(hydrosphere="brine seep")
    assert magnetism.mobile_lid_chance(seep, 1.0) == 0.0
    values = {"hydrosphere": "brine seep", "surface_temperature_k": 250.0, "water_mass_fraction": 3e-6}
    assert habitability_world.solvent(values) == "water" and habitability_world.water_source(values) == "liquid"


def test_a_brine_seep_has_no_light_class_but_keeps_the_chemical_ones():
    """GEN.908: the sunlit classes (1a, 1b, 1c, 2) need an open sea; the chemical ones need liquid water over rock."""
    from planetgen.physics import life_classes
    from tests.test_life_classes import world
    scores = life_classes.score_classes(world(hydrosphere="brine seep"))
    assert all(scores[c].score == 0 for c in ("1a", "1b", "1c", "2"))
    assert scores["3a"].score > 0 and scores["4b"].score > 0


def test_ocean_worlds_are_mostly_ice_sealed():
    draw.set_run_seed(9)
    worlds = [_body("O", 290.0, 101325.0) for _ in range(100)]
    for world in worlds:
        _generate_hydrosphere(world, 0.087)
    assert all(world.ocean_fraction == 1.0 for world in worlds)
    sealed = [world for world in worlds if world.ocean_class == "ice-sealed"]
    assert len(sealed) > 50
    assert all(world.hp_ice_km > 0 and world.phosphorus == "starved" for world in sealed)


def test_generated_bodies_store_a_consistent_hydrosphere():
    seen = set()
    for seed in range(40):
        with deterministic_entropy(seed):
            system = StarSystem(system_config=SystemConfig())
        for planet in system.planets:
            if not hasattr(planet, "moons"):
                continue
            for body in [planet, *planet.moons]:
                if body.body_type == "g":
                    assert all(getattr(body, field) is None for field in hy.HYDROSPHERE_FIELDS
                               if field != "water_mass_fraction")
                    continue
                seen.add(body.hydrosphere)
                assert body.hydrosphere in hy.HYDROSPHERE_STATES
                assert body.water_mass_fraction >= 0  # 0: its early runaway greenhouse took it all (GEN.639), or a disk moon inside its ice line was anhydrous from the start (GEN.466)
                if body.water_mass_fraction == 0:
                    assert body.hydrosphere == "dry"
                assert body.ocean_fraction + body.land_fraction == pytest.approx(1.0)
                liquid = body.hydrosphere in ("ice-covered ocean", "surface ocean", "brine seep", "hycean")
                assert (body.ocean_depth_km is not None) == liquid == (body.ocean_class is not None)
                if liquid:
                    (ph_lo, ph_hi), (aw_lo, aw_hi) = hy.OCEAN_CHEMISTRY[body.ocean_class]
                    assert ph_lo <= body.ocean_ph <= ph_hi and aw_lo <= body.water_activity <= aw_hi
                    assert body.phosphorus in ("high", "limited", "starved")
    # a runaway-greenhouse vapour world is rare (2 in 590 bodies) and is checked by hand above, so a sample need not hold one
    assert {"ice", "dry"} <= seen


def _moon(planet_class, mass_earth, radius_km, distance_au, **extra):
    return _body(planet_class, 100.0, 0.0, mass_earth=mass_earth, radius_km=radius_km, is_moon=True,
                 primary_mass_kg=constants.JUPITER_MASS_TO_KG, distance=distance_au, **extra)


def test_the_tide_stored_with_the_heat_model_thickens_the_melt_of_a_lid():
    """GEN.737: the hydrosphere reads the tide the heat model stored (`tidal_flux_w_m2`) and adds it to the body's own heat flow; it no longer
    draws an eccentricity or a resonance boost of its own."""
    for seed in range(8):
        draw.set_run_seed(seed)
        europa = _moon("D", 8.0e-3, 1560.8, 6.709e8 / constants.AU_TO_M, tidal_flux_w_m2=4.6e-3)
        _generate_hydrosphere(europa, 0.002)
        assert europa.tidal_flux_w_m2 == 4.6e-3
        draw.set_run_seed(seed)
        cold = _moon("D", 8.0e-3, 1560.8, 6.709e8 / constants.AU_TO_M)
        _generate_hydrosphere(cold, 0.002)
        assert not hasattr(cold, "tidal_flux_w_m2")
        assert (europa.ice_shell_km or 0.0) <= (cold.ice_shell_km or 0.0)  # more heat from below, a thinner lid


def test_the_tide_is_in_the_chemical_energy_of_an_ocean_under_ice():
    from planetgen.physics import habitability_world as hw
    body = types.SimpleNamespace(body_type="t", surface_temperature=100.0, tidal_flux_w_m2=0.004)
    values = {"hydrosphere": "ice-covered ocean"}
    heat = 0.002 + body.tidal_flux_w_m2
    assert hw.energy_flux_w_m2(values, None, None, heat) == pytest.approx(hw.CHEMICAL_SHARE * 0.006)
    assert hw.energy_flux_w_m2(values, None, None, heat) > 2.5 * hw.energy_flux_w_m2(values, None, None, 0.002)


def test_gen833_a_small_icy_moon_that_holds_a_thick_ice_shell_is_covered_by_it():
    """GEN.833: a cold moon of 200 km radius with a fifth of its mass as water is ice through its outer layers; the ground draw
    reads its cover and it is all of the surface, not the toy cover of the water a liquid ocean would have had."""
    for seed in range(8):
        draw.set_run_seed(seed)
        moon = _body("D", 70.0, 0.0, mass_earth=3.0e-6, radius_km=200.0, is_moon=True, water_mass_fraction=0.2,
                     primary_mass_kg=constants.JUPITER_MASS_TO_KG, distance=0.01)
        _generate_hydrosphere(moon, 0.002)
        assert moon.hydrosphere in hy.FROZEN_COVER_STATES
        assert moon.ocean_fraction == 1.0 and moon.land_fraction == 0.0


def test_gen833_the_cover_follows_the_water_the_body_holds():
    """GEN.833: a frozen body with little water is covered in part, and one with more water is covered more."""
    covers = []
    for water in (1e-5, 1e-4, 1e-3, 1e-2, 0.2):
        body = types.SimpleNamespace(mass=3.0e-6 * EARTH, radius=200.0, water_mass_fraction=water)
        gravity = 3.0e-6 / (200.0 / 6371.0) ** 2 * constants.EARTH_GRAVITY
        covers.append(hy.frozen_cover(body, gravity))
    assert covers == sorted(covers) and covers[0] < 0.5 and covers[-1] == 1.0
    assert hy.frozen_cover(types.SimpleNamespace(mass=EARTH, radius=6371.0, water_mass_fraction=0.0), constants.EARTH_GRAVITY) == 0.0


def test_gen833_below_earths_gravity_ice_relief_does_not_grow_but_a_rock_basins_does():
    layer_km = 5.0
    small = 0.05 * constants.EARTH_GRAVITY
    assert hy.ocean_fraction(layer_km, small, icy=True) == pytest.approx(min(1.0, (2 * layer_km / hy.EARTH_RELIEF_KM) ** 0.5))
    assert hy.ocean_fraction(layer_km, small, icy=True) > hy.ocean_fraction(layer_km, small)
    heavy = 3.0 * constants.EARTH_GRAVITY
    assert hy.ocean_fraction(layer_km, heavy, icy=True) == hy.ocean_fraction(layer_km, heavy)


def test_gen1222_a_thick_cold_lid_does_not_melt_its_own_base_under_the_surface_temperature():
    """GEN.1222: the pressure of a lid lowers its base melting point, but not past the ice III triple point; with no floor a cold, low-heat
    body's thick lid pulled its base under 25 K, the shell came out 0 and the body a liquid ocean at the surface."""
    gravity = 0.1
    for flux in (1e-6, 1e-5, 1e-4, 1e-3):
        shell = hy.pressure_melted_shell_km(flux, 25.0, gravity)
        assert shell > 0.0
        assert shell >= hy.ice_shell_thickness_km(flux, 25.0, 251.2)
    # the Earth-like case stays where it was: a warm top under a thin lid is untouched by the floor
    assert hy.pressure_melted_shell_km(0.06, 273.0 - 1.0, 9.8) == pytest.approx(hy.ice_shell_thickness_km(0.06, 272.0, hy.WATER_MELTING_K), rel=0.2)


def test_gen1222_a_warm_world_under_a_crushing_air_has_no_surface_ocean():
    """GEN.1222: at 278 K under 4.5 GPa the water is high-pressure ice all the way down, so the state is `ice` and no ocean column is stored (it was
    `surface ocean` with no ocean depth, chemistry or pH)."""
    draw.set_run_seed(1222)
    crushed = _body("O", 278.0, 4.5e9, mass_earth=3.6, radius_km=9000.0)
    _generate_hydrosphere(crushed, 0.087)
    assert crushed.hydrosphere == "ice" and crushed.ocean_depth_km is None and crushed.ocean_class is None
    assert crushed.ice_shell_km and crushed.ice_shell_km == crushed.hp_ice_km
