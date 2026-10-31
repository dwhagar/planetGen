# tests/test_giant_heat.py

"""
GEN.810: the internal heat of a giant (`planetgen.physics.giant_heat`) reads the stored interior layers, the helium depletion and the
helium shell: the base is Jupiter's scaled by the mass, a giant whose helium rained out radiates more, a giant whose helium went to a
shell radiates more than one whose helium stayed in a gradient, the potential of the layers is the exact one for a uniform sphere, a
giant whose helium did not rain adds nothing, and the heat is stored with the body (`interior_heat.update_heat`).
"""

import math
from types import SimpleNamespace

import pytest

from planetgen.physics import constants, giant_heat as gh, giant_interior as gi, interior_heat, rogue_surface
from planetgen.util import draw

ICES = {"water_liquid": 0.25, "nh3_ice": 0.02, "ch4_ice": 0.03}
SATURN = dict(mass_earth=95.2, radius_km=60268.0, heavy_mass_earth=25.0, envelope_z=0.04, helium_share=0.2745,
              temperature_1bar_k=140.0, ice_shares=ICES)
JUPITER = dict(mass_earth=317.8, radius_km=71492.0, heavy_mass_earth=42.0, envelope_z=0.04, helium_share=0.2745,
               temperature_1bar_k=166.0, ice_shares=ICES)
NEPTUNE = dict(mass_earth=17.1, radius_km=24764.0, heavy_mass_earth=14.6, envelope_z=0.1, helium_share=0.2745,
               temperature_1bar_k=72.0, ice_shares={"water_liquid": 0.6, "nh3_ice": 0.1, "ch4_ice": 0.1})
KG = constants.EARTH_MASS_TO_KG


def interior(spec, seed=1, **changes):
    with draw.bound(seed):
        return gi.draw_interior(gi.GiantSubject(**{**spec, **changes}))


def body_of(spec, seed=1, age=4.6, **changes):
    """A stand-in body holding what `giant_interior.body_values` stores."""
    found = interior(spec, seed, **changes)
    body = SimpleNamespace(radius=spec["radius_km"], star=SimpleNamespace(age=age))
    for name, value in gi.body_values(found).items():
        setattr(body, name, value)
    return body, found


def test_the_base_is_the_cooling_law_of_the_rogue_giants_and_falls_with_the_age():
    radius_m = JUPITER["radius_km"] * 1e3
    young, old = (gh.base_flux_w_m2(gi.EARTH_MASS_JUPITER, radius_m, age) for age in (1.0, 4.6))
    assert young > old > 0.0
    assert old == pytest.approx(rogue_surface.giant_internal_flux_w_m2(gi.EARTH_MASS_JUPITER * constants.EARTH_MASS_TO_KG, JUPITER["radius_km"], 4.6))


def test_the_potential_of_a_uniform_sphere_at_its_half_mass_radius_is_the_exact_one():
    mass, radius = 10.0 * KG, 8.0e7
    layer = gi.Layer("molecular_envelope", "molecular_fluid", 1e-4, 1.0, 10.0, {"hydrogen": 10.0}, 1.0, 0.0)
    body = gi.GiantInterior((layer,), 10.0, "gas_giant", 100.0, 1.0)
    r = radius * 0.5 ** (1.0 / 3.0)
    exact = -constants.G * mass / (2.0 * radius ** 3) * (3.0 * radius ** 2 - r ** 2)
    assert gh.potentials(body, radius)[0] == pytest.approx(exact, rel=1e-12)


def test_the_potential_is_deeper_in_the_lower_layers():
    phi = gh.potentials(interior(SATURN), SATURN["radius_km"] * 1e3)
    assert all(a > b for a, b in zip(phi, phi[1:])) and all(value < 0.0 for value in phi)


def test_a_giant_whose_helium_did_not_rain_adds_nothing_and_its_heat_is_the_base():
    for spec in (JUPITER, NEPTUNE):
        found = interior(spec, temperature_1bar_k=300.0) if spec is JUPITER else interior(spec)
        assert gh.rain_energy_j(found, spec["radius_km"] * 1e3) == (0.0, 0.0)
        body = SimpleNamespace(radius=spec["radius_km"], star=SimpleNamespace(age=4.6), helium_depletion=found.helium_depletion,
                               has_helium_shell=found.has_helium_shell, **{k: v for k, v in gi.body_values(found).items()
                                                                           if k not in ("helium_depletion", "has_helium_shell")})
        heat = gh.of(body)
        assert heat.rain_w_m2 == 0.0 and heat.rain_power_w == 0.0
        assert heat.flux_w_m2 == pytest.approx(gh.base_flux_w_m2(found.total_mass_earth, spec["radius_km"] * 1e3, 4.6))


def test_a_saturn_like_giant_with_a_helium_shell_radiates_more_than_one_without():
    body, found = body_of(SATURN)
    assert found.has_helium_shell and body.has_helium_shell and body.helium_depletion < 0.5
    warm, warm_found = body_of(SATURN, temperature_1bar_k=180.0)
    assert not warm_found.has_helium_shell and warm.helium_depletion == 1.0
    cold, with_shell = gh.of(body), gh.of(warm)
    assert cold.flux_w_m2 > with_shell.flux_w_m2 and cold.rain_w_m2 > 0.0
    assert with_shell.rain_w_m2 == 0.0 and with_shell.flux_w_m2 == pytest.approx(cold.base_w_m2)
    assert cold.flux_w_m2 == pytest.approx(cold.base_w_m2 + cold.rain_w_m2)


def test_more_depletion_sinks_more_helium_and_gives_more_heat():
    fluxes = []
    for temperature in (165.0, 150.0, 140.0):          # the rain starts under 166 K and is complete by 140 K
        body, _ = body_of(SATURN, temperature_1bar_k=temperature)
        fluxes.append(gh.of(body).rain_w_m2)
    assert fluxes[0] < fluxes[1] < fluxes[2]


def test_the_rains_heat_is_the_average_over_the_age_and_the_released_share_of_the_energy():
    body, found = body_of(SATURN, age=4.6)
    moved, energy = gh.rain_energy_j(found, SATURN["radius_km"] * 1e3)
    heat = gh.of(body)
    assert moved == pytest.approx(heat.moved_helium_earth * KG) and energy > 0.0
    assert heat.rain_power_w == pytest.approx(gh.RAIN_HEAT_SHARE * energy / (4.6 * gh.SECONDS_PER_GYR))
    older, _ = body_of(SATURN, age=9.2)
    assert gh.of(older).rain_power_w == pytest.approx(0.5 * heat.rain_power_w)
    area = 4.0 * math.pi * (SATURN["radius_km"] * 1e3) ** 2
    assert heat.rain_w_m2 == pytest.approx(heat.rain_power_w / area)


def test_what_sank_is_the_helium_the_lower_layers_hold_above_the_natal_share():
    found = interior(SATURN)
    moves = gh.helium_moved(found)
    assert math.fsum(moves) == pytest.approx(0.0, abs=1e-9 * found.total_mass_earth * KG)
    assert moves[0] < 0.0 < moves[-1] or any(move > 0.0 for move in moves[1:])
    assert gh.rain_energy_j(found, SATURN["radius_km"] * 1e3)[0] == pytest.approx(math.fsum(m for m in moves if m > 0.0))


def test_a_saturn_like_giants_heat_is_the_order_of_what_was_measured():
    """Saturn radiates 1.88 W/m2 of its own (T_int 75.9 K, Li et al. 2018); the model is not fitted to it, and lands within a factor of 2."""
    heat = gh.of(body_of(SATURN)[0])
    assert 0.9 < heat.flux_w_m2 < 3.8


def test_update_heat_stores_a_giants_heat_and_leaves_a_rocky_body_to_its_own_history():
    body, _ = body_of(SATURN)
    body.mass = SATURN["mass_earth"] * KG
    for name in interior_heat.OUTPUT_FIELDS:
        setattr(body, name, 1.0)
    assert interior_heat.update_heat(body) is None
    assert body.internal_heat_flux_w_m2 == pytest.approx(gh.of(body).flux_w_m2)
    assert body.urey_ratio is None and body.tidal_flux_w_m2 is None
    assert interior_heat.update_heat(body) is None and body.internal_heat_flux_w_m2 == pytest.approx(gh.of(body).flux_w_m2)


def test_a_body_with_no_stored_interior_has_no_giant_heat():
    assert gh.of(SimpleNamespace(radius=1000.0, star=SimpleNamespace(age=4.6))) is None
