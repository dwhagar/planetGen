# tests/test_body_strength.py

"""
The roundness ratio, yield strength and binding energy a body stores (GEN.1159, GEN.1178): the module `physics/body_strength.py`,
the fill on a planet, a moon and a comet, the draw from the body's own stream, and the columns in a saved system.
"""

import types

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import binding_energy, body_strength, roundness
from planetgen.util import draw

EARTH = dict(mass=5.972e24, radius=6371.0, density=5.51, water_mass_fraction=0.0002, ice_mass_fraction=0.0, metal_mass_fraction=0.32,
             rock_mass_fraction=0.68, body_type="t", surface_temperature=288.0, distance=1.0, is_moon=False)
MOON = dict(mass=7.3e22, radius=1737.0, density=3.34, water_mass_fraction=0.0, ice_mass_fraction=0.0, metal_mass_fraction=0.02,
            rock_mass_fraction=0.98, body_type="t", surface_temperature=250.0, distance=0.0026, is_moon=True)
SMALL = dict(mass=2.1e19, radius=135.0, density=1.2, water_mass_fraction=0.3, ice_mass_fraction=0.0, metal_mass_fraction=0.0,
             rock_mass_fraction=0.7, body_type="t", surface_temperature=100.0, distance=0.0026, is_moon=True)
GIANT = dict(mass=1.898e27, radius=69911.0, density=1.33, body_type="g", surface_temperature=130.0, distance=5.2, is_moon=False)


def _body(**fields):
    return types.SimpleNamespace(**fields)


def test_a_planet_is_round_with_a_strength_inside_its_range_and_a_binding_energy_over_the_uniform_ball():
    out = body_strength.values(_body(**EARTH))
    assert out["roundness_ratio"] > 20.0
    low, high = roundness.RANGES_PA["rock"]
    assert low * 0.5 < out["yield_strength_pa"] < high * 2.0
    uniform = binding_energy.uniform_j(EARTH["mass"], EARTH["radius"] * 1000.0)
    assert uniform < out["binding_energy_j"] < 1.3 * uniform       # a core makes it more tightly bound than a uniform ball


def test_a_small_icy_moon_is_not_round_and_a_large_rocky_one_is():
    assert body_strength.values(_body(**SMALL))["roundness_ratio"] < 1.0
    assert body_strength.values(_body(**MOON))["roundness_ratio"] > 1.0


def test_a_gas_giant_keeps_only_its_binding_energy():
    out = body_strength.values(_body(**GIANT))
    assert out["roundness_ratio"] is None and out["yield_strength_pa"] is None
    assert out["binding_energy_j"] == pytest.approx(binding_energy.binding_energy_j("gas-giant", GIANT["mass"], GIANT["radius"] * 1000.0))


def test_a_body_without_a_mass_or_a_radius_has_none_of_the_three():
    assert body_strength.values(_body(**{**EARTH, "mass": None})) == dict.fromkeys(body_strength.FIELDS)
    assert body_strength.values(_body(**{**EARTH, "radius": 0.0})) == dict.fromkeys(body_strength.FIELDS)


def test_the_quantile_comes_from_the_bodys_own_stream_and_leaves_the_run_stream_alone():
    body = _body(**MOON)
    draw.set_run_seed(5)
    first = body_strength.quantile_of(body)
    after_first = draw.random()
    draw.set_run_seed(5)
    assert body_strength.quantile_of(body) == first
    assert draw.random() == after_first                       # the run's stream is where it was
    assert 0.0 <= first < 1.0
    assert body_strength.quantile_of(_body(**{**MOON, "mass": 7.4e22})) != first


def test_the_strength_follows_the_quantile_between_the_ends_of_its_range():
    weak, strong = (body_strength.values(_body(**MOON), quantile=q)["yield_strength_pa"] for q in (0.0, 1.0))
    low, high = roundness.RANGES_PA["rock"]
    assert weak < strong
    assert low * 0.9 < weak < roundness.ROCK_STRENGTH_PA < strong < high * 1.1


def test_the_strength_is_the_cold_strength_whatever_the_surface_temperature():
    cold, hot = (body_strength.values(_body(**{**MOON, "surface_temperature": t}), quantile=0.5) for t in (30.0, 700.0))
    assert cold["yield_strength_pa"] == hot["yield_strength_pa"]


def test_a_comet_has_a_density_an_ice_share_and_the_three_values():
    from planetgen.generation.comet import Comet
    draw.set_run_seed(11)
    comet = Comet(SystemConfig(), 1.0, orbit_type="elliptical")
    assert comet.density_g_cm3 == pytest.approx(comet.mass_kg / (4.0 / 3.0 * 3.141592653589793 * (comet.nucleus_diameter_km * 500.0) ** 3) / 1000.0, rel=1e-9)
    assert 0.0 <= comet.ice_mass_fraction <= 1.0
    assert comet.roundness_ratio is not None and comet.yield_strength_pa > 0.0 and comet.binding_energy_j > 0.0
    again = Comet.from_dict(comet.to_dict(), SystemConfig())
    assert again.roundness_ratio == comet.roundness_ratio and again.binding_energy_j == comet.binding_energy_j


def test_a_generated_system_stores_the_values_wherever_the_body_has_them(mysql_config):
    config = SystemConfig()
    config.PLANETS = True
    config.MAX_PLANETS = True
    seen = {"planets": 0, "moons": 0, "comets": 0, "asteroid_belts": 0}
    for seed in range(1000, 1012):
        draw.set_run_seed(seed)
        saved = store.save_system(StarSystem(system_config=config, galactic_center_dist_ly=26000.0), config, config=mysql_config)
        with store.get_connection(mysql_config) as conn:
            for table in ("planets", "moons"):
                for row in conn.execute(f"SELECT body_type, planet_class, roundness_ratio, yield_strength_pa, binding_energy_j, density_g_cm3 FROM {table} WHERE star_system_id = ?"
                                        if table == "planets" else
                                        "SELECT m.planet_class, m.roundness_ratio, m.yield_strength_pa, m.binding_energy_j, m.density_g_cm3 FROM moons m JOIN planets p ON p.id = m.planet_id WHERE p.star_system_id = ?",
                                        (saved,)).fetchall():
                    seen[table] += 1
                    assert row["binding_energy_j"] is not None and row["binding_energy_j"] > 0.0
                    giant = row["planet_class"] in ("I", "J", "T")
                    assert (row["roundness_ratio"] is None) == giant
                    assert (row["yield_strength_pa"] is None) == giant
            for row in conn.execute("SELECT density_g_cm3, ice_mass_fraction, roundness_ratio, yield_strength_pa, binding_energy_j FROM comets WHERE star_system_id = ?", (saved,)).fetchall():
                seen["comets"] += 1
                assert all(value is not None for value in row.values())
            for row in conn.execute("SELECT density_g_cm3, ice_mass_fraction FROM asteroid_belts WHERE star_system_id = ?", (saved,)).fetchall():
                seen["asteroid_belts"] += 1
                assert row["density_g_cm3"] is not None and 0.0 <= row["ice_mass_fraction"] <= 1.0
    assert all(seen.values()), seen
