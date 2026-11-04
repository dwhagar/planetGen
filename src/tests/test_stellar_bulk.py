# tests/test_stellar_bulk.py

"""
GEN.1162: the mean density and binding energy a star, a white dwarf, a neutron star, a black hole, a rogue planet and an interstellar comet store, and the roundness ratio and
strength of a terrestrial rogue planet and of an interstellar comet (`physics/body_strength.py`).
"""

import types

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
from planetgen.generation.phenomena.rogue import InterstellarComet, RoguePlanet
from planetgen.physics import binding_energy, body_strength, constants
from planetgen.util import draw

MSUN = constants.SOLAR_MASS_TO_KG


def _star(mass_solar, radius_km, yerkes):
    return types.SimpleNamespace(mass=mass_solar * MSUN, radius=radius_km, yerkes_class=yerkes)


def test_the_kind_of_a_star_row_follows_its_class_and_mass():
    assert body_strength.stellar_kind("V", MSUN) == "star"
    assert body_strength.stellar_kind("III", 2 * MSUN) == "star"
    assert body_strength.stellar_kind("VII", 0.6 * MSUN) == "white-dwarf"
    assert body_strength.stellar_kind("V", 0.05 * MSUN) == "brown-dwarf"
    assert body_strength.stellar_kind("BH", 10 * MSUN) == "black-hole"
    assert body_strength.stellar_kind("NS", 1.4 * MSUN) == "neutron-star"


def test_the_sun_has_its_density_and_a_binding_energy_of_the_standard_model():
    density, energy = body_strength.star_values(_star(1.0, 696340.0, "V"))
    assert density == pytest.approx(1.41, rel=0.02)
    assert energy == pytest.approx(5.7e41, rel=0.05)          # 3/2 G M^2 / R, the n = 3 polytrope


def test_a_white_dwarf_is_denser_and_more_tightly_bound_than_a_neutron_free_star():
    sun_density, sun_energy = body_strength.star_values(_star(1.0, 696340.0, "V"))
    density, energy = body_strength.star_values(_star(0.6, 8700.0, "VII"))
    assert density > 1.0e5 * sun_density and energy > 10.0 * sun_energy


def test_a_black_hole_stores_minus_one_for_an_infinite_binding_energy():
    density, energy = body_strength.star_values(_star(10.0, 29.5, "BH"))
    assert energy == binding_energy.INFINITE == -1.0 and density > 1.0e13
    assert body_strength.remnant_values("black-hole", 10.0, 29.5) == (density, -1.0)


def test_a_neutron_star_has_about_a_ninth_of_its_rest_mass_as_binding_energy():
    density, energy = body_strength.remnant_values("neutron-star", 1.4, 12.0)
    assert density == pytest.approx(3.8e14, rel=0.05)
    assert energy / (1.4 * MSUN * constants.SPEED_OF_LIGHT_M_S ** 2) == pytest.approx(0.113, abs=0.01)    # Lattimer and Prakash Eq. 33 at beta = 0.172


def test_an_object_without_a_mass_or_a_radius_has_none():
    assert body_strength.star_values(types.SimpleNamespace(mass=None, radius=10.0, yerkes_class="V")) == (None, None)
    assert body_strength.remnant_values("neutron-star", 1.4, 0.0) == (None, None)


def test_a_rogue_giant_and_a_rogue_brown_dwarf_keep_a_density_and_a_binding_energy_only():
    for mass_bin in ("jupiter", "brown-dwarf"):
        draw.set_run_seed(4)
        rogue = RoguePlanet(SystemConfig(), mass_bin=mass_bin)
        assert rogue.density_g_cm3 > 0.0 and rogue.binding_energy_j > 0.0
        assert rogue.roundness_ratio is None and rogue.yield_strength_pa is None


def test_a_terrestrial_rogue_planet_is_round_and_has_a_strength():
    draw.set_run_seed(4)
    rogue = RoguePlanet(SystemConfig(), mass_bin="terrestrial")
    assert rogue.roundness_ratio > 1.0 and rogue.yield_strength_pa > 0.0 and rogue.binding_energy_j > 0.0


def test_an_interstellar_comet_has_what_a_comet_has():
    draw.set_run_seed(6)
    comet = InterstellarComet(SystemConfig())
    assert 0.0 <= comet.ice_mass_fraction <= 1.0 and comet.density_g_cm3 > 0.0
    assert comet.roundness_ratio is not None and comet.yield_strength_pa > 0.0 and comet.binding_energy_j > 0.0
    again = InterstellarComet.from_dict(comet.to_dict(), SystemConfig())
    assert again.binding_energy_j == comet.binding_energy_j and again.roundness_ratio == comet.roundness_ratio


def test_the_rows_of_every_kind_are_stored_with_the_values(mysql_config):
    config = SystemConfig()
    draw.set_run_seed(9)
    with store.get_connection(mysql_config) as conn:
        ids = {
            "black_holes": store.insert_black_hole(conn, BlackHole(config)),
            "neutron_stars": store.insert_neutron_star(conn, NeutronStar(config)),
            "rogue_planets": store.insert_rogue_planet(conn, RoguePlanet(config, mass_bin="terrestrial")),
            "interstellar_comets": store.insert_interstellar_comet(conn, InterstellarComet(config)),
        }
        conn.commit()
        rows = {table: conn.execute(f"SELECT * FROM {table} WHERE id = ?", (row_id,)).fetchone() for table, row_id in ids.items()}
    assert rows["black_holes"]["binding_energy_j"] == -1.0 and rows["black_holes"]["density_g_cm3"] > 0.0
    assert rows["neutron_stars"]["binding_energy_j"] > 0.0 and rows["neutron_stars"]["density_g_cm3"] > 1.0e14
    assert rows["rogue_planets"]["roundness_ratio"] > 1.0 and rows["rogue_planets"]["binding_energy_j"] > 0.0
    comet = rows["interstellar_comets"]
    assert all(comet[name] is not None for name in ("density_g_cm3", "ice_mass_fraction", "roundness_ratio", "yield_strength_pa", "binding_energy_j"))


def test_every_star_of_a_saved_system_keeps_a_density_and_a_binding_energy(mysql_config):
    from planetgen.generation.system import StarSystem
    config = SystemConfig()
    config.PLANETS = True
    seen = 0
    for seed in range(3000, 3006):
        draw.set_run_seed(seed)
        saved = store.save_system(StarSystem(system_config=config, galactic_center_dist_ly=26000.0), config, config=mysql_config)
        with store.get_connection(mysql_config) as conn:
            for row in conn.execute("SELECT yerkes_class, density_g_cm3, binding_energy_j FROM stars WHERE star_system_id = ?", (saved,)).fetchall():
                seen += 1
                assert row["density_g_cm3"] > 0.0 and row["binding_energy_j"] != 0.0
    assert seen >= 6
