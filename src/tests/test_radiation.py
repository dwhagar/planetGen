"""
Surface radiation dose, UV and galactic hazards (GEN.87,
`planetgen.physics.radiation`): the design doc's section 4 numbers, what
every generated body stores and the containment pass that updates the
heliosphere multiplier.
"""

import math
import types

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.system import StarSystem
from planetgen.physics import activity, constants, radiation as rad
from planetgen.util import draw

from tests.fuzz_support import deterministic_entropy
from tests.test_db_persistence import _placed_nebula, _system_point_pc


@pytest.mark.parametrize("column, expected", [(0.0, 500.0), (16.4, 237.0), (100.0, 110.0), (1033.0, 0.37)])
def test_the_cosmic_ray_curve_passes_through_its_anchors(column, expected):
    assert rad.gcr_dose_msv_yr(column) == pytest.approx(expected, rel=0.02)


def test_an_earth_dipole_cuts_cosmic_rays_to_six_tenths():
    cutoff = rad.cutoff_rigidity_gv(rad.EARTH_DIPOLE_A_M2, constants.EARTH_RADIUS_KM)
    assert cutoff == pytest.approx(14.9)
    assert rad.magnetic_factor(cutoff) == pytest.approx(0.6, abs=0.01)
    assert rad.magnetic_factor(0.0) == 1.0
    assert rad.gcr_dose_msv_yr(rad.column_g_cm2(5066.0, 1.0), cutoff) == pytest.approx(88, rel=0.1)


def test_proxima_bs_thin_air_breaks_the_flare_shield():
    free = rad.free_sep_msv_yr(2100.0 * 0.0485 ** 2, 0.0485)
    column = rad.column_g_cm2(10132.5, 1.07)
    assert rad.gcr_dose_msv_yr(column) == pytest.approx(108, rel=0.1)
    assert rad.sep_dose_msv_yr(free, column, 0.0) == pytest.approx(2300, rel=0.5)
    dense = rad.column_g_cm2(101325.0, 1.07)
    assert rad.sep_dose_msv_yr(free, dense, 0.0) < 10.0
    assert rad.sep_dose_msv_yr(free, dense, 2.2) < rad.sep_dose_msv_yr(free, dense, 0.0)


def test_the_heliosphere_only_matters_for_thin_air_and_a_squeezed_star():
    assert rad.helio_multiplier(0.0, 1.0) == pytest.approx(2.5)
    assert rad.helio_multiplier(0.0, 0.0) == 1.0
    assert rad.helio_multiplier(1033.0, 1.0) == pytest.approx(1.0, abs=1e-3)
    assert rad.heliosphere_compression(85.0, 90.0) == 0.0
    assert rad.heliosphere_compression(85.0, 1.0) == 1.0
    assert 0.0 < rad.heliosphere_compression(85.0, 10.0) < 1.0


def test_the_crust_ages_like_the_radiogenic_heat():
    earth = constants.EARTH_MASS_TO_KG
    young = rad.ground_dose_msv_yr(0.5, earth, 6371.0, 1.0, False)
    now = rad.ground_dose_msv_yr(4.5, earth, 6371.0, 1.0, False)
    assert now == pytest.approx(0.48, rel=0.02)
    assert young == pytest.approx(0.48 * 3.5, rel=0.25)
    assert rad.ground_dose_msv_yr(4.5, earth, 6371.0, 1.0, True) == pytest.approx(0.98, rel=0.02)


def test_uv_follows_the_star_and_the_ozone_layer():
    assert rad.uv_ratio(5772.0) == pytest.approx(1.0)
    assert rad.uv_ratio(7200.0) == pytest.approx(2.8, rel=0.5)
    assert rad.uv_ratio(2800.0) == rad.UV_FLOOR
    sun = types.SimpleNamespace(temperature=5772.0, luminosity=constants.SOLAR_LUMINOSITY)
    assert rad.uv_index(sun, 1.0, 21.2) == pytest.approx(1.0, rel=0.05)
    assert rad.uv_index(sun, 1.0, 0.0) == pytest.approx(rad.UV_MAX)
    assert rad.uv_index(sun, 1.0, 5.0) > rad.uv_index(sun, 1.0, 21.2)
    assert rad.uv_index(types.SimpleNamespace(temperature=None, luminosity=None), 1.0, 21.2) is None


def test_supernova_hazard_rises_toward_the_galactic_centre():
    sun = activity.lethal_event_rate_per_gyr(8000 * constants.PARSEC_M / constants.LY_TO_M)
    assert sun == pytest.approx(1.5, rel=0.01)
    assert activity.lethal_event_rate_per_gyr(0.0) == pytest.approx(1.5 * 250, rel=0.01)
    assert activity.lethal_event_rate_per_gyr(None) is None
    outer = activity.lethal_event_rate_per_gyr(12000 * constants.PARSEC_M / constants.LY_TO_M)
    assert outer < sun


def _bodies(system):
    for planet in system.planets:
        if planet.body_type == "a":
            continue
        yield planet
        yield from planet.moons


def test_every_generated_body_has_a_consistent_dose():
    seen = {"rocky": 0, "giant": 0, "ozone": 0}
    for seed in range(12):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(SystemConfig(), galactic_center_dist_ly=26000.0)
        for body in _bodies(system):
            if body.body_type == "g":
                seen["giant"] += 1
                assert all(getattr(body, field) is None for field in rad.BODY_FIELDS)
                continue
            seen["rocky"] += 1
            assert body.dose_helio_mult == 1.0
            assert body.surface_dose_msv_yr == pytest.approx(
                body.dose_gcr_msv_yr + body.dose_sep_msv_yr + body.dose_ground_msv_yr)
            assert body.dose_gcr_msv_yr > 0 and body.dose_ground_msv_yr > 0
            assert isinstance(body.ozone_loss_flag, bool)
            seen["ozone"] += body.ozone_loss_flag
            assert body.uv_surface_index is None or body.uv_surface_index > 0
    assert seen["rocky"] > 10 and seen["giant"] > 0


def test_a_star_nearer_the_centre_flags_ozone_loss_on_a_world_with_air():
    with deterministic_entropy(5):
        draw.set_run_seed(5)
        system = StarSystem(SystemConfig(), galactic_center_dist_ly=2000.0)
    star = system.stars[0] if hasattr(system, "stars") else system.primary_star
    assert star.lethal_event_rate_per_gyr > activity.SN_RATE_AT_SUN_PER_GYR * 50
    for body in _bodies(system):
        if body.body_type == "t" and (body.p_o2_kpa or 0.0) >= rad.O3_MIN_PO2_KPA:
            assert body.ozone_loss_flag


def _system_with_rocky_worlds():
    for seed in range(40):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.BINARY_SYSTEM = False
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(system_config=cfg, galactic_center_dist_ly=26000.0)
        if any(body.body_type == "t" for body in _bodies(system)):
            return cfg, system
    raise AssertionError("no seed gave a rocky world")


def test_entering_and_leaving_a_cloud_moves_the_cosmic_ray_multiplier(mysql_config):
    cfg, system = _system_with_rocky_worlds()
    sector = SpaceSector("Dosed", edge_ly=11.5)
    sector.add_system(system, position=(0.0, 0.0, 0.0), system_config=cfg)
    sector_id = store.save_sector(sector, config=mysql_config, galaxy_position={
        "center_x_pc": 100.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 100.0,
    })
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            # An airless world is the case the multiplier is built for.
            conn.execute("UPDATE planets SET atmospheric_pressure_pa = 0, dose_gcr_msv_yr = 500,"
                         " dose_sep_msv_yr = 1, dose_ground_msv_yr = 0.5, surface_dose_msv_yr = 501.5,"
                         " dose_helio_mult = 1 WHERE body_type = 't'")
            cloud = _placed_nebula(conn, sector_id, _system_point_pc(conn, sector_id), radius_ly=3.0, nebula_class="C")
            conn.execute("UPDATE nebulae SET density_cm3 = 3000, temperature_k = 20 WHERE id = ?", (cloud,))
            home = conn.execute("SELECT center_x_pc FROM nebulae WHERE id = ?", (cloud,)).fetchone()["center_x_pc"]
            conn.execute("UPDATE nebulae SET center_x_pc = center_x_pc + 50 WHERE id = ?", (cloud,))
            store.refresh_containment(conn, [sector_id])
            conn.execute("UPDATE nebulae SET center_x_pc = ? WHERE id = ?", (home, cloud))
            store.refresh_containment(conn, [sector_id])
        rows = conn.execute("SELECT dose_helio_mult, surface_dose_msv_yr FROM planets WHERE body_type = 't'").fetchall()
        assert rows
        for row in rows:
            assert row["dose_helio_mult"] == pytest.approx(2.5, abs=0.01)
            assert row["surface_dose_msv_yr"] == pytest.approx(500 * row["dose_helio_mult"] + 1.5)
        with conn:
            conn.execute("UPDATE nebulae SET center_x_pc = center_x_pc + 50 WHERE id = ?", (cloud,))
            store.refresh_containment(conn, [sector_id])
        rows = conn.execute("SELECT dose_helio_mult, surface_dose_msv_yr FROM planets WHERE body_type = 't'").fetchall()
        assert all(row["dose_helio_mult"] == 1.0 and row["surface_dose_msv_yr"] == pytest.approx(501.5)
                   for row in rows)
    finally:
        conn.close()
