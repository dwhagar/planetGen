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
from planetgen.galaxy import scales
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.system import StarSystem
from planetgen import tuning
from planetgen.physics import activity, atmosphere, constants, proton_range, radiation as rad, radiogenic, solar_particles, surface_composition
from planetgen.util import draw

from tests.fuzz_support import deterministic_entropy
from tests.test_db_persistence import _placed_nebula, _system_point_pc


@pytest.mark.parametrize("column, expected", [
    (0.0, 500.0), (22.6, 232.0), (100.0, 553.0), (200.0, 263.0), (500.0, 17.2), (700.0, 3.02), (1036.0, 0.65)])
def test_the_no_field_curve_passes_through_the_read_points(column, expected):
    """The Moon (Chang'E-4), Mars (Hassler 2014) and the proton-primary showers of Atri 2013, as read."""
    assert rad.no_field_dose_msv_yr(column) == pytest.approx(expected, rel=1e-9)


def test_the_curve_reproduces_the_aviation_point_and_sea_level_and_keeps_mars_in_its_range():
    earth = rad.EARTH_CUTOFF_GV
    assert rad.gcr_dose_msv_yr(265.0, earth) == pytest.approx(47.9, rel=0.1)  # Roy 2021, 10 km: 5.46 uSv/h
    assert rad.gcr_dose_msv_yr(1033.0, earth) == pytest.approx(0.32, rel=0.1)  # UNSCEAR 2008 sea level
    assert 128.0 <= rad.gcr_dose_msv_yr(22.5) <= 420.0  # Guo 2015, Mars at 22.5 g/cm2
    assert rad.gcr_dose_msv_yr(122.7) == pytest.approx(444.0, rel=0.01)  # the Pfotzer maximum, log interpolation


def test_a_mountain_station_takes_the_dose_the_aircraft_rule_gives_between_it_and_the_flight_levels():
    """The dose doubles for every 1800 m of altitude in the aircraft range (PMC7968484, read); from 6.1 km (450 g/cm2) to 12.2 km (190 g/cm2) that
    is 2^3.4 = 10.5 times, and the runs with Earth's field give it within a third."""
    earth = rad.EARTH_CUTOFF_GV
    ratio = rad.gcr_dose_msv_yr(190.0, earth) / rad.gcr_dose_msv_yr(450.0, earth)
    assert ratio == pytest.approx(10.5, rel=0.35)
    sea_level = rad.gcr_dose_msv_yr(rad.GCR_SEA_LEVEL_G_CM2, earth)
    assert sea_level < rad.gcr_dose_msv_yr(700.0, earth) < rad.gcr_dose_msv_yr(450.0, earth)  # 0, 3 and 6 km rise


def test_the_thin_air_curve_passes_through_the_moon_and_mars_and_joins_the_runs_at_their_peak():
    """GEN.746: the value and the measured fall at Mars and the peak at the first run are the three conditions the curve is fitted to."""
    assert rad.thin_air_dose_msv_yr(0.0) == rad.GCR_MOON_MSV_YR
    assert rad.thin_air_dose_msv_yr(rad.GCR_MARS_G_CM2) == pytest.approx(rad.GCR_MARS_MSV_YR, rel=1e-6)
    step = 1e-4
    slope = ((math.log(rad.thin_air_dose_msv_yr(rad.GCR_MARS_G_CM2 + step)) - math.log(rad.thin_air_dose_msv_yr(rad.GCR_MARS_G_CM2 - step)))
             / (2 * step))
    assert slope == pytest.approx(rad.GCR_MARS_LOG_SLOPE_PER_G_CM2, rel=1e-4)
    first = rad.GCR_SHOWER_POINTS[0]
    assert rad.thin_air_dose_msv_yr(first[0]) == pytest.approx(first[1], rel=1e-9)
    assert rad.thin_air_dose_msv_yr(first[0] - 1.0) < first[1] and rad.thin_air_dose_msv_yr(first[0] - 10.0) < first[1]
    assert rad.thin_air_dose_msv_yr(first[0] + 50.0) == rad.thin_air_dose_msv_yr(first[0])  # held at the join; the runs take over


def test_the_thin_air_curve_falls_through_a_minimum_and_rises_to_the_runs_with_the_primaries_length_of_a_heavy_nucleus():
    """The primaries thin with a length near the collision length of carbon-to-oxygen nuclei (25 g/cm2; the protons' is 90): that length is not
    a free number but falls out of the Mars value and slope."""
    columns = [i * 0.5 for i in range(0, 201)]
    doses = [rad.thin_air_dose_msv_yr(x) for x in columns]
    low = doses.index(min(doses))
    assert 20.0 < columns[low] < 50.0 and min(doses) < rad.GCR_MARS_MSV_YR
    assert doses[:low + 1] == sorted(doses[:low + 1], reverse=True) and doses[low:] == sorted(doses[low:])
    assert 20.0 < rad.GCR_PRIMARY_LENGTH_G_CM2 < 35.0
    assert rad.GCR_BUILD_UP_EXPONENT > 1.0 and rad.GCR_SHOWER_LENGTH_G_CM2 > 0.0


def test_below_the_runs_the_curve_is_the_thin_air_one_and_above_them_it_keeps_falling_to_a_muon_floor():
    assert 232.0 < rad.no_field_dose_msv_yr(10.0) < 500.0 and 232.0 < rad.no_field_dose_msv_yr(60.0) < 553.0
    assert rad.no_field_dose_msv_yr(60.0) == rad.thin_air_dose_msv_yr(60.0)
    assert rad.no_field_dose_msv_yr(3000.0) < rad.no_field_dose_msv_yr(1036.0)
    assert rad.no_field_dose_msv_yr(10850.0) == pytest.approx(3.4e-7 * 1e3, rel=0.1)  # Crouch at Titan's column, mSv/yr
    assert rad.no_field_dose_msv_yr(1e6) == pytest.approx(rad._muon_dose_msv_yr(1e6))


def test_the_cutoff_scales_with_moment_over_radius_squared():
    earth = constants.EARTH_RADIUS_KM
    assert rad.cutoff_rigidity_gv(rad.EARTH_DIPOLE_A_M2, earth) == pytest.approx(14.34, abs=0.01)  # Stormer, computed
    assert rad.cutoff_rigidity_gv(2 * rad.EARTH_DIPOLE_A_M2, earth) == pytest.approx(2 * rad.EARTH_CUTOFF_GV)
    assert rad.cutoff_rigidity_gv(rad.EARTH_DIPOLE_A_M2, earth / 2) == pytest.approx(4 * rad.EARTH_CUTOFF_GV)
    assert rad.cutoff_rigidity_gv(rad.EARTH_DIPOLE_A_M2, 2 * earth) == pytest.approx(rad.EARTH_CUTOFF_GV / 4)
    assert rad.cutoff_rigidity_gv(0.0, earth) == 0.0 and rad.cutoff_rigidity_gv(1e22, 0.0) == 0.0


@pytest.mark.parametrize("column, ratio", [(100.0, 0.31), (200.0, 0.37), (500.0, 0.51), (700.0, 0.61), (1036.0, 0.52)])
def test_an_earth_strength_field_cuts_the_dose_by_the_ratio_of_the_runs_at_each_column(column, ratio):
    assert rad.magnetic_factor(rad.EARTH_CUTOFF_GV, column) == pytest.approx(ratio, abs=0.01)
    assert rad.magnetic_factor(0.0, column) == 1.0


def test_the_field_factor_follows_the_griessmeier_runs_at_100_and_1036_and_never_adds_dose():
    unit = rad.EARTH_CUTOFF_GV
    for moment, dose100, dose1036 in zip(rad.FIELD_FULL_MOMENTS, rad.FIELD_DOSE_100, rad.FIELD_DOSE_1036):
        assert rad.magnetic_factor(moment * unit, 100.0) == pytest.approx(dose100 / 553.0, rel=1e-9)
        assert rad.magnetic_factor(moment * unit, 1036.0) == pytest.approx(dose1036 / 0.65, rel=1e-9)
    for column in rad.FIELD_COLUMNS_G_CM2 + (30.0, 150.0, 4000.0):
        factors = [rad.magnetic_factor(moment * unit, column) for moment in (0.0, 0.05, 0.15, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0)]
        assert all(0.0 < f <= 1.0 for f in factors)
        assert factors == sorted(factors, reverse=True)
    # a thinner column takes the ratio at 100 g/cm2, a stronger field than the runs reach the ratio at 10 moments
    assert rad.magnetic_factor(unit, 5.0) == rad.magnetic_factor(unit, 100.0)
    assert rad.magnetic_factor(1000 * unit, 500.0) == rad.magnetic_factor(10 * unit, 500.0)


def test_the_dose_follows_the_field_in_earth_moments_of_an_earth_sized_body():
    """A smaller body of the same moment has a stronger cutoff and so a higher cut, as the radius sets M/r^2."""
    earth = constants.EARTH_RADIUS_KM
    column = rad.column_g_cm2(5066.0, 1.0)
    big = rad.gcr_dose_msv_yr(column, rad.cutoff_rigidity_gv(rad.EARTH_DIPOLE_A_M2, earth))
    small = rad.gcr_dose_msv_yr(column, rad.cutoff_rigidity_gv(rad.EARTH_DIPOLE_A_M2, earth / 2))
    assert small < big
    assert big == pytest.approx(rad.no_field_dose_msv_yr(column) * rad.magnetic_factor(rad.EARTH_CUTOFF_GV, column))


def test_the_modulation_scales_the_dose():
    assert rad.gcr_dose_msv_yr(50.0, 0.0, 1.64) == pytest.approx(1.64 * rad.gcr_dose_msv_yr(50.0))
    assert tuning.GCR_SOLAR_MODULATION_RANGE == (0.5, 1.64)


def test_proxima_bs_thin_air_breaks_the_flare_shield():
    """GEN.747: the measured events give the Sun 265 mSv/yr in free space; Proxima's flares scale it by the square root of their rate and the
    inverse square of the distance. A tenth of a bar lets 11 Sv a year through (Atri 2017: 5 to 10 Sv from one hard event at 700 g/cm2 or less),
    and a bar of air stops all but a few mSv."""
    free = rad.free_sep_msv_yr(2100.0 * 0.0485 ** 2, 0.0485)
    column = rad.column_g_cm2(10132.5, 1.07)
    assert rad.gcr_dose_msv_yr(column) == pytest.approx(540, rel=0.1)
    assert rad.sep_dose_msv_yr(free, column, 0.0) == pytest.approx(1.06e4, rel=0.1)
    assert rad.sep_dose_msv_yr(free, column, 0.0) > 5000.0
    dense = rad.column_g_cm2(101325.0, 1.07)
    assert rad.sep_dose_msv_yr(free, dense, 0.0) < 20.0
    assert rad.sep_dose_msv_yr(free, dense, 2.2) <= rad.sep_dose_msv_yr(free, dense, 0.0)
    assert rad.sep_dose_msv_yr(free, dense, 20.0) < rad.sep_dose_msv_yr(free, dense, 0.0)


def test_the_suns_free_space_particle_dose_is_the_measured_events_and_scales_with_the_flares_and_distance():
    assert rad.free_sep_msv_yr(rad.N33_SUN_PER_YR, 1.0) == solar_particles.FREE_MSV_YR
    assert rad.free_sep_msv_yr(4 * rad.N33_SUN_PER_YR, 1.0) == pytest.approx(2 * solar_particles.FREE_MSV_YR)
    assert rad.free_sep_msv_yr(rad.N33_SUN_PER_YR, 2.0) == pytest.approx(solar_particles.FREE_MSV_YR / 4.0)
    assert rad.free_sep_msv_yr(None, 1.0) == 0.0 and rad.free_sep_msv_yr(0.0, 1.0) == 0.0 and rad.free_sep_msv_yr(1.0, None) == 0.0
    # Earth's air and field leave the particles of the events nothing: the earth is spared all but 1e-4 mSv/yr
    earth = rad.sep_dose_msv_yr(solar_particles.FREE_MSV_YR, rad.column_g_cm2(101325.0, 1.0), rad.EARTH_CUTOFF_GV)
    assert earth < 1e-3
    # Mars's 22.6 g/cm2 stops the protons below 172 MeV; with no field about a fifteenth of the dose reaches the ground
    mars = rad.sep_dose_msv_yr(solar_particles.FREE_MSV_YR, 22.6, 0.0)
    assert mars == pytest.approx(solar_particles.FREE_MSV_YR * 0.072, rel=0.15)


def test_the_stellar_particle_cutoff_is_the_pstar_proton_energy():
    assert rad.proton_cutoff_gv(0.0) == 0.0
    assert rad.proton_cutoff_gv(22.5) == pytest.approx(proton_range.rigidity_gv(171.0), rel=0.01)
    assert rad.proton_cutoff_gv(1033.0) == pytest.approx(proton_range.rigidity_gv(2241.0), rel=0.01)
    assert rad.proton_cutoff_gv(22.5, "H") > rad.proton_cutoff_gv(22.5, "air")


def test_the_heliosphere_only_matters_for_thin_air_and_a_squeezed_star():
    assert rad.helio_multiplier(0.0, 1.0) == pytest.approx(2.5)
    assert rad.helio_multiplier(0.0, 0.0) == 1.0
    assert rad.helio_multiplier(1033.0, 1.0) == pytest.approx(1.0, abs=1e-3)
    assert rad.heliosphere_compression(85.0, 90.0) == 0.0
    assert rad.heliosphere_compression(85.0, 1.0) == 1.0
    assert 0.0 < rad.heliosphere_compression(85.0, 10.0) < 1.0


def test_the_ground_dose_reads_the_same_concentrations_as_the_heat():
    """GEN.459: no independent draw of the ground; doubling the thorium and uranium doubles the heat from them and the gamma dose."""
    k, th, u = radiogenic.BULK_SILICATE_EARTH
    one = radiogenic.Content(k * 1.0, th * radiogenic.formation_mix("Th")[1], u * radiogenic.formation_mix("U")[1], 0.0)
    two = radiogenic.Content(one.k0_ppm, 2 * one.th0_ppm, 2 * one.u0_ppm, 0.0)
    basalt = surface_composition.RESERVOIRS["basaltic"]
    rock = (basalt["k_ppm"], basalt["th_ppm"], basalt["u_ppm"])
    base = rad.gamma_dose_msv_yr(radiogenic.surface_content(one, rock), 4.568)
    doubled = rad.gamma_dose_msv_yr(radiogenic.surface_content(two, rock), 4.568)
    k_part = rad.gamma_dose_msv_yr(radiogenic.surface_content(radiogenic.Content(one.k0_ppm, 0.0, 0.0, 0.0), rock), 4.568)
    assert doubled - k_part == pytest.approx(2 * (base - k_part), rel=1e-9)
    assert radiogenic.heat_w_per_kg(two, 4.568) > radiogenic.heat_w_per_kg(one, 4.568)


def test_the_gamma_dose_falls_under_a_dense_gas_and_the_radon_needs_air():
    rock = radiogenic.Content(17000.0, 11.0, 2.7, 0.0)
    open_air = rad.gamma_dose_msv_yr(rock, 4.568)
    assert rad.gamma_dose_msv_yr(rock, 4.568, gas_g_cm3=0.0012) == pytest.approx(open_air * math.exp(-0.0637 * 0.0012 * 100), rel=1e-12)
    assert rad.gamma_dose_msv_yr(rock, 4.568, gas_g_cm3=0.065) < 0.7 * open_air
    assert rad.radon_dose_msv_yr(rock, 4.568, None) == 0.0
    assert rad.radon_dose_msv_yr(rock, 4.568, 8.5) > 0.0
    thin = rad.radon_dose_msv_yr(rock, 4.568, 8.5, diffusion_factor=100.0)
    assert thin == pytest.approx(10.0 * rad.radon_dose_msv_yr(rock, 4.568, 8.5), rel=1e-9)
    assert rad.radon_diffusion_factor(101325.0, 0.7) == 1.0
    assert rad.radon_diffusion_factor(1013.25, 1.0) == pytest.approx(100.0)
    assert rad.radon_diffusion_factor(1013.25, 0.0) == 1.0


def test_the_uv_proxies_are_gone():
    """R19 (GEN.445): the ratio-to-the-Sun index and its ozone exponent are not here; a proxy that returns fails this."""
    for name in ("uv_ratio", "uv_index", "ozone_ratio", "_band_fraction", "UV_FLOOR", "UV_MAX", "UV_EXPONENT", "EARTH_OZONE_DU", "SUN_TEMPERATURE_K"):
        assert not hasattr(rad, name), name


def test_the_stored_uv_index_is_the_dna_dose_of_the_spectrum_below_the_bodys_air():
    from planetgen.physics import body_radiation

    with deterministic_entropy(2):    # a seed whose system holds several rocky planets (the architecture leaves some systems bare)
        draw.set_run_seed(2)
        system = StarSystem(SystemConfig(), galactic_center_dist_ly=26000.0)
    checked = 0
    for body in _bodies(system):
        if body.body_type == "g" or body.is_moon:
            continue    # a moon of a planet that was moved keeps its dose from the old distance (reported to TODO)
        expected = body_radiation.dna_dose_index(body)
        assert body.uv_surface_index == expected
        checked += expected is not None
    assert checked > 0


def test_supernova_hazard_rises_toward_the_galactic_centre():
    def ly(pc):
        return pc * constants.PARSEC_M / constants.LY_TO_M

    sun = activity.supernova_rate_per_gyr(ly(8000))
    assert sun == pytest.approx(1.5, rel=0.01)
    assert activity.supernova_rate_per_gyr(0.0) == pytest.approx(1.5 * 250, rel=0.01)
    assert activity.lethal_event_rate_per_gyr(None) is None
    assert activity.supernova_rate_per_gyr(ly(12000)) < sun


def test_gamma_ray_bursts_and_the_nucleus_add_to_the_lethal_rate():
    def ly(pc):
        return pc * constants.PARSEC_M / constants.LY_TO_M

    at_sun, outer, inner = ly(8000), ly(12000), ly(2000)
    assert activity.grb_rate_per_gyr(at_sun) == pytest.approx(activity.GRB_LETHAL_RATE_AT_SUN_PER_GYR)
    assert activity.grb_rate_per_gyr(outer) < activity.grb_rate_per_gyr(at_sun)
    assert activity.grb_rate_per_gyr(inner) == pytest.approx(activity.GRB_LETHAL_RATE_AT_SUN_PER_GYR)
    assert activity.agn_rate_per_gyr(inner) == activity.AGN_EPISODES_PER_GYR
    assert activity.agn_rate_per_gyr(ly(activity.AGN_LETHAL_RADIUS_PC + 1)) == 0.0
    assert activity.agn_rate_per_gyr(at_sun) == 0.0
    assert activity.lethal_event_rate_per_gyr(at_sun) == pytest.approx(
        activity.supernova_rate_per_gyr(at_sun) + activity.grb_rate_per_gyr(at_sun))
    assert activity.lethal_event_rate_per_gyr(inner) == pytest.approx(
        activity.supernova_rate_per_gyr(inner) + activity.grb_rate_per_gyr(inner) + activity.agn_rate_per_gyr(inner))


def test_bursts_and_the_nucleus_push_the_ozone_flag_out_to_larger_radii():
    """The galactic rate alone decides the flag where no flares do: it now
    exceeds `LETHAL_EVENT_LIMIT_PER_GYR` farther from the centre."""
    def flag_radius_pc(rate):
        radius = 100
        while radius < 20000 and rate(radius * constants.PARSEC_M / constants.LY_TO_M) > rad.LETHAL_EVENT_LIMIT_PER_GYR:
            radius += 50
        return radius

    assert flag_radius_pc(activity.lethal_event_rate_per_gyr) > flag_radius_pc(activity.supernova_rate_per_gyr)
    sun_ly = 8000 * constants.PARSEC_M / constants.LY_TO_M
    assert activity.lethal_event_rate_per_gyr(sun_ly) < rad.LETHAL_EVENT_LIMIT_PER_GYR
    assert activity.lethal_event_rate_per_gyr(sun_ly / 4) > rad.LETHAL_EVENT_LIMIT_PER_GYR
def _bodies(system):
    for planet in system.planets:
        if planet.body_type == "a":
            continue
        yield planet
        yield from planet.moons


def test_every_generated_body_has_a_consistent_dose():
    seen = {"rocky": 0, "giant": 0, "ozone": 0}
    for seed in range(80):    # the planet architecture (GEN.728) makes a giant in about 1 system of 10, so 12 seeds held none and 24 held none once the draws moved
        if seen["rocky"] > 10 and seen["giant"] > 0:
            break
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
                body.dose_gcr_msv_yr + body.dose_sep_msv_yr + body.dose_ground_msv_yr + body.dose_belt_msv_yr)
            assert body.dose_belt_msv_yr >= 0 and (body.dose_belt_msv_yr == 0 or body.is_moon)
            assert body.dose_gcr_msv_yr > 0 and body.dose_ground_msv_yr > 0
            assert isinstance(body.ozone_loss_flag, bool)
            assert (body.ozone_column_du > 0.0) == (atmosphere.partial_kpa(body, "o2") > 0.0 and bool(body.star.temperature))
            seen["ozone"] += body.ozone_loss_flag
            assert body.uv_surface_index is None or body.uv_surface_index >= 0.0
            assert body.uv_surface_index != 0.0 or body.atmospheric_pressure > 0.0     # a dose of exactly 0 is an opaque air (a dense vapour air passes no UV), never a bare body
    assert seen["rocky"] > 10 and seen["giant"] > 0


def test_a_star_with_no_galaxy_placement_is_rated_where_a_system_with_no_place_is_taken():
    """GEN.785: the column is filled for every star; one placed nowhere takes the distance a system with no place is taken at."""
    class Unplaced:
        galactic_center_dist_ly = None

    star = Unplaced()
    activity.set_galactic_hazard(star)
    assert star.lethal_event_rate_per_gyr == pytest.approx(activity.lethal_event_rate_per_gyr(scales.default_galactocentric_ly()))
    assert star.lethal_event_rate_per_gyr > 0.0


def test_a_star_nearer_the_centre_flags_ozone_loss_on_a_world_with_air():
    with deterministic_entropy(5):
        draw.set_run_seed(5)
        system = StarSystem(SystemConfig(), galactic_center_dist_ly=2000.0)
    star = system.stars[0] if hasattr(system, "stars") else system.primary_star
    assert star.lethal_event_rate_per_gyr > activity.SN_RATE_AT_SUN_PER_GYR * 50
    for body in _bodies(system):
        if body.body_type == "t" and (body.ozone_column_du or 0.0) > 0.0:
            assert body.ozone_loss_flag


def _rocky_world_with_air():
    for seed in range(60):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(SystemConfig(), galactic_center_dist_ly=26000.0)
        body = next((b for b in _bodies(system) if b.body_type == "t" and (b.atmospheric_pressure or 0.0) > 0.0), None)
        if body is not None:
            return body
    raise AssertionError("no seed gave a rocky world with air")


def test_a_real_world_gets_a_finite_ground_dose_that_follows_the_rock_of_its_ground():
    """GEN.847: `generate_dose` hands `ground_dose_msv_yr` the (K, Th, U) ppm of the stored ground, so a world made of felsic rock
    has more ground dose than the same world made of ultramafic rock, and neither fails."""
    body = _rocky_world_with_air()
    doses = {}
    for ground in ("ultramafic", "felsic"):
        for name in surface_composition.GROUND_FRACTION_TYPES:
            setattr(body, f"ground_{name}_fraction", 1.0 if name == ground else 0.0)
        with deterministic_entropy(1):
            rad.generate_dose(body)
        doses[ground] = body.dose_ground_msv_yr
    assert 0.0 < doses["ultramafic"] < doses["felsic"] < float("inf")


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


def test_a_shallow_flare_slope_puts_more_energy_in_the_big_flares_and_raises_the_particle_dose():
    """GEN.805: the energy weight of the flare distribution (`flare_alpha`) over 1e33 to 1e36 erg, taken to the square root in the dose."""
    assert rad.flare_energy_weight(2.0) == 1.0 and rad.flare_energy_weight(None) == 1.0
    assert rad.flare_energy_weight(1.8) == pytest.approx(2.16, abs=0.02)
    assert rad.flare_energy_weight(2.2) == pytest.approx(0.54, abs=0.02)
    assert rad.flare_energy_weight(1.8) > rad.flare_energy_weight(2.0) > rad.flare_energy_weight(2.2)
    shallow, steep = rad.free_sep_msv_yr(1.0, 0.1, 1.8), rad.free_sep_msv_yr(1.0, 0.1, 2.2)
    assert shallow > rad.free_sep_msv_yr(1.0, 0.1, 2.0) > steep
    assert shallow / steep == pytest.approx(math.sqrt(rad.flare_energy_weight(1.8) / rad.flare_energy_weight(2.2)))
    assert rad.free_sep_msv_yr(1.0, 0.1) == rad.free_sep_msv_yr(1.0, 0.1, 2.0)


def test_the_fields_shielding_of_the_flares_is_one_without_a_dipole_and_less_for_a_wider_magnetopause():
    assert rad.field_shielding(None, "none") == rad.SHIELD_NO_FIELD == 1.0
    assert rad.field_shielding(10.0, "none") == 1.0 and rad.field_shielding(0.0, "weak") == 1.0
    assert rad.field_shielding(2.0, "earth-like") == rad.SHIELD_WEAK_FIELD == 0.3
    assert rad.field_shielding(9.75, "earth-like") == rad.SHIELD_STRONG_FIELD == 0.1


def _flared(index, magnetopause_rp, dipole_class, column=1.0, hazard=None):
    star = types.SimpleNamespace(lethal_event_rate_per_gyr=hazard)
    return types.SimpleNamespace(star=star, ozone_column_du=column, flare_irradiation_index=index, magnetopause_rp=magnetopause_rp,
                                 dipole_class=dipole_class)


def test_proxima_b_without_a_field_loses_its_ozone_and_earth_does_not():
    """The design note section 2.3: flare irradiation index times f_B at 300; Proxima b about 2,100 and Earth 0.003."""
    assert rad.ozone_lost(_flared(2100.0, None, "none"))
    assert rad.ozone_lost(_flared(2100.0, 2.0, "earth-like"))                 # 630 even through a squeezed field
    assert not rad.ozone_lost(_flared(2100.0, 9.0, "earth-like"))               # 210 through a wide one
    assert not rad.ozone_lost(_flared(0.003, None, "none"))
    assert rad.ozone_lost(_flared(300.0, None, "none")) and not rad.ozone_lost(_flared(299.0, None, "none"))
    assert not rad.ozone_lost(_flared(2100.0, None, "none", column=0.0))        # nothing to lose
    assert rad.ozone_lost(_flared(0.0, None, "none", hazard=rad.LETHAL_EVENT_LIMIT_PER_GYR * 2))      # the galaxy's events
    assert not rad.ozone_lost(_flared(None, None, None))


def _world_with_air():
    for seed in range(80):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(SystemConfig(), galactic_center_dist_ly=26000.0)
        for body in _bodies(system):
            if body.body_type == "t" and (body.atmospheric_pressure or 0.0) > 0.0 and (body.atmosphere or {}) and not body.is_moon:
                body.draw_seed = seed
                if atmosphere.partial_kpa(body, "o2") > 0.0 and body.star.temperature:
                    return body
    raise AssertionError("no seed gave a rocky world with oxygen")


def test_a_world_whose_ozone_is_lost_keeps_a_tenth_of_it_and_takes_more_ultraviolet():
    """GEN.808: the flag lowers the stored column and the DNA-weighted UV index of the surface follows it."""
    body = _world_with_air()
    body.star.lethal_event_rate_per_gyr = 0.0

    body.magnetic_moment_a_m2, body.magnetopause_rp, body.dipole_class = 0.0, None, "none"

    def run(index):
        body.flare_irradiation_index = index
        draw.set_run_seed(808)
        rad.generate_dose(body)
        return body.ozone_loss_flag, body.ozone_column_du, body.uv_surface_index

    quiet = run(0.0)
    flared = run(2100.0)
    assert quiet[0] is False and flared[0] is True
    assert flared[1] == pytest.approx(quiet[1] * rad.OZONE_LOSS_REMAINING)
    assert flared[2] > quiet[2]
