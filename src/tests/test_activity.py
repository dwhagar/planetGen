"""
Stellar activity and planetary magnetic fields (GEN.86,
`planetgen.physics.activity` and `planetgen.physics.magnetism`).
"""

import math
import types

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
from planetgen.generation.system import StarSystem
from planetgen.physics import activity, constants, magnetism
from planetgen.util import draw

from tests.fuzz_support import deterministic_entropy

SUN = constants.SOLAR_LUMINOSITY


def test_the_sun_today_matches_the_design_validation():
    t_sat, decay = activity._interp(activity.ACTIVITY_TABLE, 1.0)
    l_x = activity.lx_lbol_at(4.57, activity.saturated_lx_lbol(1.0), t_sat, decay) * SUN
    assert 1e27 <= l_x * activity.ERG_PER_J <= 2e27
    flux = activity.coronal_xuv_w(l_x) / (4 * math.pi * constants.AU_TO_M ** 2) / constants.EARTH_XUV_FLUX_W_M2
    assert flux == pytest.approx(1.2, abs=0.2)
    # Saturated XUV of a G star is 3.3e-3 of its luminosity (section 2.1).
    assert activity.coronal_xuv_w(7.4e-4 * SUN) / SUN == pytest.approx(3.3e-3, rel=0.03)


def test_saturation_falls_through_the_f_stars():
    assert activity.saturated_lx_lbol(0.5) == 7.4e-4
    assert activity.saturated_lx_lbol(1.4) == pytest.approx(5e-5)
    assert activity.saturated_lx_lbol(5.0) == pytest.approx(1e-7)
    assert activity._interp(activity.ACTIVITY_TABLE, 0.1) == (4.0, 1.0)
    assert activity._interp(activity.ACTIVITY_TABLE, 2.0) == (0.1, 2.0)


def test_a_hot_photosphere_emits_its_own_xuv():
    assert activity.blackbody_fraction_above(13.6, 5772) < 1e-8
    assert activity.blackbody_fraction_above(13.6, 40000) == pytest.approx(0.41, abs=0.02)
    assert activity.blackbody_fraction_above(0.0001, 5772) == pytest.approx(1.0, abs=1e-3)


def test_the_flare_table_contains_proxima_and_the_sun():
    assert activity.flare_log_n33(0.12, 4.85) == 0.8
    assert activity.flare_log_n33(1.0, 4.6) == -2.5
    assert activity.flare_log_n33(2.0, 0.5) == activity.HOT_STAR_LOG_N33


def _star(mass_sol, age_gy, yerkes="V", temperature=4000.0):
    return types.SimpleNamespace(mass=mass_sol * constants.SOLAR_MASS_TO_KG, luminosity=0.01 * SUN, age=age_gy,
                                 temperature=temperature, yerkes_class=yerkes)


def test_m_dwarfs_stay_saturated_longer():
    draw.set_run_seed(86)
    young_m = [_star(0.15, 1.0) for _ in range(200)]
    old_g = [_star(1.0, 1.0) for _ in range(200)]
    for star in young_m + old_g:
        activity.generate_activity(star)
    assert sum(star.xuv_saturated for star in young_m) > 150
    assert sum(star.xuv_saturated for star in old_g) < 20
    assert all(1.8 <= star.flare_alpha <= 2.2 for star in young_m)


def test_white_dwarfs_and_giants_do_not_flare():
    draw.set_run_seed(1)
    dwarf = _star(0.6, 3.0, "VII", 30000.0)
    giant = _star(1.5, 3.0, "III")
    for star in (dwarf, giant):
        activity.generate_activity(star)
        assert star.flare_n33_per_yr == 0.0
    assert dwarf.log_lx_lbol is None and dwarf.l_xuv_w > 0.0
    assert giant.log_lx_lbol == pytest.approx(math.log10(activity.GIANT_LX_LBOL))


def test_compact_remnants_shine_in_xuv():
    draw.set_run_seed(2)
    pulsar = NeutronStar(SystemConfig())
    assert pulsar.l_xuv_w >= pulsar.luminosity and pulsar.flare_n33_per_yr == 0.0
    hole = BlackHole(SystemConfig())
    assert hole.l_xuv_w == pytest.approx(hole.luminosity) and hole.log_lx_lbol is None


def _body(mass_earth, density, period_hours, body_type="t", is_moon=False, age=4.5):
    star = types.SimpleNamespace(age=age, mass=constants.SOLAR_MASS_TO_KG, radius=695700.0, luminosity=SUN,
                                 log_lx_lbol=math.log10(magnetism.SOLAR_LX_W / SUN), l_xuv_w=None,
                                 xuv_fluence_j=None, flare_n33_per_yr=None)
    return types.SimpleNamespace(mass=mass_earth * constants.EARTH_MASS_TO_KG, density=density,
                                 rotation_period_hours=period_hours, body_type=body_type, is_moon=is_moon,
                                 star=star)


def test_core_fractions_follow_density():
    assert magnetism.core_mass_fraction(1.0, 5.51) == pytest.approx(0.325)
    assert magnetism.core_mass_fraction(0.0123, 3.34) == 0.0  # the Moon
    assert magnetism.core_mass_fraction(0.055, 5.43) > 0.4    # Mercury


def test_a_mercury_like_planet_keeps_its_field_at_any_age():
    """GEN.213: the lifetime test used to run first and always zeroed a 0.055 Earth-mass iron world at 4.5 Gyr."""
    draw.set_run_seed(13)
    for age in (0.5, 4.5, 10.0):
        for _ in range(50):
            body = _body(0.055, 5.43, 1407.6, age=age)
            magnetism.generate_field(body)
            assert 1e18 < body.magnetic_moment_a_m2 < 3e21  # 3e19 A m2, log-normal 0.5 dex
    small_rocky = _body(0.03, 5.43, 24.0)  # under 0.05 Earth masses: too small for a dynamo
    magnetism.generate_field(small_rocky)
    assert small_rocky.magnetic_moment_a_m2 == 0.0


def test_earth_mercury_and_jupiter_get_their_fields():
    draw.set_run_seed(4)
    earths = []
    for _ in range(300):
        body = _body(1.0, 5.51, 24.0)
        body.hydrosphere = "surface ocean"
        magnetism.generate_field(body)
        earths.append(body)
    alive = [body.magnetic_moment_a_m2 for body in earths if body.magnetic_moment_a_m2]
    # The dynamo lasts 2 to 6 Gyr at one Earth mass with plate tectonics (60% of watery Earth-mass worlds)
    # and 0.8 to 2.4 Gyr under a stagnant lid: about 0.6 x 3/8 + 0.4 x 0 of them are alive at 4.5 Gyr.
    assert 40 < len(alive) < 110
    median = sorted(alive)[len(alive) // 2]
    assert 0.1 * magnetism.EARTH_MOMENT_A_M2 < median < 3 * magnetism.EARTH_MOMENT_A_M2
    jupiter = _body(318.0, 1.33, 10.0, body_type="g")
    magnetism.generate_field(jupiter)
    assert jupiter.dipole_class == "strong" and jupiter.magnetic_moment_a_m2 > 1e26
    locked = _body(1.0, 5.51, 24.0 * 11)
    locked.hydrosphere = "surface ocean"
    magnetism.generate_field(locked)
    while not locked.magnetic_moment_a_m2:
        magnetism.generate_field(locked)
    assert locked.dipole_class == "multipolar"


def test_only_a_planet_with_a_surface_ocean_can_have_a_mobile_lid():
    """GEN.177: plate tectonics needs a liquid surface ocean and half an Earth mass."""
    earth = _body(1.0, 5.51, 24.0)
    earth.hydrosphere = "surface ocean"
    assert magnetism.mobile_lid_chance(earth, 1.0) == magnetism.MOBILE_LID_CHANCE
    assert magnetism.mobile_lid_chance(earth, 8.0) == magnetism.MOBILE_LID_HEAVY_CHANCE
    assert magnetism.mobile_lid_chance(earth, 0.3) == 0.0
    for state in ("dry", "vapour", "ice", "ice-covered ocean", "hycean", None):
        venus = _body(0.8, 5.24, 24.0)
        venus.hydrosphere = state
        assert magnetism.mobile_lid_chance(venus, 0.8) == 0.0


def test_a_stagnant_lid_shortens_the_dynamo_and_a_mobile_lid_keeps_it():
    """GEN.177: the same planet at 4.5 Gyr keeps its field far more often with plate tectonics."""
    def alive(hydrosphere, mass=1.0):
        draw.set_run_seed(11)
        count = 0
        for _ in range(400):
            body = _body(mass, 5.51, 24.0)
            body.hydrosphere = hydrosphere
            magnetism.generate_field(body)
            count += bool(body.magnetic_moment_a_m2)
        return count
    mobile, stagnant = alive("surface ocean"), alive("dry")
    assert stagnant == 0  # 0.4 x (2 to 6 Gyr) is under 4.5 Gyr for every draw
    assert mobile > 60    # 0.6 x about 3 in 8, with room for the seed
    # A watery planet that is 4.5 Gyr old has a lid that decides its field; a young one has a field either way.
    young = _body(1.0, 5.51, 24.0, age=0.3)
    young.hydrosphere = "dry"
    magnetism.generate_field(young)
    assert young.magnetic_moment_a_m2 > 0


def test_europas_ocean_gives_the_induced_dipole_galileo_saw():
    """GEN.178: Europa's 100 km MgSO4 ocean in Jupiter's field turning every 11.2 h answers at about 100 nT."""
    synodic = magnetism.synodic_period_s(9.925, 85.228 / 24.0 / 365.25)
    assert synodic / 3600.0 == pytest.approx(11.23, abs=0.05)
    moment, response = magnetism.induced_moment(1.56e27, 6.711e8, synodic, 1.5608e6, 2.75 * 100e3, 1.5608e6 - 20e3)
    assert response > 0.95
    equatorial_nt = magnetism.VACUUM_PERMEABILITY / (4.0 * math.pi) * moment / 1.5608e6 ** 3 * 1e9
    assert 70 < equatorial_nt < 160
    # Callisto, 2.8x further out in the same field, is weaker by the cube of the distance and a thinner response.
    callisto, _ = magnetism.induced_moment(1.56e27, 1.8827e9, magnetism.synodic_period_s(9.925, 400.5 / 24 / 365.25),
                                           2.41e6, 1.0 * 100e3, 2.41e6 - 100e3)
    assert 0.02 * moment < callisto < 0.4 * moment


def test_induction_needs_an_ocean_a_field_and_a_changing_field():
    """GEN.178: no ocean, no parent field, a locked parent or a dry moon give none; a thin ocean gives less."""
    moment = magnetism.induced_moment
    assert moment(0.0, 6.7e8, 4.0e4, 1.5e6, 2.75e5, 1.4e6)[0] == 0.0     # no parent field
    assert moment(1.56e27, 6.7e8, None, 1.5e6, 2.75e5, 1.4e6)[0] == 0.0  # the parent turns with the moon
    assert moment(1.56e27, 6.7e8, 4.0e4, 1.5e6, 0.0, 1.4e6)[0] == 0.0    # no conducting layer
    thick = moment(1.56e27, 6.7e8, 4.0e4, 1.5e6, 2.75e5, 1.4e6)
    thin = moment(1.56e27, 6.7e8, 4.0e4, 1.5e6, 2.75e1, 1.4e6)
    assert thin[1] < thick[1] and thin[0] < thick[0]
    parent = types.SimpleNamespace(magnetic_moment_a_m2=1.56e27, rotation_period_hours=9.925, dipole_class="strong",
                                   body_type="g", is_moon=False)
    europa = types.SimpleNamespace(hydrosphere="ice-covered ocean", ocean_depth_km=100.0, ice_shell_km=20.0,
                                   ocean_class="neutral", period=85.228 / 24.0 / 365.25, distance=6.711e8 / 1.496e11,
                                   radius=1560.8)
    magnetism.generate_induced_field(europa, parent)
    assert 2e18 < europa.induced_moment_a_m2 < 8e18
    for state in ("dry", "ice", "vapour"):
        europa.hydrosphere = state
        magnetism.generate_induced_field(europa, parent)
        assert europa.induced_moment_a_m2 == 0.0


def test_earths_magnetopause_is_ten_radii():
    body = _body(1.0, 5.51, 24.0)
    body.magnetic_moment_a_m2 = magnetism.EARTH_MOMENT_A_M2
    share = magnetism.dipole_share(24.0)
    assert magnetism.magnetopause_rp(body, 1.0) == pytest.approx(9.75 * share ** (1 / 3), rel=0.02)
    assert magnetism.dipole_share(10.0) > 0.95 and magnetism.dipole_share(24.0 * 30) < 0.2


def test_generated_bodies_store_their_field_and_exposure():
    seen = set()
    for seed in range(10):
        with deterministic_entropy(seed):
            system = StarSystem(system_config=SystemConfig())
        star = system.star
        assert star.l_xuv_w is not None and star.xuv_fluence_j is not None
        for planet in system.planets:
            if not hasattr(planet, "moons"):
                continue
            for body in [planet, *planet.moons]:
                assert body.dipole_class in magnetism.DIPOLE_CLASSES
                seen.add(body.dipole_class)
                assert (body.magnetopause_rp is None) == (body.magnetic_moment_a_m2 == 0.0)
                assert body.xuv_flux_earth > 0.0 and body.xuv_exposure_index > 0.0
                assert body.flare_irradiation_index >= 0.0
            assert planet.xuv_flux_earth == pytest.approx(
                activity.xuv_flux_earth(planet.star, planet.distance))
    assert {"none", "strong"} <= seen


def test_a_generated_moon_has_an_induced_field_figure_and_a_planet_has_none():
    """GEN.178: `generate_moons` sets the induced moment (0 without an ocean) on each moon; a planet keeps None."""
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from tests.fuzz_support import deterministic_entropy

    moons = 0
    for seed in range(6):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(SystemConfig(), galactic_center_dist_ly=26000.0)
        for planet in system.planets:
            if not hasattr(planet, "moons"):
                continue
            assert planet.induced_moment_a_m2 is None
            for moon in planet.moons:
                moons += 1
                assert moon.induced_moment_a_m2 is not None and moon.induced_moment_a_m2 >= 0.0
                if moon.hydrosphere not in magnetism.INDUCTION_HYDROSPHERES:
                    assert moon.induced_moment_a_m2 == 0.0
    assert moons > 5
