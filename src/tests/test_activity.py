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


def test_earth_mercury_and_jupiter_get_their_fields():
    draw.set_run_seed(4)
    earths = []
    for _ in range(300):
        body = _body(1.0, 5.51, 24.0)
        magnetism.generate_field(body)
        earths.append(body)
    alive = [body.magnetic_moment_a_m2 for body in earths if body.magnetic_moment_a_m2]
    # The dynamo lasts 2 to 6 Gyr at one Earth mass, so at 4.5 Gyr about 3 in 8 are alive.
    assert 80 < len(alive) < 150
    median = sorted(alive)[len(alive) // 2]
    assert 0.1 * magnetism.EARTH_MOMENT_A_M2 < median < 3 * magnetism.EARTH_MOMENT_A_M2
    jupiter = _body(318.0, 1.33, 10.0, body_type="g")
    magnetism.generate_field(jupiter)
    assert jupiter.dipole_class == "strong" and jupiter.magnetic_moment_a_m2 > 1e26
    locked = _body(1.0, 5.51, 24.0 * 11)
    magnetism.generate_field(locked)
    while not locked.magnetic_moment_a_m2:
        magnetism.generate_field(locked)
    assert locked.dipole_class == "multipolar"


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
