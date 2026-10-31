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
from planetgen.physics import activity, constants, dynamo, giant_dynamo, interior_heat, magnetism, orbit_elements, radiogenic, spin
from planetgen.util import draw

from tests.fuzz_support import deterministic_entropy

SUN = constants.SOLAR_LUMINOSITY


def test_the_sun_today_matches_the_design_validation():
    t_sat, decay = activity.saturation_age_gy(1.0), activity.decay_exponent(1.0)
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
    assert activity.decay_exponent(0.1) == 1.0 and activity.decay_exponent(2.0) == 1.45
    assert activity.decay_exponent(0.7) == pytest.approx((1.8 + 1.45) / 2)


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


def test_a_stars_saturation_age_is_the_one_its_rotation_percentile_gives():
    """GEN.493: the age a star leaves saturation is `saturation_time` at the star's mass and its own uniformly drawn percentile."""
    from planetgen.physics import saturation_time

    assert activity.saturation_age_gy(0.3) == pytest.approx(saturation_time.saturation_time_myr(0.3, 50.0) / 1000.0)
    assert activity.saturation_age_gy(0.3, 5.0) < activity.saturation_age_gy(0.3) < activity.saturation_age_gy(0.3, 95.0)
    ages = []
    for seed in range(400):
        star = _star(1.0, 0.1)  # 100 Myr old: past the median saturation (23 Myr), inside the fast tail (297 Myr at the 95th percentile)
        draw.set_run_seed(seed)
        activity.generate_activity(star)
        ages.append(star.xuv_saturated)
    share = sum(ages) / len(ages)
    # saturated at 100 Myr when the percentile is above where the curve crosses 100 Myr (about the 75th of the 2 to 98 draw)
    assert 0.15 < share < 0.40


def test_the_stream_after_a_stars_activity_draws_does_not_depend_on_its_mass():
    """The percentile draw is made for every star, so a star's mass moves no later draw of the run."""
    after = []
    for mass in (0.15, 1.0, 3.0):
        draw.set_run_seed(17)
        activity.generate_activity(_star(mass, 2.0))
        after.append(draw.uniform(0.0, 1.0))
    assert after[0] == after[1] == after[2]


def test_the_decay_exponents_hold_the_old_table_level_at_the_reference_age():
    """AB1: with the median saturation age of Johnstone et al. 2021 the exponents give the `L_X/L_bol` at 5 Gyr the earlier table gave."""
    earlier = {0.45: 4.43e-5, 0.6: 4.74e-6, 0.8: 6.66e-7, 1.0: 2.96e-7}
    for mass, ratio in earlier.items():
        now = activity.lx_lbol_at(5.0, activity.saturated_lx_lbol(mass), activity.saturation_age_gy(mass), activity.decay_exponent(mass))
        assert now == pytest.approx(ratio, rel=0.08), mass
    assert 1.35 <= activity.decay_exponent(1.0) <= 1.5  # the Skumanich band


def _spinning(mass_sol, age_gy, temperature, period_over_nominal):
    star = _star(mass_sol, age_gy, "V", temperature)
    nominal = spin.nominal_rotation_period_hours(star.mass, temperature, age_gy)
    star.rotation_period_hours = nominal * period_over_nominal
    return star


def test_a_fast_rotator_is_more_active_and_flares_more_than_a_slow_one_of_the_same_mass_and_age():
    """GEN.807: the same draws, a sun-like star at 3 Gyr that spins in half its nominal period and one that spins in twice it."""
    made = {}
    for name, ratio in (("fast", 0.5), ("nominal", 1.0), ("slow", 2.0)):
        draw.set_run_seed(807)
        star = _spinning(1.0, 3.0, 5772.0, ratio)
        activity.generate_activity(star)
        made[name] = star
    assert made["fast"].l_xuv_w > made["nominal"].l_xuv_w > made["slow"].l_xuv_w
    assert made["fast"].flare_n33_per_yr > made["nominal"].flare_n33_per_yr > made["slow"].flare_n33_per_yr
    assert made["fast"].xuv_fluence_j > made["slow"].xuv_fluence_j
    # the unsaturated law: L_X/L_bol goes as the spin's -2.7 power (at b = 2 a half-period star is 2^2.7 times as active as it would be at the
    # same saturation age, but the decay is (t / t_sat)^-b, so the ratio between stars is (nominal period / period)^2.7 while both are unsaturated)
    ratio = 10 ** (made["fast"].log_lx_lbol - made["nominal"].log_lx_lbol)
    assert ratio == pytest.approx(2.0 ** activity.ROTATION_SLOPE, rel=0.02)
    # the flare rate follows the activity
    assert made["fast"].flare_n33_per_yr / made["nominal"].flare_n33_per_yr == pytest.approx(ratio, rel=0.02)


def test_a_fast_rotator_stays_saturated_longer_and_a_slow_one_leaves_saturation_sooner():
    saturated = {"fast": 0, "slow": 0}
    for seed in range(60):
        for name, ratio in (("fast", 0.3), ("slow", 3.0)):
            draw.set_run_seed(900 + seed)
            star = _spinning(1.0, 0.3, 5772.0, ratio)
            activity.generate_activity(star)
            saturated[name] += bool(star.xuv_saturated)
    # the saturation age is also drawn from the initial-rotation percentile (GEN.493), so the fast rotator is saturated at 0.3 Gyr in a share of the draws, not in all of them
    assert saturated["fast"] > 15 > saturated["slow"] + 10, saturated


def test_a_star_without_a_spin_or_with_a_law_that_does_not_cover_it_keeps_the_table():
    plain, hot = _star(1.0, 3.0, "V", 5772.0), _star(1.6, 1.0, "V", 8000.0)
    hot.rotation_period_hours = 20.0
    draw.set_run_seed(810)
    activity.generate_activity(plain)
    draw.set_run_seed(810)
    spun = _spinning(1.0, 3.0, 5772.0, 1.0)
    activity.generate_activity(spun)
    assert spun.log_lx_lbol == pytest.approx(plain.log_lx_lbol) and spun.flare_n33_per_yr == pytest.approx(plain.flare_n33_per_yr)
    draw.set_run_seed(811)
    activity.generate_activity(hot)
    assert hot.log_lx_lbol < -5.0
    assert activity.rotation_saturation_scale(None, 10.0, 2.0) == 1.0 and activity.rotation_saturation_scale(100.0, 0.0, 2.0) == 1.0


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
    mass = mass_earth * constants.EARTH_MASS_TO_KG
    radius_km = (3.0 * mass / (4.0 * math.pi * density * 1000.0)) ** (1.0 / 3.0) / 1000.0
    return types.SimpleNamespace(mass=mass, density=density, radius=radius_km,
                                 rotation_period_hours=period_hours, body_type=body_type, is_moon=is_moon,
                                 star=star)


def test_core_fractions_follow_density():
    assert magnetism.core_mass_fraction(1.0, 5.51) == pytest.approx(0.325)
    assert magnetism.core_mass_fraction(0.0123, 3.34) == 0.0  # the Moon
    assert magnetism.core_mass_fraction(0.055, 5.43) > 0.4    # Mercury


def _core(mass_earth, radius_km, metal, age=4.5, lid=interior_heat.STAGNANT):
    """The `core_heat.Core` of a rocky planet (`interior_heat.heat`, bulk silicate Earth's radiogenic content)."""
    k, th, u = radiogenic.BULK_SILICATE_EARTH
    content = radiogenic.Content(k * radiogenic.formation_mix("K")[1], th * radiogenic.formation_mix("Th")[1],
                                 u * radiogenic.formation_mix("U")[1], radiogenic.SOLAR_AL26_RATIO)
    rock = radiogenic.Rock(1.0 - metal, metal, 0.0, radiogenic.BULK_SILICATE_EARTH_AL_MASS_FRACTION)
    return interior_heat.heat_flow(content, interior_heat.Interior(1.0e12, 250.0), mass_earth * constants.EARTH_MASS_TO_KG, radius_km,
                                   age, rock, None, lid).core


def test_a_mercury_like_planet_keeps_its_field_at_any_age():
    """GEN.213: an iron world of 0.055 Earth masses keeps Mercury's field, whatever the core model says of a pure-iron core."""
    draw.set_run_seed(13)
    for age in (0.5, 4.5, 10.0):
        for _ in range(50):
            body = _body(0.055, 5.43, 1407.6, age=age)
            magnetism.generate_field(body, _core(0.055, 2440.0, 0.68, age=age))
            assert 1e17 < body.magnetic_moment_a_m2 < 3e21  # 3e19 A m2, log-normal 0.5 dex (four sigma each way)
    small_rocky = _body(0.03, 5.43, 24.0)  # under 0.05 Earth masses: too small for a dynamo
    magnetism.generate_field(small_rocky, None)
    assert small_rocky.magnetic_moment_a_m2 == 0.0


def test_earth_mercury_and_jupiter_get_their_fields():
    draw.set_run_seed(4)
    core = _core(1.0, 6371.0, 0.325, age=4.568)
    assert core.running
    moments = []
    for _ in range(300):
        body = _body(1.0, 5.51, 24.0, age=4.568)
        magnetism.generate_field(body, core)
        moments.append(body.magnetic_moment_a_m2)
    median = sorted(moments)[len(moments) // 2]
    assert all(moment > 0.0 for moment in moments)
    field = dynamo.field(core.dynamo_w * magnetism.CORE_POWER_SCALE, core.radius_m, core.inner_fraction * core.radius_m, core.density_kg_m3, giant_dynamo.spin_rad_s(24.0))
    assert median == pytest.approx(field.moment_a_m2, rel=0.05)           # log-normal 0.07 dex about the dynamo's moment
    assert 0.5 * magnetism.EARTH_MOMENT_A_M2 < median < 2 * magnetism.EARTH_MOMENT_A_M2
    locked = _body(1.0, 5.51, 24.0 * 11)
    magnetism.generate_field(locked, core)
    assert locked.dipole_class == "multipolar"


def _giant(spec, seed, period_hours):
    """A giant that stores the layers of `tests.test_giant_interior.interior` (the structure's radii placed) and what it radiates."""
    from planetgen.physics import giant_interior
    from tests.test_giant_interior import interior
    layers = interior(spec, seed)
    body = _body(spec["mass_earth"], 1.33, period_hours, body_type="g")
    body.radius = spec["radius_km"]
    for name, value in giant_interior.body_values(layers).items():
        setattr(body, name, value)
    return body


def test_a_giants_field_comes_from_its_convecting_layer_and_its_heat():
    """GEN.967: Jupiter-like and Saturn-like giants are strong dipoles of the right size, the ice giant is weaker, none is drawn from anchors."""
    from tests.test_giant_interior import JUPITER, SATURN
    draw.set_run_seed(7)
    jupiter = _giant(JUPITER, 1, 9.93)
    magnetism.generate_field(jupiter)
    assert 3.0e26 < jupiter.magnetic_moment_a_m2 < 6.0e27             # Jupiter's 1.5e27 within a factor of four; its Rossby number 0.13 sits on the gate (0.12), so the class may be either
    assert jupiter.dipole_class in ("strong", "multipolar", "earth-like")
    saturn = _giant(SATURN, 5, 10.66)
    magnetism.generate_field(saturn)
    assert 0.0 < saturn.magnetic_moment_a_m2 < jupiter.magnetic_moment_a_m2
    bare = _body(318.0, 1.33, 10.0, body_type="g")                                 # no stored layers: no region, no field
    magnetism.generate_field(bare)
    assert bare.magnetic_moment_a_m2 == 0.0 and bare.dipole_class == "none"


def test_a_moons_core_is_run_like_a_planets():
    """GEN.967: a moon with a running core has a field from `dynamo.field`; no Ganymede chance is drawn."""
    draw.set_run_seed(8)
    young = _core(0.025, 2630.0, 0.3, age=0.3)
    moon = _body(0.025, 1.94, 171.7, is_moon=True)
    magnetism.generate_field(moon, young)
    expected = dynamo.field(young.dynamo_w * magnetism.CORE_POWER_SCALE, young.radius_m, young.inner_fraction * young.radius_m, young.density_kg_m3,
                            giant_dynamo.spin_rad_s(171.7)).moment_a_m2
    if young.running:
        assert moon.magnetic_moment_a_m2 == pytest.approx(expected, rel=0.5) and moon.magnetic_moment_a_m2 > 0.0
    dead = _body(0.025, 1.94, 171.7, is_moon=True)
    magnetism.generate_field(dead, None)
    assert dead.magnetic_moment_a_m2 == 0.0


def test_the_core_decides_whether_a_planet_has_a_field():
    """GEN.838: a dead core (Mars) has no field at 4.5 Gyr and a young one has; a bigger planet's stronger core makes a stronger field."""
    draw.set_run_seed(5)
    mars_old, mars_young = _core(0.107, 3390.0, 0.22), _core(0.107, 3390.0, 0.22, age=0.3)
    assert not mars_old.running and mars_young.running
    old, young = _body(0.107, 3.93, 24.6), _body(0.107, 3.93, 24.6, age=0.3)
    magnetism.generate_field(old, mars_old)
    magnetism.generate_field(young, mars_young)
    assert old.magnetic_moment_a_m2 == 0.0 and old.dipole_class == "none"
    assert young.magnetic_moment_a_m2 > 0.0
    heavy = _core(2.0, 7800.0, 0.33)
    assert heavy.dipole_scale > _core(1.0, 6371.0, 0.325).dipole_scale > mars_young.dipole_scale > 0.0


def test_a_planet_with_no_core_to_read_has_no_field():
    draw.set_run_seed(6)
    body = _body(1.0, 5.51, 24.0)
    magnetism.generate_field(body, None)
    assert body.magnetic_moment_a_m2 == 0.0 and body.dipole_class == "none"


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


def test_the_lid_reaches_the_field_through_the_mantle_the_core_sits_under():
    """GEN.838 (GEN.177 before it): a stagnant lid keeps the mantle, and so the base of it the core is measured against, hotter; Earth's
    core makes its field under either lid."""
    stagnant, mobile = _core(1.0, 6371.0, 0.325, lid=interior_heat.STAGNANT), _core(1.0, 6371.0, 0.325, lid=interior_heat.MOBILE)
    assert stagnant.boundary_k > mobile.boundary_k
    draw.set_run_seed(11)
    for core in (stagnant, mobile):
        body = _body(1.0, 5.51, 24.0)
        magnetism.generate_field(body, core)
        assert body.magnetic_moment_a_m2 > 0.0


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
    assert magnetism.magnetopause_rp(body, 1.0) == pytest.approx(9.75, rel=0.02)
    body.magnetic_moment_a_m2 = magnetism.EARTH_MOMENT_A_M2 / 8.0            # the stored moment is the dipole's: an eighth of it, half the standoff
    assert magnetism.magnetopause_rp(body, 1.0) == pytest.approx(9.75 / 2.0, rel=0.02)


def test_the_stored_magnetopause_is_in_radii_of_the_body_itself():
    """GEN.843: the column's unit is planet radii (`column_units` R_planet), so a distance in metres divided by the body's own radius:
    Earth, Mars and Jupiter hold the same field and wind to a standoff that falls as the radius grows, and it agrees with the metres
    the belts work from."""
    from planetgen.physics import radiation_belts
    cases = (("Earth", 1.0, 5.51, constants.EARTH_RADIUS_KM, magnetism.EARTH_MOMENT_A_M2),
             ("Mars", 0.107, 3.93, 3389.5, magnetism.EARTH_MOMENT_A_M2),
             ("Jupiter", 317.8, 1.33, 69911.0, 1.55e27))
    stored = {}
    for name, mass_earth, density, radius_km, moment in cases:
        body = _body(mass_earth, density, 24.0)
        body.radius = radius_km
        body.magnetic_moment_a_m2 = moment
        stored[name] = magnetism.magnetopause_rp(body, 1.0)
        metres = radiation_belts.standoff_m(moment,
                                            magnetism.wind_pressure_pa(body.star, 1.0))
        assert stored[name] == pytest.approx(metres / (radius_km * 1000.0), rel=1e-9), name
    assert stored["Earth"] == pytest.approx(9.75, rel=0.02)
    assert stored["Jupiter"] * 69911.0 > stored["Earth"] * constants.EARTH_RADIUS_KM     # a Jupiter-strong field reaches farther in metres
    assert stored["Mars"] > stored["Earth"]                                                                # a smaller body has more radii in the same distance


def test_generated_bodies_store_their_field_and_exposure():
    seen = set()
    for seed in range(30):         # GEN.967: a giant's field is the dynamo's, so a strong one needs a giant that turns fast; seed 25 has one
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
            # GEN.727: the light a planet gets is that of the mean flux of its ellipse, the distance of `flux_distance_au`
            assert planet.xuv_flux_earth == pytest.approx(
                activity.xuv_flux_earth(planet.star, orbit_elements.flux_distance_au(planet.distance, planet.eccentricity)))
    assert {"none", "strong"} <= seen


def test_a_generated_moon_has_an_induced_field_figure_and_a_giant_has_none():
    """GEN.178: `generate_moons` sets the induced moment (0 without an ocean) on each moon. GEN.959: the tide pulls the moons
    of a planet locked to its star in, so the seeds that give moons are the ones with a planet that turns freely (seeds 6, 10 and 11 of the first twelve)."""
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from tests.fuzz_support import deterministic_entropy

    moons = 0
    for seed in range(24):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(SystemConfig(), galactic_center_dist_ly=26000.0)
        for planet in system.planets:
            if not hasattr(planet, "moons"):
                continue
            assert (planet.induced_moment_a_m2 is None) == (planet.body_type == "g")      # GEN.981: a rocky planet has its star's induction, a giant none
            for moon in planet.moons:
                moons += 1
                assert moon.induced_moment_a_m2 is not None and moon.induced_moment_a_m2 >= 0.0
                if moon.hydrosphere not in magnetism.INDUCTION_HYDROSPHERES and interior_heat.core_of(moon) is None:
                    assert moon.induced_moment_a_m2 == 0.0           # GEN.981: a moon with no ocean answers with its metal core, if it has one
    assert moons > 5


# --- GEN.982: the far-ultraviolet flux a star gives at a distance, relative to the quiet Sun at 1 AU ---------------------------------

def _plain_star(temperature, luminosity_solar, mass_solar, radius_solar, log_lx_lbol):
    return types.SimpleNamespace(temperature=temperature, luminosity=luminosity_solar * constants.SOLAR_LUMINOSITY, mass=mass_solar * constants.SOLAR_MASS_TO_KG,
                           radius=radius_solar * constants.SOLAR_RADIUS_M / 1.0e3, log_lx_lbol=log_lx_lbol)


def test_the_far_ultraviolet_flux_falls_with_the_square_of_the_distance_and_is_one_for_the_quiet_sun_at_one_au():
    sun_log_ratio = math.log10(activity.lx_lbol_at(activity.SUN_AGE_GY, activity.saturated_lx_lbol(1.0), activity.saturation_age_gy(1.0), activity.decay_exponent(1.0)))
    sun = _plain_star(activity.SUN_TEMPERATURE_K, 1.0, 1.0, 1.0, sun_log_ratio)
    assert activity.far_uv_flux_ratio(sun, 1.0) == pytest.approx(1.0, rel=0.05)
    assert activity.far_uv_flux_ratio(sun, 2.0) == pytest.approx(activity.far_uv_flux_ratio(sun, 1.0) / 4.0)


def test_a_more_active_star_gives_more_far_ultraviolet_at_the_same_luminosity_and_distance():
    quiet = _plain_star(3400.0, 0.03, 0.3, 0.3, -5.5)
    active = _plain_star(3400.0, 0.03, 0.3, 0.3, -3.5)
    assert activity.far_uv_flux_ratio(active, 0.1) > 10.0 * activity.far_uv_flux_ratio(quiet, 0.1)


def test_a_star_without_its_quantities_or_a_distance_has_no_ratio():
    star = _plain_star(5772.0, 1.0, 1.0, 1.0, -6.0)
    assert activity.far_uv_flux_ratio(star, 0.0) is None and activity.far_uv_flux_ratio(star, None) is None
    star.luminosity = None
    assert activity.far_uv_flux_ratio(star, 1.0) is None
