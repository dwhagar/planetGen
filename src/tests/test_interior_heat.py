"""
The heat flowing out of a rocky body: radiogenic plus primordial, the Urey ratio an output (GEN.699, `planetgen.physics.interior_heat`).
"""

import math
import time
from types import SimpleNamespace

import pytest

from planetgen import tuning
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants, giant_heat, interior_heat as ih, planets, radiogenic
from planetgen.util import draw

from tests.fuzz_support import deterministic_entropy

EARTH = constants.EARTH_MASS_TO_KG
EARTH_AGE_GY = 4.568
EARTH_ROCK = radiogenic.Rock(0.675, 0.325, 0.0, radiogenic.BULK_SILICATE_EARTH_AL_MASS_FRACTION)
MID = ih.Interior(1.0e12, 250.0)
_EARTH_G = constants.G * EARTH / (6371.0e3) ** 2
EARTH_CEILING_K = ih.mantle_ceiling_k(_EARTH_G, 6371.0e3 * (1 - 0.325 ** 0.5))


def _bulk_silicate_earth():
    k, th, u = radiogenic.BULK_SILICATE_EARTH
    return radiogenic.Content(k * radiogenic.formation_mix("K")[1], th * radiogenic.formation_mix("Th")[1],
                              u * radiogenic.formation_mix("U")[1], radiogenic.SOLAR_AL26_RATIO)


def _earth(interior=MID, age=EARTH_AGE_GY, mass=EARTH, radius=6371.0, rock=EARTH_ROCK):
    return ih.heat_flow(_bulk_silicate_earth(), interior, mass, radius, age, rock)


CORNERS = [ih.Interior(eta0, activation) for eta0 in (2.6e10, 5.3e13) for activation in (200.0, 300.0)]


def test_the_assembly_energy_of_an_earth_is_the_accretion_plus_the_iron_sinking():
    mass, radius = EARTH, 6.371e6
    uniform = ih.binding_energy_j(mass, radius, 0.0, 0.0)
    assert uniform == pytest.approx(0.6 * constants.G * mass ** 2 / radius, rel=1e-12)
    differentiated = ih.binding_energy_j(mass, radius, 0.325 * mass, 0.55 * radius)
    # the iron sinking releases 1 to 3e31 J on the Earth, a tenth of its accretion [R1388]
    assert 1.0e31 < differentiated - uniform < 3.0e31
    assert differentiated > uniform


def test_the_heat_out_is_the_radiogenic_heat_plus_the_cooling_of_the_interior():
    """(Mm cm + Mc cc) dTm/dt = radiogenic heat - heat out, so the primordial heat is minus the heat capacity times the cooling rate."""
    capacity = EARTH * (0.675 * ih.MANTLE_SPECIFIC_HEAT + 0.325 * ih.CORE_SPECIFIC_HEAT)
    before, after = _earth(age=4.4), _earth(age=4.6)
    cooling = -(after.mantle_temperature_k - before.mantle_temperature_k) / (0.2 * ih.GY_S)
    assert cooling > 0
    assert _earth(age=4.5).primordial_w == pytest.approx(capacity * cooling, rel=0.02)
    found = _earth()
    assert found.total_w == pytest.approx(found.radiogenic_w + found.primordial_w)
    assert found.urey_ratio == pytest.approx(found.radiogenic_w / found.total_w)
    assert found.flux_w_m2 == pytest.approx(found.total_w / (4 * math.pi * (6371e3) ** 2))


def test_the_radiogenic_part_is_the_stored_contents_heat():
    found = _earth()
    assert found.radiogenic_w == pytest.approx(radiogenic.heat_w(_bulk_silicate_earth(), EARTH * 0.675, EARTH_AGE_GY), rel=1e-12)
    assert found.radiogenic_w / 1e12 == pytest.approx(21.0 * 0.675 * EARTH / 4.042e24 / 1.0, rel=0.04)


def test_an_earth_twin_gives_a_heat_and_a_urey_ratio_in_the_published_ranges_without_drawing_the_ratio():
    """Whatever the rheology within the paper's ranges: the geoneutrino total 38.2 +13.6 -12.7 TW and a Urey ratio between 0.1 and 0.9; the stagnant lid is under the Earth's measured 42 to 47 TW (a mobile lid, not modelled)."""
    for interior in CORNERS + [MID]:
        found = _earth(interior)
        assert 25.5 <= found.total_w / 1e12 <= 51.8
        assert 0.1 <= found.urey_ratio <= 0.9
        assert 1200.0 <= found.mantle_temperature_k <= EARTH_CEILING_K
        assert found.primordial_w > 0


def test_the_urey_ratio_is_about_the_stagnant_lids_of_guimond_et_al_2022():
    """Their present-day Urey ratios are about 0.65 for planets of 0.1 to 5 Earth masses ([R1716] section 3.2.1)."""
    for mass_earth, radius_km in ((0.1, 3000.0), (1.0, 6371.0), (5.0, 9300.0)):
        found = _earth(mass=mass_earth * EARTH, radius=radius_km)
        assert 0.5 <= found.urey_ratio <= 0.85


def test_a_younger_body_gives_more_heat_and_its_radiogenic_share_starts_small():
    ages = [0.02, 0.05, 0.1, 0.3, 0.6, 1.0, 2.0, 3.0, 4.5, 7.0, 10.0, 13.8]
    for interior in CORNERS + [MID]:
        found = [_earth(interior, age=age) for age in ages]
        totals = [f.total_w for f in found]
        assert totals == sorted(totals, reverse=True)
        # the stiffest rheology sits at the liquidus ceiling for its first half gigayear, radiating all it keeps: no primordial surplus yet
        assert all(f.primordial_w >= 0 for f in found) and all(f.primordial_w > 0 for f in found[ages.index(1.0):])
        # once the first gigayear is behind it the interior holds between half and four fifths radiogenic heat, whatever the rheology
        assert all(0.5 <= f.urey_ratio <= 0.85 for f in found[ages.index(2.0):])
    assert _earth(age=0.1).total_w > 10 * _earth(age=4.5).total_w
    assert _earth(age=0.05).urey_ratio < 0.2


def test_the_mantle_cools_from_the_magma_ocean_and_never_holds_more_than_its_cap():
    hot = _earth(age=0.001)
    assert hot.mantle_temperature_k == pytest.approx(EARTH_CEILING_K)
    # held at the cap by radiating the 26Al heat away: the heat out is at least what the rock makes
    assert hot.total_w >= hot.radiogenic_w
    temperatures = [_earth(age=age).mantle_temperature_k for age in (0.001, 0.1, 1.0, 4.5)]
    assert temperatures == sorted(temperatures, reverse=True)


def test_a_small_body_starts_cold_a_big_one_hot():
    """E / C: a planetesimal's assembly energy warms it by a kelvin, a planet's by thousands (capped at the magma ocean)."""
    asteroid_rock = radiogenic.Rock(0.9, 0.1, 0.0, radiogenic.BULK_SILICATE_EARTH_AL_MASS_FRACTION)
    kg, radius_m = 1.0e18, 5.0e4
    energy = ih.binding_energy_j(kg, radius_m, 0.1 * kg, radius_m * math.sqrt(0.1))
    capacity = kg * (0.9 * ih.MANTLE_SPECIFIC_HEAT + 0.1 * ih.CORE_SPECIFIC_HEAT)
    assert energy / capacity < 5.0
    content = radiogenic.Content(0.0, 0.0, 0.0, 0.0)  # no heat sources at all: nothing but the start to cool
    small = ih.heat_flow(content, MID, kg, 50.0, 1.0e-6, asteroid_rock)
    assert small.mantle_temperature_k < 300.0
    big = ih.heat_flow(content, MID, EARTH, 6371.0, 1.0e-6, EARTH_ROCK)
    assert big.mantle_temperature_k > 2500.0


def test_the_flux_grows_with_mass_and_the_heat_is_a_flux_over_the_surface():
    fluxes = [_earth(mass=m * EARTH, radius=r).flux_w_m2 for m, r in ((0.1, 3000.0), (1.0, 6371.0), (5.0, 9300.0))]
    assert fluxes == sorted(fluxes)
    assert 1.2 < fluxes[2] / fluxes[1] < 3.5


def test_the_step_does_not_change_the_answer():
    reference = _earth()
    saved = ih.STEP_LIMIT_GY, ih.FIRST_STEP_GY, ih.STEP_GROWTH
    try:
        ih.STEP_LIMIT_GY, ih.FIRST_STEP_GY, ih.STEP_GROWTH = 0.002, 1.0e-6, 1.03
        fine = _earth()
    finally:
        ih.STEP_LIMIT_GY, ih.FIRST_STEP_GY, ih.STEP_GROWTH = saved
    assert reference.total_w == pytest.approx(fine.total_w, rel=1e-4)


def test_a_small_cold_body_is_integrated_stably():
    """The time step follows the cooling rate: a body of 50 km would blow up an explicit step sized for an Earth."""
    asteroid_rock = radiogenic.Rock(0.9, 0.1, 0.0, radiogenic.BULK_SILICATE_EARTH_AL_MASS_FRACTION)
    for age in (0.01, 1.0, 13.8):
        found = ih.heat_flow(_bulk_silicate_earth(), ih.Interior(2.6e10, 300.0), 1.0e18, 50.0, age, asteroid_rock)
        assert 273.0 <= found.mantle_temperature_k <= 3000.0 and found.total_w >= 0.0 and math.isfinite(found.urey_ratio)


def test_a_body_with_no_silicate_has_no_heat():
    found = ih.heat_flow(_bulk_silicate_earth(), MID, EARTH, 6371.0, 4.5, radiogenic.Rock(0.0, 0.0, 1.0, 0.0))
    assert found.total_w == 0.0 and found.flux_w_m2 == 0.0 and found.urey_ratio == 0.0


def test_water_over_the_rock_shrinks_the_rock_the_heat_goes_through():
    wet_rock = radiogenic.Rock(0.45, 0.15, 0.4, radiogenic.BULK_SILICATE_EARTH_AL_MASS_FRACTION)
    assert ih.dry_radius_km(3.0 * EARTH, 9000.0, EARTH_ROCK) == 9000.0
    dry = ih.dry_radius_km(3.0 * EARTH, 9000.0, wet_rock)
    assert 3000.0 < dry < 9000.0
    # outside both fits: the dry share's cube root of the radius
    tiny = ih.dry_radius_km(1.0e18, 400.0, radiogenic.Rock(0.3, 0.1, 0.6, 0.0238))
    assert tiny == pytest.approx(400.0 * 0.4 ** (1 / 3))
    assert _earth(mass=3.0 * EARTH, radius=9000.0, rock=wet_rock).total_w > 0.0


def test_an_age_before_formation_is_refused():
    with pytest.raises(ih.InteriorHeatError):
        _earth(age=-1.0)


def test_the_draws_are_two_in_their_ranges_and_the_stream_stays_in_step():
    stream_a, stream_b = draw.Stream("interior-a"), draw.Stream("interior-a")
    drawn = [ih.draw_interior(stream_a) for _ in range(200)]
    for _ in range(200):
        stream_b.uniform(0, 1)
        stream_b.uniform(0, 1)
    assert stream_a.random() == stream_b.random()
    low, high = tuning.INTERIOR_VISCOSITY_PREFACTOR_RANGE_PA_S
    assert all(low <= found.eta0_pa_s <= high for found in drawn)
    low, high = tuning.INTERIOR_ACTIVATION_RANGE_KJ_MOL
    assert all(low <= found.activation_kj_mol <= high for found in drawn)
    assert max(f.activation_kj_mol for f in drawn) - min(f.activation_kj_mol for f in drawn) > 80.0


def test_a_body_with_no_stored_interior_has_no_heat_and_a_giant_none_at_all():
    old = SimpleNamespace(star=SimpleNamespace(age=4.5), mass=EARTH, radius=6371.0, c_k0_ppm=None, c_th0_ppm=None, c_u0_ppm=None,
                          al26_ratio=None, mantle_eta0_pa_s=None, mantle_activation_kj_mol=None)
    assert ih.update_heat(old) is None and old.urey_ratio is None and old.internal_heat_flux_w_m2 is None
    assert planets.internal_heat_flux_w_m2(old) == 0.0
    giant = SimpleNamespace(body_type="g", star=SimpleNamespace(age=4.5), c_k0_ppm=1.0, mantle_eta0_pa_s=1.0)
    ih.generate_interior(giant)
    assert all(getattr(giant, name) is None for name in (*ih.BODY_FIELDS, *ih.OUTPUT_FIELDS))


def test_every_generated_rocky_body_stores_its_draws_and_an_urey_ratio_that_is_an_output():
    seen = seen_convecting = 0
    for seed in range(8):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(SystemConfig(), galactic_center_dist_ly=26000.0)
        for planet in system.planets:
            for body in [planet, *getattr(planet, "moons", [])]:
                if body.body_type == "a":
                    continue
                if body.body_type == "g":
                    # a giant has no rock to draw and no Urey ratio or tide; its own heat is what it radiates (GEN.810)
                    assert all(getattr(body, name) is None for name in (*ih.BODY_FIELDS, "urey_ratio", "tidal_flux_w_m2"))
                    planets.internal_heat_flux_w_m2(body)         # runs `update_heat` again, at the star's age now
                    assert body.internal_heat_flux_w_m2 == giant_heat.of(body).flux_w_m2 > 0.0
                    continue
                seen += 1
                assert ih.interior_of(body) is not None
                again = planets.internal_heat_flux_w_m2(body)  # no draws: the same answer twice
                assert planets.internal_heat_flux_w_m2(body) == again == body.internal_heat_flux_w_m2
                # the ratio is an output with no cap, and a mantle that is not convecting loses nothing
                assert again >= 0.0 and body.urey_ratio >= 0.0
                seen_convecting += again > 0.0
    assert seen > 8 and seen_convecting > 0.8 * seen


def _earth_lid(interior, lid, age=EARTH_AGE_GY):
    return ih.heat_flow(_bulk_silicate_earth(), interior, EARTH, 6371.0, age, EARTH_ROCK, None, lid)


def test_a_mobile_lid_cools_the_mantle_further_and_spends_its_heat_sooner():
    """GEN.721: at one mantle temperature a mobile lid loses ten times the heat of a stagnant one, so the mantle ends up cooler (an Earth twin
    near the Earth's potential temperature of 1,600 K, not the stagnant lid's 1,950 K) and the heat is spent sooner."""
    for interior in CORNERS + [MID]:
        stagnant, mobile = _earth_lid(interior, ih.STAGNANT), _earth_lid(interior, ih.MOBILE)
        assert mobile.mantle_temperature_k < stagnant.mantle_temperature_k
        assert mobile.urey_ratio > stagnant.urey_ratio
        assert 22.0 <= mobile.total_w / 1e12 <= 36.0
    assert _earth_lid(MID, ih.MOBILE, age=0.1).total_w > _earth_lid(MID, ih.STAGNANT, age=0.1).total_w  # early, at the middle rheology
    assert _earth_lid(MID, ih.MOBILE).total_w < _earth_lid(MID, ih.STAGNANT).total_w  # late
    assert 1300.0 <= _earth_lid(MID, ih.MOBILE).mantle_temperature_k <= 1700.0


def test_the_default_lid_is_the_stagnant_one_and_a_body_reads_its_stored_regime():
    assert _earth_lid(MID, ih.STAGNANT).total_w == _earth().total_w
    body = SimpleNamespace(lid_regime=ih.MOBILE)
    assert ih.lid_of(body) == ih.MOBILE
    assert ih.lid_of(SimpleNamespace(lid_regime=None)) == ih.STAGNANT == ih.lid_of(SimpleNamespace())


def _watery(mass_earth, hydrosphere):
    return SimpleNamespace(mass=mass_earth * EARTH, hydrosphere=hydrosphere, c_k0_ppm=1.0, c_th0_ppm=1.0, c_u0_ppm=1.0, al26_ratio=1.0)


def test_the_regime_is_one_draw_at_the_chance_of_the_old_coin_and_a_giant_has_none():
    from planetgen.physics import magnetism
    draw.set_run_seed(5)
    bodies = [_watery(1.0, "surface ocean") for _ in range(600)]
    for body in bodies:
        ih.generate_lid(body)
    mobile = sum(body.lid_regime == ih.MOBILE for body in bodies)
    assert abs(mobile / 600 - magnetism.MOBILE_LID_CHANCE) < 0.07
    for hydrosphere, mass in (("dry", 1.0), ("ice", 1.0), ("surface ocean", 0.3)):
        dry = _watery(mass, hydrosphere)
        ih.generate_lid(dry)
        assert dry.lid_regime == ih.STAGNANT
    giant = SimpleNamespace(mass=100 * EARTH, hydrosphere=None)
    ih.generate_lid(giant)
    assert giant.lid_regime is None
    draw.set_run_seed(5)
    first = draw.random()
    draw.set_run_seed(5)
    ih.generate_lid(_watery(1.0, "surface ocean"))
    assert draw.random() != first  # exactly one draw was taken
    draw.set_run_seed(5)
    draw.random()
    after_one = draw.random()
    draw.set_run_seed(5)
    ih.generate_lid(_watery(1.0, "surface ocean"))
    assert draw.random() == after_one


def test_the_liquidus_is_manskes_fit_and_the_ceiling_follows_the_pressure_in_the_mantle():
    """[R2857: Table 1]: 1,874 K at the surface, 4,717 K at the Earth's core-mantle boundary (135 GPa). The Earth's mantle
    ceiling is the liquidus at about 47 GPa, a Moon-sized body's is near the surface liquidus."""
    assert ih.liquidus_k(0.0) == pytest.approx(1874.0, abs=3.0)
    assert ih.liquidus_k(135.0) == pytest.approx(4717.0, abs=10.0)
    assert 3000.0 <= EARTH_CEILING_K <= 3300.0
    moon_g = constants.G * 7.3e22 / (1737e3) ** 2
    assert 1900.0 <= ih.mantle_ceiling_k(moon_g, 1737e3 * 0.7) <= 2000.0
    assert ih.mantle_ceiling_k(2 * _EARTH_G, 6371.0e3) > EARTH_CEILING_K


def test_the_surface_temperature_a_mobile_lid_measures_its_drop_against_is_the_bodys_own():
    """A hotter surface leaves a smaller drop, so a mobile lid loses less heat and its mantle stays warmer; a stagnant lid's boundary layer
    does not read the surface at all but its birth temperature does."""
    cold = ih.heat_flow(_bulk_silicate_earth(), MID, EARTH, 6371.0, EARTH_AGE_GY, EARTH_ROCK, None, ih.MOBILE, 200.0)
    hot = ih.heat_flow(_bulk_silicate_earth(), MID, EARTH, 6371.0, EARTH_AGE_GY, EARTH_ROCK, None, ih.MOBILE, 737.0)
    assert hot.mantle_temperature_k > cold.mantle_temperature_k
    stagnant = [ih.heat_flow(_bulk_silicate_earth(), MID, EARTH, 6371.0, EARTH_AGE_GY, EARTH_ROCK, None, ih.STAGNANT, ts) for ts in (200.0, 737.0)]
    assert stagnant[0].mantle_temperature_k == pytest.approx(stagnant[1].mantle_temperature_k, rel=0.01)


def test_the_stored_heat_is_the_heat_at_the_stars_final_age_not_the_provisional_one():
    """GEN.763: the star's age only rises after its planets are made, so a heat stored at the provisional age is wrong for the finished
    system (many bodies sat at the ceiling with a Urey ratio of 1.0, and a small moon showed 4.8). The life pass runs the model again."""
    seen = 0
    for seed in (17, 12):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(SystemConfig(), galactic_center_dist_ly=26000.0)
        for planet in system.planets:
            for body in [planet, *getattr(planet, "moons", [])]:
                if body.body_type in ("a", "g"):
                    continue
                stored = (body.urey_ratio, body.internal_heat_flux_w_m2)
                ih.update_heat(body)
                assert stored == (body.urey_ratio, body.internal_heat_flux_w_m2)
                seen += 1
    assert seen > 10


SWEEP_MASSES_EARTH = (0.0004, 0.01, 0.1, 0.5, 1.0, 3.0, 10.0)


def test_the_urey_ratio_of_a_bulk_silicate_earth_body_of_any_size_sits_between_the_published_extremes_once_it_is_a_gigayear_old():
    """GEN.763: a table of small and large bodies (design note, section 3.2). The bounds of the sweep over 0.0004 to 10 Earth masses, 1, 4.5
    and 10 Gy, both lids and the four rheology corners: a Urey ratio from 0.3 (a large stagnant lid, whose heat is mostly the
    cooling) to 1.0 (the mantle held at its ceiling, losing exactly what is made). [R1716] have about 0.65 for 0.1 to 5
    Earth masses; for the Earth the geochemical estimate is about 0.3 and the geophysical view close to one."""
    found = []
    for mass_earth in SWEEP_MASSES_EARTH:
        for lid in (ih.STAGNANT, ih.MOBILE):
            for interior in CORNERS + [MID]:
                for age in (1.0, 4.5, 10.0):
                    heat = ih.heat_flow(_bulk_silicate_earth(), interior, mass_earth * EARTH, 6371.0 * mass_earth ** 0.3, age, EARTH_ROCK, None, lid, 273.0)
                    assert heat.total_w > 0.0
                    found.append((heat.urey_ratio, mass_earth, lid, age))
    assert 0.3 <= min(f[0] for f in found) and max(f[0] for f in found) <= 1.0 + 1e-9
    # the larger the body the more of its heat is the cooling: at 4.5 Gy the median falls with the mass
    medians = []
    for mass_earth in (0.0004, 1.0, 10.0):
        ratios = sorted(f[0] for f in found if f[1] == mass_earth and f[3] == 4.5)
        medians.append(ratios[len(ratios) // 2])
    assert medians[0] > medians[1] > medians[2]
    assert 0.55 <= medians[0] and medians[2] <= 0.8


def test_a_generated_body_has_a_urey_ratio_above_one_only_while_its_mantle_is_still_warming():
    """A mantle that is warming loses less than it makes, so the ratio passes 1 only then (a moon of 0.0008 Earth masses at 3 Gy gives
    1.25), and a body past 8 Gy has settled. Over the sampled systems every body over 5e22 kg of a host past 2 Gy lies between 0.2 and 10 (a small moon whose drawn
    viscosity is near the top of the range, 3e13 Pa s and 299 kJ/mol, is still heating up at 4.7 Gy and gives 7.8, with or without its tide). A small body
    (GEN.476: most stars are born with little aluminium-26, so a small body never had the early heat) can still be warming at any age: its ratio is only positive and finite; so is a body whose surface is so cold (a moon of 570 km at 36 K) that its mantle only starts to convect late."""
    seen = 0
    for seed in range(12):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(SystemConfig(), galactic_center_dist_ly=26000.0)
        for planet in system.planets:
            for body in [planet, *getattr(planet, "moons", [])]:
                if body.body_type in ("a", "g") or not body.internal_heat_flux_w_m2 or body.star.age < 2.0:
                    continue
                if body.mass < 5.0e22:
                    assert 0.0 < body.urey_ratio < float("inf")
                    continue
                seen += 1
                assert 0.2 <= body.urey_ratio <= 10.0
                if body.star.age >= 8.0 and ih.surface_of(body) > 100.0:
                    assert body.urey_ratio <= 1.0 + 1e-6
    assert seen > 10


@pytest.mark.timeout(60)
def test_a_hot_surface_never_collapses_the_step_below_it():
    """GEN.866: a 4.7e24 kg wet planet whose surface is 2304 K (an ecosphere orbit of an F5V star, the greenhouse run away) started at the
    surface temperature and a step overshot below zero kelvin; the viscosity law on that temperature gave a flux of 1e59 W/m2 and a step of 1e-65 Gy, so
    the draw never returned. The mantle stays between the surface and the ceiling, and the run takes well under a second."""
    terms = ((4.06712528459501e-12, 0.5554063946794433), (1.1519457618003578e-12, 0.04951051289713895),
             (1.9827998614583646e-12, 0.15530969763834757), (3.744530060985561e-12, 0.9845840632953767),
             (4.075177246711262e-07, 966.7324694002026))
    rock = radiogenic.Rock(0.7545618533566459, 0.0, 0.2454381466433541, 0.018143235516770093)
    surface_k = 2303.74544273136
    started = time.monotonic()
    found = ih.heat(terms, ih.Interior(382324606739.6106, 218.6088709840769), 4.666239116198571e24, 7011.034245444521,
                    4.086548725622917, rock, None, ih.STAGNANT, surface_k)
    assert time.monotonic() - started < 20.0
    assert surface_k <= found.mantle_temperature_k
    assert math.isfinite(found.total_w) and math.isfinite(found.released_j) and found.total_w > 0.0


def _earth_under(surface_k, lid, forcing=None):
    return ih.heat_flow(_bulk_silicate_earth(), MID, EARTH, 6371.0, EARTH_AGE_GY, EARTH_ROCK, forcing, lid, surface_k)


@pytest.mark.parametrize("lid", [ih.STAGNANT, ih.MOBILE])
def test_a_surface_above_the_mantles_liquidus_leaves_a_mantle_at_the_surface_which_makes_its_heat_and_ends(lid):
    """GEN.856: the integrator ran the mantle's temperature down past zero and never ended (a stagnant lid), or divided by a zero drop
    (a mobile one), for a surface hotter than the liquidus. A mantle can be no cooler than its surface: it sits at it and the melt
    carries out what the silicate made."""
    surface_k = EARTH_CEILING_K + 300.0
    found = _earth_under(surface_k, lid)
    assert found.mantle_temperature_k == surface_k
    assert found.total_w == pytest.approx(found.radiogenic_w) and found.primordial_w == 0.0 and found.tidal_w == 0.0
    assert found.flux_w_m2 == pytest.approx(found.total_w / (4 * math.pi * (6371e3) ** 2))
    cool = _earth_lid(MID, lid)
    assert found.radiogenic_w == pytest.approx(cool.radiogenic_w)                 # the same decay heat as a cool surface
    assert 0.0 < found.released_j < cool.released_j * 10 and found.released_j > found.radiogenic_w * ih.GY_S


def test_a_surface_exactly_at_the_liquidus_is_that_body_too():
    assert _earth_under(EARTH_CEILING_K, ih.STAGNANT).mantle_temperature_k == EARTH_CEILING_K


def test_a_molten_body_under_a_hot_surface_is_flexed_as_a_molten_one():
    """A magma ocean has the fluid Love number and almost no loss (GEN.955): the tide it feels is the rheology's at its surface's temperature."""
    from planetgen.physics import tidal_heating
    forcing = tidal_heating.Forcing(1.9e27, 4.2e8, 0.05)
    surface_k = EARTH_CEILING_K + 300.0
    found = _earth_under(surface_k, ih.STAGNANT, forcing)
    mantle = ih.mantle_of(EARTH, 6371.0, EARTH_ROCK)
    assert found.tidal_w == pytest.approx(ih._tide_w(forcing, mantle, surface_k, MID)) and found.tidal_w > 0.0
    assert found.tidal_w < 1.0e-2 * tidal_heating.peale_power_w(1.5, forcing, mantle.dry_m)
    assert found.total_w == pytest.approx(found.radiogenic_w + found.tidal_w) and found.primordial_w == 0.0


def test_a_mobile_lid_with_no_drop_carries_no_heat():
    assert ih._flux_w_m2(300.0, MID, 9.8, 3.0e6, ih.MOBILE, 300.0) == 0.0
    assert ih._flux_w_m2(250.0, MID, 9.8, 3.0e6, ih.MOBILE, 300.0) == 0.0
    assert ih._flux_w_m2(1500.0, MID, 9.8, 3.0e6, ih.MOBILE, 300.0) > 0.0
