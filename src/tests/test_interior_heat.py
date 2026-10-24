"""
The heat flowing out of a rocky body: radiogenic plus primordial, the Urey ratio an output (GEN.699, `planetgen.physics.interior_heat`).
"""

import math
from types import SimpleNamespace

import pytest

from planetgen import tuning
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants, interior_heat as ih, planets, radiogenic
from planetgen.util import draw

from tests.fuzz_support import deterministic_entropy

EARTH = constants.EARTH_MASS_TO_KG
EARTH_AGE_GY = 4.568
EARTH_ROCK = radiogenic.Rock(0.675, 0.325, 0.0, radiogenic.BULK_SILICATE_EARTH_AL_MASS_FRACTION)
MID = ih.Interior(1.0e12, 250.0)


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
    # the iron sinking releases 1 to 3e31 J on the Earth, a tenth of its accretion (Flasar and Birch 1973; Solomatov 2000)
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
    """Whatever the rheology within the paper's ranges: the geoneutrino total 38.2 +13.6 -12.7 TW (Agostini et al. 2019 via Nimmo et al.
    2020) and a Urey ratio between 0.1 and 0.9; the stagnant lid is under the Earth's measured 42 to 47 TW (a mobile lid, not modelled)."""
    for interior in CORNERS + [MID]:
        found = _earth(interior)
        assert 25.5 <= found.total_w / 1e12 <= 51.8
        assert 0.1 <= found.urey_ratio <= 0.9
        assert 1200.0 <= found.mantle_temperature_k <= tuning.INTERIOR_INITIAL_TEMPERATURE_CAP_K
        assert found.primordial_w > 0


def test_the_urey_ratio_is_about_the_stagnant_lids_of_guimond_et_al_2022():
    """Their present-day Urey ratios are about 0.65 for planets of 0.1 to 5 Earth masses (arXiv:2201.05636 section 3.2.1)."""
    for mass_earth, radius_km in ((0.1, 3000.0), (1.0, 6371.0), (5.0, 9300.0)):
        found = _earth(mass=mass_earth * EARTH, radius=radius_km)
        assert 0.5 <= found.urey_ratio <= 0.85


def test_a_younger_body_gives_more_heat_and_its_radiogenic_share_starts_small():
    ages = [0.02, 0.05, 0.1, 0.3, 0.6, 1.0, 2.0, 3.0, 4.5, 7.0, 10.0, 13.8]
    for interior in CORNERS + [MID]:
        found = [_earth(interior, age=age) for age in ages]
        totals = [f.total_w for f in found]
        assert totals == sorted(totals, reverse=True)
        assert all(f.primordial_w > 0 for f in found)
        # once the first gigayear is behind it the interior holds between half and four fifths radiogenic heat, whatever the rheology
        assert all(0.5 <= f.urey_ratio <= 0.85 for f in found[ages.index(2.0):])
    assert _earth(age=0.1).total_w > 10 * _earth(age=4.5).total_w
    assert _earth(age=0.05).urey_ratio < 0.2


def test_the_mantle_cools_from_the_magma_ocean_and_never_holds_more_than_its_cap():
    hot = _earth(age=0.001)
    assert hot.mantle_temperature_k == pytest.approx(tuning.INTERIOR_INITIAL_TEMPERATURE_CAP_K)
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
    seen = 0
    for seed in range(8):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(SystemConfig(), galactic_center_dist_ly=26000.0)
        for planet in system.planets:
            for body in [planet, *getattr(planet, "moons", [])]:
                if body.body_type == "a":
                    continue
                if body.body_type == "g":
                    assert all(getattr(body, name) is None for name in (*ih.BODY_FIELDS, *ih.OUTPUT_FIELDS))
                    continue
                seen += 1
                assert ih.interior_of(body) is not None
                again = planets.internal_heat_flux_w_m2(body)  # no draws: the same answer twice
                assert planets.internal_heat_flux_w_m2(body) == again == body.internal_heat_flux_w_m2
                # A small cold-start body (a tiny moon) is heated by its radiogenic decay for most of its life, so what
                # flows out can be a third of what the decay makes: the ratio is above 1 for it (TEST.353), not capped.
                assert again > 0.0 and 0.0 < body.urey_ratio < math.inf
    assert seen > 8
