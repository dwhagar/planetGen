# tests/test_radiogenic.py

"""
GEN.459: the one radiogenic module. The isotope-wise sum, a chondrite and a bulk silicate case against their calibration targets
(test targets, not presets), that the heat and the ground dose of one body read the same stored concentrations, that the old
constants are gone, and the draws (a star's offset, the potassium depletion, the Urey ratio) in their ranges.
"""

import math
import pathlib

import pytest

from planetgen import tuning
from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import elements, radiation, radiogenic, rogue_surface, surface_composition
from planetgen.util import draw

from tests.fuzz_support import deterministic_entropy

EARTH_AGE_GY = 4.568


def _content(k, th, u, al26=0.0):
    """A body whose silicate has these K, Th, U today (ppm), as the formation `Content` they come from."""
    return radiogenic.Content(k * radiogenic.formation_mix("K")[1], th * radiogenic.formation_mix("Th")[1],
                              u * radiogenic.formation_mix("U")[1], al26)


def test_the_heat_is_the_isotope_wise_sum_of_the_four_parents():
    content = _content(260.0, 0.0834, 0.0218, al26=radiogenic.SOLAR_AL26_RATIO)
    masses = radiogenic.isotope_mass_fractions(content, EARTH_AGE_GY)
    assert set(masses) >= set(radiogenic.HEAT_PARENTS)
    by_hand = sum(masses[name] * elements.heat_per_kg_of_isotope(name) for name in radiogenic.HEAT_PARENTS)
    assert radiogenic.heat_w_per_kg(content, EARTH_AGE_GY) == pytest.approx(by_hand, rel=1e-6)
    assert by_hand > 0


def test_a_bulk_silicate_earth_gives_the_calibration_heat():
    """[R3416] edition: 21.0 TW for 4.042e24 kg (radiation-model.md section 6), within 3 percent; the older 19.8 TW
    (McDonough and Sun) and the repository's former 29.8 TW (Turcotte and Schubert) are outside what the model is calibrated to."""
    heat_tw = radiogenic.heat_w(_content(*radiogenic.BULK_SILICATE_EARTH), 4.042e24, EARTH_AGE_GY) / 1e12
    assert heat_tw == pytest.approx(21.0, rel=0.03)


def test_a_chondrite_gives_the_calibration_heat():
    """CI chondrite K 546, Th 0.030, U 0.0081 ppm: 3.5 pW/kg today (part E 2.8b), within 5 percent."""
    heat = radiogenic.heat_w_per_kg(_content(546.0, 0.030, 0.0081), EARTH_AGE_GY)
    assert heat * 1e12 == pytest.approx(3.5, rel=0.05)


def test_the_heat_falls_with_age_as_the_calibration_says():
    content = _content(*radiogenic.BULK_SILICATE_EARTH)
    now = radiogenic.heat_w_per_kg(content, EARTH_AGE_GY)
    # 3.5 times higher 4 Gyr ago ([R3293]; the lane's 3.67 at 0.567 Gyr)
    assert radiogenic.heat_w_per_kg(content, 0.567) / now == pytest.approx(3.67, rel=0.05)
    ages = [0.1, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0]
    heats = [radiogenic.heat_w_per_kg(content, age) for age in ages]
    assert heats == sorted(heats, reverse=True)
    # potassium dominates early, thorium late (part E 2.8c)
    early = radiogenic.isotope_mass_fractions(content, 0.5)
    late = radiogenic.isotope_mass_fractions(content, 10.0)
    assert early["K-40"] * elements.heat_per_kg_of_isotope("K-40") > early["Th-232"] * elements.heat_per_kg_of_isotope("Th-232")
    assert late["K-40"] * elements.heat_per_kg_of_isotope("K-40") < 0.05 * radiogenic.heat_w_per_kg(content, 10.0)


def test_aluminium_26_heats_a_young_body_and_nothing_older_than_a_few_million_years():
    with_al = _content(*radiogenic.BULK_SILICATE_EARTH, al26=radiogenic.SOLAR_AL26_RATIO)
    without = _content(*radiogenic.BULK_SILICATE_EARTH, al26=0.0)
    # part E 2.8f: about 4.3e-7 W/kg at formation for the bulk silicate Al and the Solar System ratio
    young = radiogenic.heat_w_per_kg(with_al, 0.0) - radiogenic.heat_w_per_kg(without, 0.0)
    assert young == pytest.approx(4.27e-7, rel=0.1)
    assert radiogenic.heat_w_per_kg(with_al, 0.1) == pytest.approx(radiogenic.heat_w_per_kg(without, 0.1), rel=1e-12)


def test_the_surface_gamma_dose_reaches_the_calibration_values():
    primitive = _content(*radiogenic.BULK_SILICATE_EARTH)
    assert radiogenic.gamma_ngy_per_hour(primitive, EARTH_AGE_GY) == pytest.approx(0.67, rel=0.05)
    felsic_ppm = surface_composition.RESERVOIRS["felsic"]
    felsic = radiogenic.surface_content(primitive, (felsic_ppm["k_ppm"], felsic_ppm["th_ppm"], felsic_ppm["u_ppm"]))
    assert radiogenic.gamma_ngy_per_hour(felsic, EARTH_AGE_GY) == pytest.approx(71.5, rel=0.05)
    # the UNSCEAR soil means, 33, 45 and 420 Bq/kg, give 60 nGy/h; the activities come from the specific activities of GEN.438
    soil = _content(420.0 / 314.0 * 1e4, 45.0 / 4.04, 33.0 / 12.35)
    assert radiogenic.gamma_ngy_per_hour(soil, EARTH_AGE_GY) == pytest.approx(60.0, rel=0.06)
    # ordering: the older the rock, the less K-40 and U-235 it has; the same rock has less gamma than a young one
    assert radiogenic.gamma_ngy_per_hour(soil, 1.0) > radiogenic.gamma_ngy_per_hour(soil, EARTH_AGE_GY)
    # scaling: twice the concentrations, twice the dose
    twice = radiogenic.Content(2 * soil.k0_ppm, 2 * soil.th0_ppm, 2 * soil.u0_ppm, 0.0)
    assert radiogenic.gamma_ngy_per_hour(twice, 3.0) == pytest.approx(2 * radiogenic.gamma_ngy_per_hour(soil, 3.0), rel=1e-12)
    # ultramafic below basaltic below felsic
    reservoirs = [surface_composition.RESERVOIRS[name] for name in ("ultramafic", "basaltic", "felsic")]
    doses = [radiogenic.gamma_ngy_per_hour(radiogenic.surface_content(primitive, (r["k_ppm"], r["th_ppm"], r["u_ppm"])), EARTH_AGE_GY)
             for r in reservoirs]
    assert doses == sorted(doses)


def test_the_earth_soil_gives_the_world_mean_radon_and_ground_totals():
    """Radon J 0.016 to 0.026 Bq/m2/s and 10 Bq/m3 outdoors (UNSCEAR), 0.49 mSv/yr (radiation-model.md section 11 world A)."""
    soil = _content(420.0 / 314.0 * 1e4, 45.0 / 4.04, 33.0 / 12.35)
    flux = radiation.radon_flux_bq_m2_s(soil, EARTH_AGE_GY)
    assert 0.014 < flux < 0.03
    assert radiation.radon_dose_msv_yr(soil, EARTH_AGE_GY, 8.5) == pytest.approx(0.488, rel=0.1)
    gamma = radiation.gamma_dose_msv_yr(soil, EARTH_AGE_GY)
    assert gamma == pytest.approx(0.368, rel=0.07)
    # rebuilt as UNSCEAR does it (0.2 outdoors, indoor gamma 1.4 outdoors, 40 Bq/m3 indoors, 0.29 ingestion, GCR 0.34)
    external = gamma * (0.2 + 0.8 * 1.4)
    radon = 40 * 0.4 * 9e-6 * 0.8 * 8766 + radiation.radon_dose_msv_yr(soil, EARTH_AGE_GY, 8.5) * 0.2
    assert 0.34 + external + radon + 0.29 == pytest.approx(2.4, rel=0.1)


def test_the_old_constants_and_the_separate_draw_are_gone():
    for name in ("ROGUE_RADIOGENIC_ISOTOPES", "ROGUE_RADIOGENIC_ABUNDANCE_RANGE", "ROGUE_UREY_RATIO", "ROGUE_MANTLE_MASS_FRACTION"):
        assert not hasattr(tuning, name), name
    assert not hasattr(rogue_surface, "radiogenic_flux_w_m2")
    for name in ("GROUND_MSV_YR", "RADON_MSV_YR", "RADON_MIN_PRESSURE_PA"):
        assert not hasattr(radiation, name), name
    source = pathlib.Path(radiation.__file__).read_text(encoding="utf-8")
    assert "ROGUE_RADIOGENIC_ABUNDANCE_RANGE" not in source
    assert "abundance = " not in source


def test_a_stars_offset_is_its_own_and_stays_inside_the_truncation():
    first = radiogenic.host_offset_dex(0.0, 4.5, "star-a")
    assert radiogenic.host_offset_dex(0.0, 4.5, "star-a") == first
    offsets = [radiogenic.host_offset_dex(fe_h, age, f"star-{i}")
               for i, (fe_h, age) in enumerate([(0.0, 4.5), (-0.5, 9.0), (0.3, 2.0), (-1.2, 11.0), (0.1, 6.0)] * 20)]
    limit = tuning.RADIOGENIC_OFFSET_LIMIT_DEX
    assert all(-limit <= value <= limit for value in offsets)
    assert max(offsets) - min(offsets) > 0.05  # stars differ


def test_the_truncated_offset_has_the_stated_scatter():
    values = [radiogenic.truncated_offset_dex((i + 0.5) / 2000) for i in range(2000)]
    assert min(values) >= -0.5 and max(values) <= 0.5
    mean = sum(values) / len(values)
    sigma = math.sqrt(sum((v - mean) ** 2 for v in values) / len(values))
    assert abs(mean) < 1e-9
    assert 0.17 < sigma < 0.22  # 0.22 dex truncated at 0.5 (2.3 sigma) is 0.20 dex


def test_the_draws_stay_in_their_ranges_and_always_use_the_same_number_of_draws():
    stream_a, stream_b = draw.Stream("radiogenic-a"), draw.Stream("radiogenic-a")
    contents = []
    for i in range(300):
        contents.append(radiogenic.draw_content(0.2 if i % 2 else None, rng=stream_a))
        radiogenic.draw_content(-0.3 if i % 3 else None, rng=stream_b)
    # whatever the branch, the streams stay in step: the next draw is the same
    assert stream_a.random() == stream_b.random()
    ratios = [c.k0_ppm / (radiogenic.formation_mix("K")[1]) / (c.u0_ppm / radiogenic.formation_mix("U")[1]) / radiogenic.CI_CHONDRITE_K_OVER_U
              for c in contents]
    assert all(tuning.RADIOGENIC_K_DEPLETION_RANGE[0] * 0.999 <= r <= tuning.RADIOGENIC_K_DEPLETION_RANGE[1] * 1.001 for r in ratios)
    assert all(c.al26_ratio >= 0.0 for c in contents)
    assert min(r for r in ratios) < 0.1 and max(r for r in ratios) > 0.5  # the depletion spans the range


def test_an_enriched_silicate_has_more_thorium_and_uranium():
    plain = radiogenic.draw_content(0.0, rng=draw.Stream("e"))
    rich = radiogenic.draw_content(0.0, rng=draw.Stream("e"), refractory_enrichment=2.0)
    assert rich.th0_ppm == pytest.approx(2 * plain.th0_ppm)
    assert rich.u0_ppm == pytest.approx(2 * plain.u0_ppm)
    assert rich.al26_ratio == plain.al26_ratio


def _bodies(system):
    for planet in system.planets:
        if planet.body_type == "a":
            continue
        yield planet
        yield from planet.moons


def test_every_generated_body_stores_one_set_of_concentrations_that_its_heat_and_dose_read():
    seen = 0
    for seed in range(10):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(SystemConfig(), galactic_center_dist_ly=26000.0)
        for body in _bodies(system):
            if body.body_type == "g":
                assert all(getattr(body, name) is None for name in radiogenic.BODY_FIELDS)
                continue
            seen += 1
            content = radiogenic.content_of(body)
            assert content is not None
            assert content.th0_ppm > 0 and content.u0_ppm > 0 and content.k0_ppm > 0
            from planetgen.physics import interior_heat
            from planetgen.physics import planets as planet_physics
            rock = radiogenic.rock_of(body)
            from planetgen.physics import tidal_heating
            forcing = None if tidal_heating.is_icy(body) else tidal_heating.forcing_of(body)  # GEN.737: the tide of the stored orbit
            found = interior_heat.heat_flow(content, interior_heat.interior_of(body), body.mass, body.radius, body.star.age, rock, forcing,
                                            interior_heat.lid_of(body), interior_heat.surface_of(body))
            assert planet_physics.internal_heat_flux_w_m2(body) == pytest.approx(found.own_flux_w_m2)
            assert body.internal_heat_flux_w_m2 == pytest.approx(found.own_flux_w_m2) and body.urey_ratio == pytest.approx(found.urey_ratio)
            # the heat of the interior is the radiogenic heat of the stored concentrations and the cooling of the rest
            silicate_kg = body.mass * rock.silicate
            assert found.radiogenic_w == pytest.approx(radiogenic.heat_w(content, silicate_kg, body.star.age, rock.aluminium))
            assert found.total_w == pytest.approx(found.radiogenic_w + found.primordial_w + found.tidal_w)
            assert body.dose_ground_msv_yr > 0
    assert seen > 8


def test_a_planets_siblings_share_their_stars_offset():
    for seed in range(8):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(SystemConfig(), galactic_center_dist_ly=26000.0)
        offsets = set()
        for body in _bodies(system):
            if body.body_type == "g":
                continue
            offset = radiogenic.star_offset_dex(body.star)
            assert offset is not None
            offsets.add(offset)
        assert len(offsets) <= 2  # a binary's planets may hang on either star; otherwise one


def test_a_stored_body_reads_back_its_content(mysql_config):
    from planetgen.galaxy.sector import SpaceSector

    cfg = SystemConfig()
    cfg.BINARY_SYSTEM = False
    system = None
    for seed in range(30):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(system_config=cfg, galactic_center_dist_ly=26000.0)
        if any(b.body_type == "t" for b in _bodies(system)):
            break
    sector = SpaceSector("Radiogenic", edge_ly=11.5)
    sector.add_system(system, position=(0.0, 0.0, 0.0), system_config=cfg)
    store.save_sector(sector, config=mysql_config, galaxy_position={
        "center_x_pc": 100.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 100.0})
    conn = store.get_connection(mysql_config)
    try:
        rows = conn.execute("SELECT c_k0_ppm, c_th0_ppm, c_u0_ppm, al26_ratio, body_type FROM planets "
                            "WHERE body_type = 't' ORDER BY id").fetchall()
        stored = sorted(rows, key=lambda r: r["c_th0_ppm"])
        made = sorted((b for b in _bodies(system) if b.body_type == "t" and not b.is_moon), key=lambda b: b.c_th0_ppm)
        assert len(stored) == len(made) > 0
        for row, body in zip(stored, made):
            assert row["c_th0_ppm"] == pytest.approx(body.c_th0_ppm)
    finally:
        conn.close()


def test_a_rocky_body_heats_with_its_own_rock_share_and_aluminium():
    """GEN.428: the rock share and the Al2O3 of the body's drawn composition replace the Earth's stand-ins."""
    from types import SimpleNamespace
    body = SimpleNamespace(metal_mass_fraction=0.3, rock_mass_fraction=0.7, water_mass_fraction=0.0, ice_mass_fraction=0.0,
                           envelope_mass_fraction=0.0, rock_feo_fraction=0.08, rock_al2o3_fraction=0.09, rock_mg_si_ratio=1.2,
                           envelope_helium_fraction=None, envelope_z_fraction=None, porosity=None, carbon_rock=False,
                           composition_origin="drawn")
    rock = radiogenic.rock_of(body)
    assert rock.silicate == pytest.approx(0.7) and rock.metal == pytest.approx(0.3) and rock.wet == 0.0
    assert rock.aluminium == pytest.approx(2 * radiogenic.BULK_SILICATE_EARTH_AL_MASS_FRACTION)
    body.metal_mass_fraction = None  # saved before GEN.428: no composition, the Earth's stand-ins
    stand_in = radiogenic.rock_of(body)
    assert stand_in == radiogenic.Rock(tuning.RADIOGENIC_SILICATE_MASS_FRACTION, 1.0 - tuning.RADIOGENIC_SILICATE_MASS_FRACTION, 0.0,
                                       radiogenic.BULK_SILICATE_EARTH_AL_MASS_FRACTION)


def test_the_aluminium_26_of_a_star_is_a_spiky_sum_of_decayed_additions_about_its_regions_mean():
    """GEN.476: the median star is born with far less than the region's mean and a few with many times more, and the mean follows the gas depletion time and the star's aluminium."""
    mean = radiogenic.al26_mean()
    stream = draw.Stream("al26-spiky")
    ratios = sorted(radiogenic.draw_content(0.0, rng=stream).al26_ratio for _ in range(1500))
    assert ratios[len(ratios) // 2] < 0.3 * mean and ratios[-1] > 5.0 * mean
    assert 0.3 * mean < sum(ratios) / len(ratios) < 3.0 * mean
    assert radiogenic.al26_mean(fe_h=0.3) == pytest.approx(radiogenic.al26_mean() / 10.0 ** 0.3)
    assert radiogenic.al26_mean(fe_h=None) == radiogenic.al26_mean()


def test_the_given_mean_sets_the_scale_and_one_seed_gives_one_ratio():
    low = [radiogenic.draw_content(0.0, rng=draw.Stream(f"m{i}"), al26_mean_ratio=1e-7).al26_ratio for i in range(200)]
    high = [radiogenic.draw_content(0.0, rng=draw.Stream(f"m{i}"), al26_mean_ratio=1e-5).al26_ratio for i in range(200)]
    assert sum(high) > 20 * sum(low)
    assert radiogenic.draw_content(0.0, rng=draw.Stream("same")).al26_ratio == radiogenic.draw_content(0.0, rng=draw.Stream("same")).al26_ratio


def test_every_body_of_a_star_shares_its_birth_aluminium_26_and_a_star_without_a_history_has_none_to_share():
    from types import SimpleNamespace
    star = SimpleNamespace(metallicity_fe_h=-0.1, age=3.0, mass=1.0, galactic_center_dist_ly=26000.0)
    first, second = radiogenic.star_al26_ratio(star), radiogenic.star_al26_ratio(SimpleNamespace(**vars(star)))
    assert first == second and first >= 0.0
    other = radiogenic.star_al26_ratio(SimpleNamespace(metallicity_fe_h=-0.1, age=3.0, mass=1.1, galactic_center_dist_ly=26000.0))
    assert other != first or first == 0.0
    assert radiogenic.star_al26_ratio(SimpleNamespace(metallicity_fe_h=None, age=3.0)) is None
    content = radiogenic.draw_content(0.0, rng=draw.Stream("shared"), al26_ratio=first)
    assert content.al26_ratio == first
    assert draw.Stream("shared").random() == draw.Stream("shared").random()
