"""
UX.208: the fill audit of the property sheet. The registry `describe.sheet.APPLIES` says where a column does not apply;
a small fixed-seed sample is generated and tallied; every column the registry rules out is checked against the
physical rule it states, and the tally has the shape the report prints.
"""

import re

from planetgen.describe import sheet

from tests import sheet_fill_support as support


def test_a_giant_has_no_surface_and_a_rocky_world_has_no_envelope_core():
    assert sheet.applies("planet", "J", "interior_kind") and not sheet.applies("planet", "M", "interior_kind")
    assert sheet.applies("planet", "M", "phi4") and not sheet.applies("planet", "J", "phi4")
    assert sheet.applies("moon", "A", "atm_density") and not sheet.applies("moon", "C", "atm_density")
    assert sheet.applies("star", "G", "mass_kg") and sheet.applies("planet", "M", "mass_kg")


def test_conditional_rules_read_the_bodys_own_columns():
    ocean = {"hydrosphere": "surface ocean"}
    assert sheet.applies("planet", "M", "ocean_depth_km", ocean)
    assert not sheet.applies("planet", "M", "ocean_depth_km", {"hydrosphere": "dry"})
    assert sheet.applies("planet", "M", "life_stage", {"life_microbial": 2})
    assert not sheet.applies("planet", "M", "life_stage", {"life_microbial": 0, "life_flora": 0})
    assert not sheet.applies("planet", "M", "magnetopause_rp", {"magnetic_moment_a_m2": 0})
    locked = {"rotation_period_hours": 8766.0, "period_years": 1.0}
    turning = {"rotation_period_hours": 24.0, "period_years": 1.0}
    assert not sheet.applies("planet", "M", "day_share_equator", locked) and sheet.applies("planet", "M", "day_share_equator", turning)
    # a rule that reads the body does not rule a column out when the body is not known
    assert sheet.applies("planet", "M", "ocean_depth_km") and sheet.applies("planet", "M", "day_share_equator")


def test_the_new_rules_read_the_physics_they_state():
    dry_air = {"atmospheric_pressure_pa": 1.0e5}
    assert sheet.applies("planet", "M", "atm_density", dry_air) and not sheet.applies("moon", "A", "atm_density", {"atmospheric_pressure_pa": 0.0})
    assert not sheet.applies("planet", "J", "lid_regime") and sheet.applies("planet", "M", "lid_regime")
    locked = {"rotation_period_hours": 8766.0, "period_years": 1.0}
    assert not sheet.applies("planet", "M", "solar_day_hours", locked) and sheet.applies("planet", "M", "solar_day_hours", {"rotation_period_hours": 24.0, "period_years": 1.0})
    assert sheet.applies("planet", "C", "shelter_column_g_cm2", {"eq_shelter": 0, "eq_shelter_rec": 2})
    assert not sheet.applies("planet", "C", "shelter_column_g_cm2", {"eq_shelter": 0, "eq_shelter_rec": 0})
    assert sheet.applies("moon", "C", "outside_hours", {"outside_hours": 800.0}) and not sheet.applies("moon", "C", "outside_hours", {"outside_hours": None})
    assert sheet.applies("moon", "C", "ice_shell_km", {"hydrosphere": "ice"}) and not sheet.applies("moon", "C", "ice_shell_km", {"hydrosphere": "dry"})
    assert sheet.applies("moon", "C", "hp_ice_km", {"hydrosphere": "ice-covered ocean", "water_base_pressure_pa": 1.0e9})
    assert not sheet.applies("moon", "C", "hp_ice_km", {"hydrosphere": "surface ocean", "water_base_pressure_pa": 1.0e9})
    assert not sheet.applies("moon", "C", "hp_ice_km", {"hydrosphere": "ice-covered ocean", "water_base_pressure_pa": 2.0e8})    # GEN.781: too shallow for ice VI
    assert sheet.applies("planet", "M", "ring_kind", {"ring_kind": "ice"}) and not sheet.applies("planet", "M", "ring_mass_kg", {"ring_kind": None})
    assert sheet.applies("moon", "C", "aquifer_depth_km", {"aquifer_depth_km": 1.0}) and not sheet.applies("moon", "C", "aquifer_depth_km", {})
    assert sheet.applies("moon", "C", "parent_shadow_share") and not sheet.applies("planet", "C", "parent_shadow_share")
    assert sheet.applies("moon", "C", "induced_moment_a_m2") and sheet.applies("planet", "C", "induced_moment_a_m2")      # GEN.981: a planet in its star's field


def test_the_ages_after_oxygen_apply_only_to_a_world_that_reached_them():
    """A world whose life waits for oxygen that never comes stores no age for the stages from complex cells on (TEST.449)."""
    for key in ("life_age_complex_cells_gyr", "life_age_multicellularity_gyr", "life_age_technological_civilization_gyr"):
        assert not sheet.applies("planet", "M", key, {"life_microbial": 2, "life_stage": "photosynthesis"})
        assert not sheet.applies("moon", "C", key, {"life_microbial": 2, "life_stage": None})
        assert sheet.applies("planet", "M", key, {"life_microbial": 2, "life_stage": "complex_cells"})
        assert sheet.applies("planet", "M", key)
    assert sheet.applies("planet", "M", "life_age_photosynthesis_gyr", {"life_microbial": 2, "life_stage": "photosynthesis"})


def test_the_tide_heat_applies_to_a_rocky_or_icy_body_on_an_orbit_with_an_eccentricity():
    """GEN.778: `interior_heat.update_heat` stores the flow for a rocky or icy body whose orbit is eccentric; a giant or a circular orbit stores none."""
    key = "tidal_flux_w_m2"
    assert sheet.applies("moon", "D", key, {"eccentricity": 0.01}) and sheet.applies("planet", "C", key, {"eccentricity": 0.2})
    assert not sheet.applies("moon", "D", key, {"eccentricity": 0.0}) and not sheet.applies("moon", "D", key, {"eccentricity": 3e-170})
    assert not sheet.applies("planet", "I", key, {"eccentricity": 0.05}) and not sheet.applies("planet", "J", key, {"eccentricity": 0.05})
    assert ("planet", key) not in sheet.OPEN_GAPS and ("moon", key) not in sheet.OPEN_GAPS


def test_a_generated_body_stores_its_tide_wherever_the_sheet_says_it_applies():
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw
    from tests.fuzz_support import deterministic_entropy
    seen = {True: 0, False: 0}
    for seed in range(7100, 7112):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(SystemConfig(), galactic_center_dist_ly=26000.0)
        for planet in system.planets:
            if not hasattr(planet, "moons"):
                continue
            for body in [planet, *planet.moons]:
                if body.body_type == "g":
                    assert body.tidal_flux_w_m2 is None
                    continue
                eccentric = (body.eccentricity or 0.0) >= sheet.NEGLIGIBLE_TIDE_ECCENTRICITY
                assert body.tidal_flux_w_m2 is not None or not eccentric, (body.name, body.eccentricity)      # empty only where the rule excuses it
                seen[eccentric] += 1
    assert seen[True] > 5 and seen[False] > 0


def test_porosity_applies_to_a_body_lighter_than_the_taper_end():
    """GEN.780: `bulk_draw` stores a porosity for a heap lighter than 1e21 kg and none above it."""
    assert sheet.applies("moon", "D", "porosity", {"mass_kg": 5e20}) and not sheet.applies("moon", "C", "porosity", {"mass_kg": 3e22})
    assert not sheet.applies("planet", "N", "porosity", {"mass_kg": 6e24})
    assert ("planet", "porosity") not in sheet.OPEN_GAPS and ("moon", "porosity") not in sheet.OPEN_GAPS
    assert ("planet", "hp_ice_km") not in sheet.OPEN_GAPS and ("moon", "hp_ice_km") not in sheet.OPEN_GAPS


def test_a_generated_body_stores_porosity_and_high_pressure_ice_wherever_the_sheet_says_they_apply():
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw
    from tests.fuzz_support import deterministic_entropy
    seen = {"light": 0, "heavy": 0, "deep": 0}
    for seed in range(7100, 7106):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(SystemConfig(), galactic_center_dist_ly=26000.0)
        for planet in system.planets:
            if not hasattr(planet, "moons"):
                continue
            for body in [planet, *planet.moons]:
                if body.body_type == "g":
                    continue
                light = body.mass < sheet.POROSITY_MAX_MASS_KG
                assert (body.porosity is not None) == light, (body.name, body.mass, body.porosity)
                seen["light" if light else "heavy"] += 1
                deep = (body.hydrosphere in sheet.ICE_SHELL_STATES and (body.water_base_pressure_pa or 0.0) >= sheet.HP_ICE_MIN_BASE_PA)
                if deep:
                    assert body.hp_ice_km, (body.name, body.water_base_pressure_pa)
                    seen["deep"] += 1
                else:
                    assert body.hydrosphere not in sheet.ICE_SHELL_STATES or not body.hp_ice_km or (body.water_base_pressure_pa or 0.0) >= sheet.HP_ICE_MIN_BASE_PA
    assert seen["light"] > 0 and seen["heavy"] > 5


def test_a_giant_keeps_no_starlight_while_its_moons_and_a_rocky_planet_do():
    """GEN.784: `insolation` and `habitability_world` leave a giant's `energy_flux_w_m2` and `flux_variation` none; every other planet and moon stores them."""
    for key in ("energy_flux_w_m2", "flux_variation"):
        assert not sheet.applies("planet", "J", key, {}) and not sheet.applies("planet", "I", key, {})
        assert sheet.applies("planet", "N", key, {}) and sheet.applies("moon", "D", key, {})
        assert ("planet", key) not in sheet.OPEN_GAPS and ("moon", key) not in sheet.OPEN_GAPS


def test_a_placed_star_stores_its_rate_of_lethal_galactic_events(mysql_config):
    """GEN.785: the rate is worked out from the star's place in the galaxy (`activity.set_galactic_hazard`) and stored with the star; a system
    built with no place is rated where such a system is taken to be (`scales.default_galactocentric_ly`), so every star has one."""
    from planetgen.db import store
    from planetgen.generation.config import SystemConfig
    from planetgen.galaxy import scales
    from planetgen.generation.system import StarSystem
    from planetgen.physics import activity
    from planetgen.util import draw
    assert ("star", "lethal_event_rate_per_gyr") not in sheet.OPEN_GAPS
    rates = {}
    for placed in (True, False):
        draw.set_run_seed(7100)
        config = SystemConfig()
        system = StarSystem(system_config=config, **({"galactic_center_dist_ly": support.SAMPLE_GALACTIC_DISTANCE_LY} if placed else {}))
        system_id = store.save_system(system, config, config=mysql_config)
        conn = store.get_connection(mysql_config)
        try:
            rates[placed] = [row["lethal_event_rate_per_gyr"] for row in conn.execute(
                "SELECT lethal_event_rate_per_gyr FROM stars WHERE star_system_id = ?", (system_id,)).fetchall()]
        finally:
            conn.close()
    assert rates[True] and all(rate is not None and rate > 0.0 for rate in rates[True])
    default_rate = activity.lethal_event_rate_per_gyr(scales.default_galactocentric_ly())
    assert rates[False] and all(rate == default_rate for rate in rates[False])


def test_star_and_comet_rules():
    key = "wide_binary_a_crit_km"
    assert sheet.applies("star", "M", key, {"role": "primary", "binary_configuration": "wide"})
    assert sheet.applies("star", "M", key, {"role": "secondary", "binary_configuration": "wide"})
    assert sheet.applies("star", "M", key, {"role": "member", "binary_configuration": "close"})      # GEN.789: a further member has its window
    assert not sheet.applies("star", "M", key, {"role": "primary", "binary_configuration": "close"})  # a close pair is orbited as one body
    assert not sheet.applies("star", "M", key, {"role": "single", "binary_configuration": None})
    assert ("star", key) not in sheet.OPEN_GAPS
    assert sheet.applies("comet", "comet", "orbital_period_years", {"orbit_type": "elliptical"})
    assert not sheet.applies("comet", "comet", "orbital_period_years", {"orbit_type": "parabolic"})
    assert sheet.applies("comet", "comet", "eccentricity", {"orbit_type": "parabolic"})


def test_every_open_gap_names_a_column_of_its_kind_and_an_item():
    for (kind, key), item in sheet.OPEN_GAPS.items():
        assert key in {spec.key for spec in sheet.SPECS[kind]}, f"{kind} has no column {key}"
        assert re.fullmatch(r"[A-Z]+\.\d+", item), item
        assert sheet.applies(kind, "M" if kind != "star" else "G", key), f"{kind} {key} is excused but ruled out by a rule"


def test_every_rule_names_columns_of_its_kind_and_gives_its_reason():
    for kind, rules in sheet.APPLIES.items():
        keys = {spec.key for spec in sheet.SPECS[kind]}
        for names, _test, reason in rules:
            assert reason and not reason.endswith(" ")
            assert set(names) <= keys, f"{kind}: {sorted(set(names) - keys)} are not columns of its sheet"


def test_the_tally_counts_every_body_and_splits_the_empty_columns(mysql_config):
    support.generate(mysql_config, 4, first_seed=7000)
    sampled = support.tally(mysql_config)
    kinds = {key[0] for key in sampled if key != "bare"}
    assert {"star", "planet"} <= kinds
    for key, entry in sampled.items():
        if key == "bare":
            continue
        assert entry["bodies"] > 0
        for filled, gap, not_applicable in entry["columns"].values():
            assert filled + gap + not_applicable == entry["bodies"]
    text = "\n".join(support.lines(sampled))
    assert "bodies" in text and "bare number" in text


def test_a_sample_has_no_empty_column_that_applies_and_is_not_filed(mysql_config):
    """UX.208: a column empty where it applies is a gap in generation; each is excused by the filed item that closes it."""
    support.generate(mysql_config, 12, first_seed=7100)
    missing = support.unexplained(support.tally(mysql_config))
    assert not missing, "empty where they apply, with no open item (file one, or add the rule that says why): " + "; ".join(
        f"{kind} {key} ({', '.join(f'{klass} {n}/{bodies}' for klass, bodies, n in rows)})" for (kind, key), rows in sorted(missing.items()))
