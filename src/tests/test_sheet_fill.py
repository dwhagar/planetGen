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
    assert sheet.applies("planet", "M", "temp_day_side_k", locked) and not sheet.applies("planet", "M", "temp_polar_k", locked)
    assert sheet.applies("planet", "M", "temp_polar_k", turning) and not sheet.applies("planet", "M", "temp_day_side_k", turning)
    # a rule that reads the body does not rule a column out when the body is not known
    assert sheet.applies("planet", "M", "ocean_depth_km") and sheet.applies("planet", "M", "temp_polar_k")


def test_the_new_rules_read_the_physics_they_state():
    dry_air = {"atmospheric_pressure_pa": 1.0e5}
    assert sheet.applies("planet", "M", "atm_density", dry_air) and not sheet.applies("moon", "A", "atm_density", {"atmospheric_pressure_pa": 0.0})
    assert not sheet.applies("planet", "J", "lid_regime") and sheet.applies("planet", "M", "lid_regime")
    locked = {"rotation_period_hours": 8766.0, "period_years": 1.0}
    assert not sheet.applies("planet", "M", "solar_day_hours", locked) and sheet.applies("planet", "M", "solar_day_hours", {"rotation_period_hours": 24.0, "period_years": 1.0})
    assert sheet.applies("planet", "C", "shelter_column_g_cm2", {"eq_shelter": 0, "eq_shelter_rec": 2})
    assert not sheet.applies("planet", "C", "shelter_column_g_cm2", {"eq_shelter": 0, "eq_shelter_rec": 0})
    assert sheet.applies("moon", "C", "outside_hours", {"eq_shelter": 2}) and not sheet.applies("moon", "C", "outside_hours", {"eq_shelter": 0})
    assert sheet.applies("moon", "C", "ice_shell_km", {"hydrosphere": "ice"}) and not sheet.applies("moon", "C", "ice_shell_km", {"hydrosphere": "dry"})
    assert sheet.applies("moon", "C", "hp_ice_km", {"hydrosphere": "ice-covered ocean"}) and not sheet.applies("moon", "C", "hp_ice_km", {"hydrosphere": "surface ocean"})
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


def test_star_and_comet_rules():
    assert sheet.applies("star", "M", "wide_binary_a_crit_km", {"role": "primary"}) and not sheet.applies("star", "M", "wide_binary_a_crit_km", {"role": "single"})
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
