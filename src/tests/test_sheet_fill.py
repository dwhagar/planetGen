"""
UX.208: the fill audit of the property sheet. The registry `describe.sheet.APPLIES` says where a column does not apply;
a small fixed-seed sample is generated and tallied; every column the registry rules out is checked against the
physical rule it states, and the tally has the shape the report prints.
"""

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
