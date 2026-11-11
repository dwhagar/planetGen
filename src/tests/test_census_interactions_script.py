# tests/test_census_interactions_script.py

"""GEN.1448: the interaction census script (`scripts/census_interactions.py`) counts the bodies a change moved, apart from the ones it did not."""

import importlib.util
import pathlib

SPEC = importlib.util.spec_from_file_location("census_interactions", pathlib.Path(__file__).resolve().parents[2] / "scripts" / "census_interactions.py")
census = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(census)


def _body(**values):
    base = {"ocean_tide_spring_m": 1.0, "tidal_flux_w_m2": 0.01, "lock_state": "free", "rotation_period_hours": 24.0, "axial_tilt_deg": 23.0, "eccentricity": 0.05}
    base.update(values)
    return base


def _system(index, *planets, star="G2V"):
    return {"star": star, "index": index, "planets": [{"body": body, "moons": moons} for body, moons in planets]}


def test_the_same_draw_twice_moves_nothing():
    run = [_system(0, (_body(), [_body()])), _system(1, (_body(), []))]
    row = census.compare(run, run)["G2V"]
    assert row["matched"] == 2 and row["bodies"] == 3 and row["moved_any"] == 0


def test_a_tilt_change_beyond_the_tolerance_counts_for_that_quantity_only():
    before = [_system(0, (_body(), []))]
    after = [_system(0, (_body(axial_tilt_deg=40.0), []))]
    row = census.compare(before, after)["G2V"]
    assert row["moved"]["axial_tilt_deg"] == 1 and row["moved"]["eccentricity"] == 0 and row["moved_any"] == 1


def test_a_change_inside_the_tolerance_is_not_a_move():
    before = [_system(0, (_body(), []))]
    after = [_system(0, (_body(axial_tilt_deg=23.0 * 1.001), []))]
    assert census.compare(before, after)["G2V"]["moved_any"] == 0


def test_a_locking_state_change_is_a_move():
    before = [_system(0, (_body(), [_body(lock_state="free")]))]
    after = [_system(0, (_body(), [_body(lock_state="locked")]))]
    assert census.compare(before, after)["G2V"]["moved"]["lock_state"] == 1


def test_a_system_with_a_different_body_count_is_reshaped_not_matched():
    before = [_system(0, (_body(), [_body()]))]
    after = [_system(0, (_body(), []))]
    row = census.compare(before, after)["G2V"]
    assert row["reshaped"] == 1 and row["matched"] == 0 and row["bodies"] == 0


def test_an_abandoned_or_missing_system_is_failed():
    before = [_system(0, (_body(), [])), _system(1, (_body(), []))]
    after = [{"star": "G2V", "index": 0, "abandoned": True}]
    row = census.compare(before, after)["G2V"]
    assert row["failed"] == 2 and row["matched"] == 0


def test_kozai_losses_are_summed_from_the_bodies_that_record_them():
    before = [_system(0, (_body(), []))]
    after = [_system(0, (_body(kozai={"kozai_moons_lost": 2}), []))]
    row = census.compare(before, after)["G2V"]
    assert row["kozai_before"] == 0 and row["kozai_after"] == 2


def test_the_report_names_each_quantity():
    text = census.report(census.compare([_system(0, (_body(), []))], [_system(0, (_body(), []))]), "a.json", "b.json")
    assert all(label in text for _key, label in census.QUANTITIES) and "| all |" in text


def test_a_seed_is_the_same_each_time_and_differs_by_index():
    assert census.seed_for("G2V", 3) == census.seed_for("G2V", 3) != census.seed_for("G2V", 4)
