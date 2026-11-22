# tests/test_earth_analog.py

"""TEST.726: the Earth analog as a standing script (`scripts/earth_analog.py`): the comparison text, the never-blocks rule, and that a trial redraws the baseline's row."""

import importlib.util
import json
import pathlib
import sys

import pytest

SCRIPTS = pathlib.Path(__file__).resolve().parents[2] / "scripts"


def _load(name):
    sys.path.insert(0, str(SCRIPTS))
    try:
        spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(SCRIPTS))
    return module


@pytest.fixture(scope="module")
def analog():
    return _load("earth_analog")


def _row(i, cls, **extra):
    row = {"i": i, "w": 0.001, "cls": cls, "hyd": "surface ocean", "T": 288.0, "P": 101325.0, "oc": 0.7, "life_class": 1, "fails": []}
    row.update(extra)
    return row


def test_the_summary_counts_classes_bins_medians_and_failing_groups(analog):
    rows = [_row(0, "M"), _row(1, "F", T=400.0, fails=["O2", "T_mean"]), _row(2, "F", T=200.0, fails=["O2"]), {"i": 3, "w": 0.1, "error": "ValueError: x"}]
    summary = analog.summarize(rows)
    assert summary["trials"] == 4 and summary["errors"] == 1 and summary["classes"] == {"M": 1, "F": 2}
    assert summary["temperature"] == {"258-343 K": 1, ">343 K": 1, "<258 K": 1} and summary["median_temperature_k"] == 288
    assert summary["m_fails"]["O2"] == 2 and summary["m_fails"]["T_mean"] == 1 and summary["m_fails"]["gravity"] == 0


def test_the_comparison_is_paired_and_names_the_class_moves(analog):
    baseline = {"made_on": "abc1234", "block": 1, "age_gy": 4.6, "water": "sample", "rows": [_row(0, "F"), _row(1, "F"), _row(2, "O"), _row(3, "F")]}
    same = analog.compare_text(baseline, [_row(0, "F"), _row(1, "F")])                           # the first two trials only
    assert "Trials that changed class: none." in same and "2 trials" in same
    moved = analog.compare_text(baseline, [_row(0, "M"), _row(1, "F"), _row(2, "O"), _row(3, "F")])
    assert "F to M 1" in moved and "| class M | 0 | 1 (+1) |" in moved


def test_a_missing_baseline_or_a_failure_is_a_sentence_never_an_error(analog, monkeypatch, tmp_path):
    monkeypatch.setattr(analog, "BASELINE_PATH", str(tmp_path / "none.json"))
    assert analog.report(10).startswith("Earth analog: not run (")


def test_pr_body_prints_the_section_and_can_leave_it_out(monkeypatch, capsys):
    train = _load("train")
    analog = sys.modules["earth_analog"] if "earth_analog" in sys.modules else _load("earth_analog")
    monkeypatch.setattr(analog, "report", lambda n: f"REPORT {n}")
    train.pr_body(type("Args", (), {"earth_analog": 40, "no_earth_analog": False})())
    out = capsys.readouterr().out
    assert "## Known issues" in out and "## Earth analog" in out and "REPORT 40" in out
    train.pr_body(type("Args", (), {"earth_analog": 40, "no_earth_analog": True})())
    assert "## Earth analog" not in capsys.readouterr().out


def test_the_stored_baseline_is_well_formed_and_its_first_trials_redraw_the_same(analog):
    baseline = analog.load_baseline()
    assert baseline["rows"] and baseline["rows"][0]["i"] == 0 and len(baseline["rows"]) == 800
    again = analog.trials(baseline["block"], 0, 2, baseline["age_gy"], baseline["water"])
    assert [json.dumps(r, sort_keys=True) for r in again] == [json.dumps(r, sort_keys=True) for r in baseline["rows"][:2]]


def test_the_stored_earth_water_rows_are_the_baselines_and_redraw_the_same(analog):
    """TEST.820: the Earth-water row (the Earth's own sea, not the drawn water) is stored beside the sampled rows and its first trials redraw the same."""
    baseline = analog.load_baseline()
    assert baseline["earth_water"] == analog.EARTH_WATER and len(baseline["earth_rows"]) == analog.EARTH_WATER_TRIALS
    assert all(abs(row["w"] - analog.EARTH_WATER) < 1e-12 for row in baseline["earth_rows"])
    again = analog.trials(baseline["block"], 0, 2, baseline["age_gy"], analog.EARTH_WATER)
    assert [json.dumps(r, sort_keys=True) for r in again] == [json.dumps(r, sort_keys=True) for r in baseline["earth_rows"][:2]]


def test_the_earth_water_text_sets_the_rows_against_the_earths_numbers(analog):
    rows = [_row(0, "M", T=290.0, P=118000.0, oc=0.55), _row(1, "F", T=292.0, P=116000.0, oc=0.55), {"i": 2, "w": analog.EARTH_WATER, "error": "ValueError: x"}]
    summary = analog.summarize_earth(rows)
    assert summary["trials"] == 3 and summary["errors"] == 1 and summary["class_m"] == 1
    assert summary["median_temperature_k"] == 291.0 and summary["median_pressure_kpa"] == 117.0 and summary["median_ocean_cover"] == 0.55
    assert summary["temperature_gap_k"] == 3.0 and summary["pressure_gap_kpa"] == 15.7
    text = analog.earth_water_text({"earth_rows": [_row(0, "M", T=288.0, P=101000.0)]}, rows)
    assert "288 K and 101 kPa" in text and "| median temperature K | 288.0 | 291.0 | 288 |" in text
    assert "| median pressure kPa | 101.0 | 117.0 | 101 |" in text
    assert "Paired change over the 1 trials both have: temperature +2.0 K, pressure +17.0 kPa" in text
    assert "none" in analog.earth_water_text({}, rows)                    # a baseline from before the row: the stored column says none, nothing raises


def test_the_analog_sets_the_veneer_the_air_reads_and_pins_what_the_generated_planet_would_bring(analog, monkeypatch):
    """TEST.820: the line that set `planet._late_accreted_percent` was read by nothing (GEN.1829), so the planet kept its system's own veneer (3.0 percent, then 0.28) and its mantle law and radiogenic content. The veneer is on the build record the air reads,
    the interior fields are fixed, and the veneer moves the air (so the pin is the one the generator reads)."""
    from planetgen import tuning
    from planetgen.generation.build_record import record_of
    from planetgen.physics import impact_history, radiogenic

    planet, _class_m = analog._setup(analog.AGE_GY)
    assert record_of(planet).late_accreted_percent == impact_history.LATE_ACCRETED_PERCENT_EARTH
    assert not hasattr(planet, "_late_accreted_percent")
    assert planet.mantle_eta0_pa_s == sum(tuning.INTERIOR_VISCOSITY_PREFACTOR_RANGE_PA_S) / 2.0
    assert planet.mantle_activation_kj_mol == sum(tuning.INTERIOR_ACTIVATION_RANGE_KJ_MOL) / 2.0
    assert planet.c_k0_ppm == radiogenic.BULK_SILICATE_EARTH[0] * radiogenic.formation_mix("K")[1] and planet.al26_ratio == radiogenic.SOLAR_AL26_RATIO
    earth = analog.trials(1, 0, 1, analog.AGE_GY, analog.EARTH_WATER)[0]
    monkeypatch.setattr(impact_history, "LATE_ACCRETED_PERCENT_EARTH", 10.0 * impact_history.LATE_ACCRETED_PERCENT_EARTH)
    assert analog.trials(1, 0, 1, analog.AGE_GY, analog.EARTH_WATER)[0]["P"] != earth["P"]


def _calls_in_order(source, names):
    """The `module.function` calls of `source` (a function's text) whose last two parts are in `names`, in the order they appear."""
    import ast

    alias = {"planetPhysics": "planets", "planetLife": "life"}
    import textwrap

    found = []
    for node in sorted((n for n in ast.walk(ast.parse(textwrap.dedent(source))) if isinstance(n, ast.Call)), key=lambda n: (n.lineno, n.col_offset)):
        func = node.func
        if isinstance(func, ast.Attribute) and isinstance(func.value, ast.Name):
            name = f"{alias.get(func.value.id, func.value.id)}.{func.attr}"
            if name in names:
                found.append(name)
    return found


def test_the_script_runs_the_generators_steps_in_the_generators_order(analog):
    """TEST.739: the trial calls the system's own steps in the order `generate_spin`'s tail and `StarSystem._finish_planet` run them, so it cannot drift from the generator again."""
    import inspect

    from planetgen.generation.planet import Planet
    from planetgen.generation.system import StarSystem
    from planetgen.physics import planets

    tail = ["daylight.generate", "magnetism.generate_field", "radiation.generate_dose", "body_colour.generate", "habitability_world.generate"]
    finish = ["interior_heat.update_heat", "planets.refresh_at_final_age", "planets.prepare_class", "magnetism.generate_stellar_induction", "life.apply_life_data"]
    air = ["planets.calculate_atmospheric_conditions", "planets.prepare_class"]
    assert _calls_in_order(inspect.getsource(Planet.__init__), set(air)) == air
    assert _calls_in_order(inspect.getsource(planets.generate_spin), set(tail)) == tail
    assert _calls_in_order(inspect.getsource(StarSystem._finish_planet), set(finish)) == finish
    # the script's `_generate` is the three in turn (the class is read once, after the life, by `apply_life_data`'s own `apply_class`, GEN.1745)
    script = _calls_in_order(inspect.getsource(analog._generate), set(air + tail + finish))
    assert script == air + tail + finish


def test_each_trial_starts_on_a_fresh_body(analog):
    """TEST.739: a trial leaves nothing on the body the next one draws (the light the life step reads was left over from the trial before)."""
    import inspect

    assert "copy.deepcopy(planet0)" in inspect.getsource(analog.trials)
