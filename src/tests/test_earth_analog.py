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


def _calls_in_order(source, names):
    """The `module.function` calls of `source` (a function's text) whose last two parts are in `names`, in the order they appear."""
    import ast

    alias = {"planetPhysics": "planets", "planetLife": "life"}
    found = []
    for node in sorted((n for n in ast.walk(ast.parse(source)) if isinstance(n, ast.Call)), key=lambda n: (n.lineno, n.col_offset)):
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
    finish = ["interior_heat.update_heat", "planets.refresh_at_final_age", "planets.apply_class", "magnetism.generate_stellar_induction", "life.apply_life_data"]
    air = ["planets.calculate_atmospheric_conditions", "planets.prepare_class"]
    assert _calls_in_order(inspect.getsource(Planet.__init__), set(air)) == air
    assert _calls_in_order(inspect.getsource(planets.generate_spin), set(tail)) == tail
    assert _calls_in_order(inspect.getsource(StarSystem._finish_planet), set(finish)) == finish
    # the script's `_generate` is the three in turn (the class after the life is `apply_life_data`'s own, GEN.980)
    script = _calls_in_order(inspect.getsource(analog._generate), set(air + tail + finish))
    assert script == air + tail + finish


def test_each_trial_starts_on_a_fresh_body(analog):
    """TEST.739: a trial leaves nothing on the body the next one draws (the light the life step reads was left over from the trial before)."""
    import inspect

    assert "copy.deepcopy(planet0)" in inspect.getsource(analog.trials)
