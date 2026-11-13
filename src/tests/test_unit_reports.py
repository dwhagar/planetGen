"""OPS.170: the units of `generation/run_plan.py` that used to take `report` and never call it say their phases, report their progress and, in a stage of one unit, move a bar."""

import ast
import inspect

import pytest

from planetgen.generation import run_plan
from planetgen.util import log


def _database(mysql_config):
    return {"host": mysql_config.host, "port": mysql_config.port, "user": mysql_config.user, "database": mysql_config.database}


def test_no_unit_of_the_plan_module_takes_report_and_never_uses_it():
    """The review's M2 (the standard's rule 6): a unit that is given `report` calls it, or the bar stands still for the whole unit."""
    tree = ast.parse(inspect.getsource(run_plan))
    silent = []
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name.endswith("_unit") and any(arg.arg == "report" for arg in node.args.args):
            uses = [name for name in ast.walk(node) if isinstance(name, ast.Name) and name.id == "report" and isinstance(name.ctx, ast.Load)]
            if not uses:
                silent.append(node.name)
    assert silent == [], f"units that never report: {silent}"


def test_the_phases_of_a_unit_are_said_reported_and_stamped(monkeypatch):
    lines, reports, beats = [], [], []
    monkeypatch.setattr(log, "normal", lines.append)
    monkeypatch.setattr(run_plan.unitbeat, "phase", beats.append)
    phases = run_plan._UnitPhases(lambda done, total: reports.append((done, total)), 3)
    phases.next("Emptying the thing")
    phases.next("Making the other thing...")
    phases.finish()
    assert reports == [(0, 3), (1, 3), (3, 3)]
    assert beats == ["Emptying the thing", "Making the other thing..."]
    assert len(lines) == 2 and lines[0].startswith("Emptying the thing (") and lines[1].startswith("Making the other thing (") and lines[0].endswith(" elapsed).")


def test_the_phases_of_one_of_a_stages_many_units_are_debug_records(monkeypatch):
    normal, debug = [], []
    monkeypatch.setattr(log, "normal", normal.append)
    monkeypatch.setattr(log, "debug", debug.append)
    run_plan._UnitPhases(lambda done, total: None, 1, loud=False).next("One range of ids")
    assert normal == [] and len(debug) == 1


def test_a_walk_over_things_reports_every_so_many_and_at_the_end(monkeypatch):
    reports = []
    step = run_plan._counting(lambda done, total: reports.append((done, total)), "Things, systems", 250, every=100)
    for done in range(1, 251):
        step(done)
    assert reports == [(100, 250), (200, 250), (250, 250)]


class _Bar:
    def __init__(self, *args, **kwargs):
        self.args, self.kwargs, self.updates = args, kwargs, []

    def update(self, **kwargs):
        self.updates.append(kwargs)

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        return False


def test_a_stage_of_one_unit_moves_a_bar_from_the_units_reports_and_fills_it_when_the_unit_is_back(monkeypatch):
    bars = []
    monkeypatch.setattr(run_plan.steps, "Step", lambda *a, **k: bars.append(_Bar(*a, **k)) or bars[-1])
    follower = run_plan.ClearScatterFollower.__new__(run_plan.ClearScatterFollower)
    run_plan.followed.Follower.__init__(follower, None)
    follower.key, follower.label = "mass-clear", "Clear the massive stars"
    follower.progress = object()
    follower.on_plan({"units": 1})
    assert bars[0].args[0] == "Clear the massive stars" and bars[0].kwargs["percent"] is True and bars[0].kwargs["record"] is False
    follower.on_progress(0, 1, 3)
    follower.on_progress(0, 3, 3)
    follower.on_result(0, {"lines": []})
    assert [update["completed"] for update in bars[0].updates] == [1 / 3, 1.0, 1.0]
    follower.end(None)


@pytest.mark.db
def test_the_clearing_units_say_their_phase_and_report_the_steps(mysql_config, monkeypatch):
    lines = []
    monkeypatch.setattr(log, "normal", lines.append)
    params = {"database": _database(mysql_config)}
    for name, extra in (("empty_phenomenon_scatter_unit", {}), ("drop_phenomenon_indexes_unit", {}), ("clear_bright_stars_unit", {"clear": "all"})):
        reports = []
        getattr(run_plan, name)(dict(params, **extra), lambda done, total: reports.append((done, total)))
        assert reports and reports[-1] == (1, 1) and reports[0][0] == 0, name
    assert len(lines) == 3 and all(" elapsed)." in line for line in lines)


@pytest.mark.db
def test_the_fill_and_the_stamped_system_units_report_their_ranges(mysql_config):
    from tests import rescan_proof_support as proof
    proof.small_galaxy(mysql_config)
    with proof.connection(mysql_config) as conn:
        low, high = conn.execute("SELECT MIN(id) AS low, MAX(id) AS high FROM star_systems").fetchone().values()
    params = {"low": low, "high": high, "all": True, "database": _database(mysql_config)}
    for unit in (run_plan.rescore_habitability_unit, run_plan.refill_effect_ledger_unit):
        reports = []
        result = unit(params, lambda done, total: reports.append((done, total)))
        assert reports[0][0] == 0 and reports[-1] == (result["systems"], result["systems"]) and result["systems"] > 0, unit.__name__
    reports = []
    run_plan.fill_map_luminosity_unit({"low": 1, "high": 10, "database": _database(mysql_config)}, lambda done, total: reports.append((done, total)))
    assert reports[0] == (0, 3) and reports[-1] == (3, 3)

