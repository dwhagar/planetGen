"""
OPS.196 and OPS.262: no command goes silent. Every registered command is cut into its units, as the starter that asks for it cuts it (`expand_function`, with the arguments the starter hands it), and
its first units are run on a small galaxy under the beat a worker gives a unit (`unitbeat.UnitBeat`); the longest stretch the unit's count stood still (`longest_gap`, OPS.135: the pieces done, the
second-level count, a phase or a statement's label moving it) is held to a budget. A command with no row in `ARGUMENTS` and none in `NOT_RUN` fails the guard at the foot of this file, so a new command
cannot go unmeasured.
A command that has put its work back into one loop of its own, or one piece that runs for minutes with no count, shows as a gap here before it shows as a stuck bar on the Queue page.
The registry guard is `test_commands_have_pieces.py`; the budget here is for a galaxy of a few sectors, so a gap over it is a command that waits for something that is not its work.
"""


import os

import pytest

from planetgen.cli import generate
from planetgen.generation import run_galaxy, run_plan
from planetgen.queue import commands, redisqueue, unitbeat, wrapper
from planetgen.web import warmup
from tests.test_api import _place_sector
from tests.test_galaxy_gen import EDGE_PC, _seed_skeleton

pytestmark = pytest.mark.fresh_database

BUDGET_SECONDS = 5.0
"""float: The longest a command may stand without a count, a phase or a statement label moving, on a small galaxy."""

RECORDED_SILENCE = {
    "generate-start-sector": (12.0, "OPS.265: on some runs one body's seasonal solve (`seasonal_swing.temperature_series`, in `refresh_at_final_age` of a pre-placed system's moon) stands 4 to 6.6 s inside one call, though its mean is 66 ms; measured 6.2 to 6.6 s on 4 of 7 runs, none on the others"),
    "generate-neighborhood": (12.0, "the same single seasonal solve (6.5 s)"),
    "generate-sectors": (12.0, "the same single seasonal solve (6.3 s)"),
}
"""dict: The commands that go silent longer than `BUDGET_SECONDS` today, `{command: (their budget, what was measured)}` (OPS.262). Each is a bug to fix (a count inside the piece), not a norm: the budget is about twice what was measured, so a command that gets worse still fails, and a fixed one has its row deleted."""

UNITS_RUN = 3
"""int: The most units of one command that are run (the first ones: the cut of a large stage repeats the same pieces)."""

FLOOR = 1.0e5
"""float: The luminosity floor, in solar luminosities, the scatter cases are cut at."""

SEED = 424242
"""int: The seed of the scatter cases."""


def _parsed(*argv, config=None, run_seed=7):
    """The arguments the starter `planetgen *argv` parses, with the database and the run's seed set the way the starter sets them."""
    args = generate.build_parser()[0].parse_args(list(argv))
    if args.command in ("sector", "galaxy"):
        generate._shape_like_sector(args)
    if config is not None:
        args.database = wrapper.database_args(config)
        args.run_seed = run_seed
    return args


def _step_params(config, *argv, key=None):
    """The `params` of the step of the run `planetgen *argv` builds (the one named `key`, else the first), the way the wrapper builds them."""
    args = _parsed(*argv, config=config)
    steps = commands.build_steps(args.command, args)
    return dict(next(step for step in steps if key is None or step.key == key).params)


def _plain(extra):
    return lambda config: dict(extra)


def _sector_ids(config):
    from planetgen.db import store
    conn = store.get_connection(config)
    try:
        return [row["id"] for row in conn.execute("SELECT id FROM sectors ORDER BY id").fetchall()]
    finally:
        conn.close()


def _galaxy_sectors(config, items):
    """The cut of the sectors stages: the run's arguments as JSON, the sector edge and the addresses (shell, band, wedge) of the sectors to make."""
    return {"args": run_galaxy._sector_args_json(_parsed("galaxy", "--num-systems", "1", config=config)), "edge_pc": EDGE_PC, "items": list(items)}


def _site(config):
    """The warm commands name the site's database, not the connection's arguments."""
    return {"database": config.database, "cap": None}


def _scatter(config, **more):
    return {"min_luminosity_sol": FLOOR, "mass_limit_sol": 8.0, "seed": SEED, **more}


ARGUMENTS = {
    "check-counts": [("", _plain({}))],
    "check-datasets": [("", _plain({}))],
    "check-ids": [("", _plain({}))],
    "check-math": [("", _plain({}))],
    "check-object-index": [("", _plain({}))],
    "check-orphans": [("", _plain({}))],
    "check-revision": [("", _plain({}))],
    "check-storage": [("", _plain({}))],
    "check-systems": [("", _plain({}))],
    "check-tables": [("", _plain({}))],
    "check-tables-extended": [("", _plain({}))],
    "check-values": [("", _plain({}))],
    "check-version-keys": [("", _plain({}))],
    "clean-nebula-field": [("", _plain({"apply": False}))],
    "clear-bright-stars": [("", _plain({"clear": "all"}))],
    "clear-scatter": [("phenomena", _plain({"part": "phenomena"})), ("nebula-field", _plain({"part": "nebula-field"}))],
    "draw-nebula-field": [("", _plain({}))],
    "drop-phenomenon-indexes": [("", _plain({}))],
    "empty-nebula-field": [("", _plain({}))],
    "empty-phenomenon-scatter": [("", _plain({}))],
    "end-map-rebuild-maintenance": [("", _plain({}))],
    "fill-map-luminosity": [("", _plain({}))],
    "refill-effect-ledger": [("", _plain({}))],
    "refresh-block-lights": [("", _plain({}))],
    "refresh-run-blocks": [("", _plain({}))],
    "rescore-habitability": [("", _plain({}))],
    "reset-galaxy": [("", _plain({}))],
    "run-population": [("", _plain({}))],
    "run-territories": [("", _plain({}))],
    "scatter-mass": [("", lambda config: _scatter(config))],
    "scatter-luminosity": [("", lambda config: _scatter(config))],
    "scatter-band": [("", lambda config: _scatter(config, current_luminosity_sol=FLOOR, target_luminosity_sol=FLOOR / 10))],
    "scatter-phenomena": [("", lambda config: {"min_mass_solar": 8.0, "seed": SEED})],
    "scatter-backfill": [("", lambda config: {"since": None})],
    "link-neighbors": [("", lambda config: {"since": None, "sector_ids": _sector_ids(config)})],
    "settle-sector-paths": [("", lambda config: {"since": None, "sector_ids": _sector_ids(config), "workers": 1})],
    "store-scatter-limits": [("", lambda config: {"prepare": {"min_luminosity_sol": FLOOR, "mass_limit_sol": 8.0, "seed": SEED}})],
    "plan-skeleton": [("", lambda config: run_plan.skeleton_params(_parsed("plan", "--seed", "0123456789abcdef0123456789abcdef", config=config), config))],
    "generate-unplaced-sector": [("", lambda config: _step_params(config, "sector", "--num-sectors", "1", "--num-systems", "1"))],
    "generate-system": [("", lambda config: _step_params(config, "system", "--name", "Silent Star"))],
    "generate-start-sector": [("", lambda config: _galaxy_sectors(config, [0, 0, 0]))],
    "generate-neighborhood": [("", lambda config: _galaxy_sectors(config, [1, 0, 0, 1, 0, 1]))],
    "generate-sectors": [("", lambda config: _galaxy_sectors(config, [1, 0, 0, 1, 0, 1]))],
    "write-special-phenomena": [("", lambda config: {"seed": SEED, "min_mass_solar": 8.0})],
    "build-phenomenon-indexes": [("", lambda config: {"seed": SEED, "min_mass_solar": 8.0})],
    "record-phenomenon-scatter": [("", lambda config: {"seed": SEED, "min_mass_solar": 8.0})],
    "warm-cluster-field": [("", _site)],
    "warm-object-index": [("", _site)],
    "warm-galaxy-map": [("", _site)],
    "warm-galaxy-nodes": [("", _site)],
    "warm-galaxy-deep-nodes": [("", _site)],
    "warm-lists": [("", _site)],
    "warm-sectors": [("", _site)],
    "warm-systems": [("", _site)],
    "warm-phenomena": [("", _site)],
    "warm-search": [("", _site)],
}
"""dict: Each registered command that is run, `{command: [(case id, builder)]}`: the builder is given the throwaway database's config and returns the arguments the starter hands the command's cut (beside the
database, which is added when the builder does not give one)."""

NOT_RUN = {
}
"""dict: The commands this test does not run, `{command: why}`: a reason that stands in the way, never "later"."""


CASES = [(kind, case, build) for kind, cases in ARGUMENTS.items() for case, build in cases]


@pytest.mark.parametrize("kind, case, build", CASES, ids=[f"{kind}{'-' + case if case else ''}" for kind, case, build in CASES])
def test_a_command_never_stands_without_a_count_past_the_budget(mysql_config, redis_server, monkeypatch, kind, case, build):
    real_app = warmup.build_app
    monkeypatch.setattr(warmup, "build_app", lambda config=None: real_app(mysql_config))      # the warm commands visit the pages of this database
    _seed_skeleton(mysql_config, shells=[(0, 0), (1, 1)])
    _place_sector(mysql_config, "One", (5.0, 5.0, 5.0))
    if kind.startswith("warm-"):
        warmup.keep_app()      # a worker builds the warm-up app before its unit starts (PERF.257), so the imports are not the unit's silence
    params = build(mysql_config)
    params.setdefault("database", wrapper.database_args(mysql_config))
    units = [unit for unit, _cost in commands.expand_function(kind)(params).units][:UNITS_RUN]
    gaps = []
    if os.environ.get("SILENCE_DEBUG"):
        close = unitbeat.UnitBeat._close_gap

        def noisy(self, now):
            if now - self._changed_at - self._waited > 2.0:
                print("GAP", kind, round(now - self._changed_at - self._waited, 1), self._phase, getattr(self, "_step", None), flush=True)
            close(self, now)
        monkeypatch.setattr(unitbeat.UnitBeat, "_close_gap", noisy)
    for unit in units:
        with unitbeat.UnitBeat(redisqueue.connect(redis_server), "silence", 0, kind) as report:
            commands.run_alone(kind, unit, report)
            gaps.append(report.longest_gap())
    budget = RECORDED_SILENCE.get(kind, (BUDGET_SECONDS, ""))[0]
    assert max(gaps, default=0.0) < budget, f"{kind} stood {max(gaps):.1f} s with no count, phase or statement moving (budget {budget} s)"


def test_every_registered_command_is_run_here_or_has_its_reason_not_to_be():
    unlisted = sorted(set(commands.REGISTRY) - set(ARGUMENTS) - set(NOT_RUN))
    assert not unlisted, f"commands the silence run does not know (give each a row in ARGUMENTS, the arguments its starter hands it, or a reason in NOT_RUN): {unlisted}"
    assert not set(ARGUMENTS) & set(NOT_RUN), "a command is either run or not run"
    assert not sorted((set(ARGUMENTS) | set(NOT_RUN)) - set(commands.REGISTRY)), "a row for a command that is not registered"
    assert all(len(reason.split()) >= 4 for reason in NOT_RUN.values()), "a NOT_RUN reason says what stands in the way"


def test_a_recorded_silence_is_of_a_registered_command_and_is_longer_than_the_budget():
    assert set(RECORDED_SILENCE) <= set(ARGUMENTS), "a recorded silence of a command this file does not run"
    assert all(budget > BUDGET_SECONDS and reason for budget, reason in RECORDED_SILENCE.values()), "a row that fits the ordinary budget is deleted"


def test_the_silence_meter_sees_a_unit_that_does_not_move(redis_server):
    """The meter is not a rubber stamp: a count that stands still for 30 s of a clock the test moves is a gap of 30 s."""
    now = [0.0]
    beat = unitbeat.UnitBeat(redisqueue.connect(redis_server), "silence", 0, "meter", clock=lambda: now[0])
    beat(1, 2)
    now[0] = 30.0
    beat(2, 2)
    assert beat.longest_gap() == 30.0
