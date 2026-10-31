"""PERF.182: real runs of the generator stopped partway and carried on with `planetgen finish-run` end with the galaxy an
uninterrupted run makes (the done condition of the item): the same sectors, the same links, the same paths, the same
population. Each test stops a run in a different place, finishes it in one database and compares it with the same run
done whole in a second. A second `finish-run` after it has nothing to do."""

import hashlib

import pytest

from planetgen.db import fingerprint, store
from planetgen.generation import finish_run, followed, run_galaxy, run_population
from planetgen.population import model
from planetgen.queue import redisqueue, runrecord, units
from planetgen.db import sector_paths

from tests import worker_patches
from tests.bughunt_support import mysql_argv, run_cli
from tests.test_galaxy_gen import _seed_skeleton
from tests.test_gen_resume import MODES, NUM_SYSTEMS, SHELLS, second_mysql_config  # noqa: F401 -- a fixture


@pytest.fixture
def conn(redis_server):
    return redisqueue.connect(redis_server)


def _digest(config, table, leave_out=("id",)):
    """The rows of `table` as one hash, row ids and the columns in `leave_out` left out, in an order that does not depend on
    the order they were written in."""
    conn = store.get_connection(config)
    try:
        rows = conn.execute(f"SELECT * FROM {table}").fetchall()
    finally:
        conn.close()
    lines = sorted(repr(sorted((name, value) for name, value in dict(row).items() if name not in leave_out))
                   for row in rows)
    return len(lines), hashlib.sha256("\n".join(lines).encode()).hexdigest()


GALAXY = ("sectors", "star_systems", "stars", "planets", "bright_stars")
DERIVED = {"nearest_systems": ("id",), "sector_paths": ("id", "computed_at"),
           "sector_path_knots": ("id",), "species": ("id",), "polities": ("id",), "system_owners": ("id",)}


def _galaxy(config):
    """What a galaxy holds that a test compares: the fingerprint of its content (sectors, systems, stars, bright stars) and a
    digest of the tables rebuilt from it (links, paths, population)."""
    conn = store.get_connection(config)
    try:
        content = fingerprint.region_fingerprint(conn).lines()
    finally:
        conn.close()
    return content, {table: _digest(config, table, leave_out) for table, leave_out in DERIVED.items()}


def _same_galaxy(config, other):
    """Asserts the two databases hold the same galaxy (`_galaxy`), naming the first table that differs."""
    content, derived = _galaxy(config)
    other_content, other_derived = _galaxy(other)
    assert content == other_content
    for table in DERIVED:
        assert derived[table] == other_derived[table], f"{table} differs"


def _last_run(conn):
    """The run that started last: the ids sort by the second they started in, and two runs of a test can share one."""
    return max(runrecord.recent_runs(conn, limit=50), key=lambda run_id: runrecord.read_run(conn, run_id)["created"])


def _stops(monkeypatch, tmp_path, owner, name, call_number):
    worker_patches.patch_everywhere(monkeypatch, owner, name, "tests.worker_patches:fail_on_call",
                                    counter=str(tmp_path / f"{name}.calls"), call_number=call_number)


def test_a_sector_run_stopped_before_a_save_is_carried_on_to_the_uninterrupted_galaxy(
        mysql_config, second_mysql_config, tmp_path, conn):
    _seed_skeleton(mysql_config, shells=SHELLS)
    with pytest.MonkeyPatch.context() as failing:
        _stops(failing, tmp_path, run_galaxy, "generate_and_save_sector_at", 4)
        with pytest.raises(units.RunStopped, match="Interrupted"):
            run_cli("sectors", [*MODES["band"], "--num-systems", str(NUM_SYSTEMS), *mysql_argv(mysql_config)])
    stopped = _last_run(conn)
    assert runrecord.read_run(conn, stopped)["state"] == "failed"
    saved = _digest(mysql_config, "sectors")[0]
    assert 0 < saved

    run_cli("finish-run", mysql_argv(mysql_config))

    finished = runrecord.read_run(conn, _last_run(conn))
    assert finished["state"] == "done" and finished["finishes"] == stopped
    _seed_skeleton(second_mysql_config, shells=SHELLS)
    run_cli("sectors", [*MODES["band"], "--num-systems", str(NUM_SYSTEMS), *mysql_argv(second_mysql_config)])
    assert _digest(mysql_config, "sectors")[0] == _digest(second_mysql_config, "sectors")[0] > saved
    _same_galaxy(mysql_config, second_mysql_config)
    with pytest.raises(SystemExit):         # nothing is left to carry on: the newest run is done
        run_cli("finish-run", mysql_argv(mysql_config))


def _tail(args, started_at):
    """The stages after the sectors that need no scatter: link, the sector paths, the population and the territories."""
    return [run_galaxy.LinkFollower(args, started_at), run_galaxy.SettleFollower(args, started_at),
            run_population.population_follower(args), run_population.territories_follower(args)]


def _sectors_then_tail(config, stops=None, tmp_path=None):
    """Fills a shell of `config` and runs the tail on it; `stops` is `(owner, name)`: the call that stops the run."""
    from planetgen.generation import run_common
    _seed_skeleton(config, shells=SHELLS)
    conn = store.get_connection(config)
    try:
        started_at = store.database_now(conn)
    finally:
        conn.close()
    run_cli("sectors", [*MODES["band"], "--num-systems", str(NUM_SYSTEMS), *mysql_argv(config)])
    args = run_galaxy._default_generation_args(config=config)
    args.workers = 1
    if stops is None:
        followed.run_followed(args, _tail(args, started_at), "Tail")
        return
    with pytest.MonkeyPatch.context() as failing:
        _stops(failing, tmp_path, *stops, 1)
        with pytest.raises(units.RunStopped, match="Interrupted"):
            followed.run_followed(args, _tail(args, started_at), "Tail")
    run_common._RUN_STATS.clear()


@pytest.mark.parametrize("stage, stops", [
    ("Sector paths", (sector_paths, "settle_sectors")),
    ("Population pass", (model, "run_population")),
])
def test_a_run_of_the_stages_after_the_sectors_stopped_in_a_stage_is_carried_on_to_the_uninterrupted_galaxy(
        mysql_config, second_mysql_config, tmp_path, conn, stage, stops):
    _sectors_then_tail(mysql_config, stops, tmp_path)
    stopped = runrecord.read_run(conn, _last_run(conn))
    assert stopped["state"] == "failed" and stopped["title"] == "Tail"
    where = stopped["stages"][stopped["current"]]["label"]
    assert where == stage and stopped["galaxy_seed"] and stopped["build"]
    assert finish_run.last_incomplete(conn) == stopped["id"]

    run_cli("finish-run", mysql_argv(mysql_config))

    finished = runrecord.read_run(conn, _last_run(conn))
    assert finished["state"] == "done" and finished["finishes"] == stopped["id"]
    assert [stage_["units"] for stage_ in finished["stages"][:stopped["current"]]] == [0] * stopped["current"], \
        "the stages the stopped run finished were carried over, not done again"
    _sectors_then_tail(second_mysql_config)
    _same_galaxy(mysql_config, second_mysql_config)


def test_a_job_stopped_before_its_later_scatter_opened_is_carried_on_to_the_uninterrupted_galaxy(
        mysql_config, second_mysql_config, tmp_path, conn):
    """PERF.203 step 5: the job is one run, so the stages it never started are the stopped run's own. Finishing it plans them from
    the run record, and a scatter whose stage had not opened takes its mass limit when it opens, not from the stars the stopped
    stage left behind."""
    from planetgen.generation import job, run_common, run_plan
    from tests.test_bright_star_scatter import _seed_galaxy
    from tests.test_run_job import FLOOR, _job_args, _stars

    spec = [("scatter", ["--mass", *FLOOR]), ("scatter", ["--luminosity", *FLOOR])]
    _seed_galaxy(mysql_config)
    with pytest.MonkeyPatch.context() as failing:
        _stops(failing, tmp_path, run_plan, "scatter_shells_unit", 1)
        with pytest.raises(units.RunStopped, match="Interrupted"):
            job.run_job(_job_args(mysql_config, spec))
    run_common._RUN_STATS.clear()
    stopped = runrecord.read_run(conn, _last_run(conn))
    assert stopped["state"] == "failed" and stopped["title"] == "Generate job"
    assert stopped["stages"][stopped["current"]]["label"] != stopped["stages"][-1]["label"]

    run_cli("finish-run", mysql_argv(mysql_config))

    finished = runrecord.read_run(conn, _last_run(conn))
    assert finished["state"] == "done" and finished["finishes"] == stopped["id"]
    _seed_galaxy(second_mysql_config)
    job.run_job(_job_args(second_mysql_config, spec))
    assert _stars(mysql_config) == _stars(second_mysql_config) != []
