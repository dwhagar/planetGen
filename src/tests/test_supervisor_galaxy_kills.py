"""OPS.213, the half that needs a galaxy: real processes, a real kill in the conclusion of a stage, and the galaxy the disturbed run makes against the undisturbed one.

These run a small galaxy through a real supervisor and real workers twice: once untouched and once with the process that is ending a stage (the one that holds the lease `runrecord.ending_key`, the
unit that completed the stage) killed with SIGKILL in the middle of it. The run must end `done` and the galaxy must be the one the undisturbed run made (OPS.186's resume and the supervisor's
reconciler, OPS.210, are what make it so).

They need a machine that can finish a galaxy run (the worker of a small box is killed by its memory limit), so they are skipped unless `PLANETGEN_ENGINE_KILL_TESTS=1`:

    cd src && PLANETGEN_ENGINE_KILL_TESTS=1 python -m pytest tests/test_supervisor_galaxy_kills.py -q -p no:cacheprovider -n 1

Not run by their author (the box could not finish a galaxy): the first run on the gate machine is their first run.
"""

import os
import signal
import subprocess
import sys
import threading
import time

import pytest

from planetgen.queue import redisqueue, runrecord, supervisor, wrapper

pytestmark = [pytest.mark.real_queue, pytest.mark.fresh_database, pytest.mark.slow,
              pytest.mark.skipif(os.environ.get("PLANETGEN_ENGINE_KILL_TESTS") != "1", reason="needs a machine that can finish a galaxy run (PLANETGEN_ENGINE_KILL_TESTS=1)")]

SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
KILL_WITHIN_SECONDS = 1500


@pytest.fixture
def host(redis_server, tmp_path, monkeypatch):
    """A real supervisor for the test's runs (the suite's own `ask_supervisor` stand-in answers one turn of the starts only and would start no worker), a state folder of its own, and the
    processes it starts, which end with the test. As in `test_supervisor_kill.py` the processes know only the real key prefix."""
    monkeypatch.setenv("PLANETGEN_STATE_DIR", str(tmp_path))
    monkeypatch.setattr(runrecord, "PREFIX", "planetgen")
    env = {**os.environ, "PLANETGEN_STATE_DIR": str(tmp_path), "PYTHONPATH": os.pathsep.join([SRC, *sys.path])}
    processes = []
    conn = redisqueue.connect(redis_server)
    for key in conn.scan_iter(match="planetgen:*"):
        conn.delete(key)

    def start(argv):
        processes.append(subprocess.Popen(argv, env=env, cwd=SRC))

    monkeypatch.setattr(wrapper, "ask_supervisor", lambda connection, url: supervisor.ensure_started(connection, url or redis_server, start=start))
    yield conn
    pid = conn.hget(supervisor.alive_key(), "pid")
    if pid:
        try:
            os.kill(int(runrecord._s(pid)), signal.SIGKILL)
        except ProcessLookupError:
            pass
    for process in processes:
        if process.poll() is None:
            process.kill()
        process.wait()
    for key in conn.scan_iter(match="planetgen:*"):
        conn.delete(key)


def _digest(config):
    """What the galaxy is, in rows that do not depend on the order the sectors were made in: the sectors charted, the paths between them and the systems."""
    from planetgen.db import store
    db = store.get_connection(config)
    try:
        sectors = [tuple(row.values()) for row in db.execute(
            "SELECT s.shell_index, s.band_index, s.wedge_index FROM sectors s JOIN chart_ledger c ON c.shell_index = s.shell_index"
            " AND c.band_index = s.band_index AND c.wedge_index = s.wedge_index ORDER BY 1, 2, 3").fetchall()]
        paths = db.execute("SELECT COUNT(*) AS n FROM sector_paths").fetchone()["n"]
        systems = db.execute("SELECT COUNT(*) AS n FROM star_systems").fetchone()["n"]
    finally:
        db.close()
    return sectors, paths, systems


def _galaxy(config):
    from tests.bughunt_support import mysql_argv, run_cli
    from tests.test_galaxy_gen import _seed_skeleton
    _seed_skeleton(config, shells=[(0, 0), (1, 1), (2, 1)])
    run_cli("galaxy", ["--shell", "2", "--num-systems", "2", *mysql_argv(config)])
    return _digest(config)


def _killer(conn, killed, stop):
    """Kills, once, the process that holds the lease of a stage's end: the unit that completed the stage, in the middle of its conclusion."""
    pattern = runrecord.ending_key("*", "*")
    deadline = time.monotonic() + KILL_WITHIN_SECONDS
    while not stop.is_set() and time.monotonic() < deadline and not killed:
        for key in conn.scan_iter(match=pattern):
            pid = conn.get(key)
            if pid and int(runrecord._s(pid)) != os.getpid():
                try:
                    os.kill(int(runrecord._s(pid)), signal.SIGKILL)
                except ProcessLookupError:
                    continue
                killed.append((runrecord._s(key), int(runrecord._s(pid))))
                return
        time.sleep(0.02)


def test_a_galaxy_whose_stage_end_was_killed_is_the_galaxy_an_undisturbed_run_makes(mysql_config, host):
    from tests.db_schema_support import scratch_database
    undisturbed = _galaxy(mysql_config)
    assert undisturbed[0] and undisturbed[1] > 0, "the undisturbed run made a galaxy"
    killed, stop = [], threading.Event()
    thread = threading.Thread(target=_killer, args=(host, killed, stop), daemon=True)
    thread.start()
    try:
        with scratch_database(mysql_config) as other:
            disturbed = _galaxy(other)
    finally:
        stop.set()
        thread.join(timeout=10)
    assert killed, "a process was found in the middle of ending a stage and killed"
    runs = [runrecord.read_run(host, run_id) for run_id in runrecord.recent_runs(host)]
    assert runs and all(record["state"] == "done" for record in runs), [(record["state"], record.get("error")) for record in runs]
    assert disturbed == undisturbed, "the galaxy of the disturbed run is the undisturbed run's"
