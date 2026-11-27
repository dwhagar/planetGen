# tests/test_interrupted_galaxy_run.py

"""
TEST.21, canceling a run: SIGTERM, Ctrl+C and `SystemExit` during a
parallel run end it as canceled and leave no half-written sector.
"""

import os
import signal
import subprocess
import sys
import time

import pytest

from planetgen.db import store as _db

from tests.test_galaxy_gen import _mysql_argv, _plan_wide_galaxy
from tests.test_work_queue import _run_states, control_config  # noqa: F401 -- fixture

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


START_WAIT_SECONDS = 480
"""How long a test waits for a run to save its first sector before it calls the run hung (TEST.141). It guards a hang,
not a speed: a gate that saturates the machine took longer than 120 s, and a run that dies fails at once whatever this is."""

GALAXY_SHELL = 5
# The test is about interrupting a run, not about generation size, so a sector holds two systems (more than one, so the
# signal can land between two systems of a sector): the first sector is saved in an eighth of the time (TEST.822 and
# TEST.851: with 16 systems at about 27 s each the first sector took 425 s on a quiet box, and with 4 systems one of three
# cases took 350 s and two passed 480 s while the box was shared, too close to START_WAIT_SECONDS). The run still cannot end on its own between the first sector being
# saved and the test's signal landing (TEST.101): shell 5 plans 296 sectors for 2 workers, so after the first sector
# the run has about 147 more rounds of sectors to go, hours of work, and the signal lands within seconds of the save.
NUM_SYSTEMS = 2


def _sector_counts(config):
    conn = _db.get_connection(config)
    try:
        sectors = conn.execute("SELECT COUNT(*) AS n FROM sectors").fetchone()["n"]
        short = conn.execute(
            "SELECT s.id, COUNT(ss.id) AS n FROM sectors s LEFT JOIN star_systems ss ON ss.sector_id = s.id"
            " JOIN chart_ledger c ON c.shell_index = s.shell_index AND c.band_index = s.band_index AND c.wedge_index = s.wedge_index"
            " GROUP BY s.id HAVING n < ?", (NUM_SYSTEMS,)).fetchall()      # a charted sector only: a rogue planet's landing cell has a row no fill charted (GEN.1111 part 7)   # a first fill adds bright stars' systems (GEN.633)
        orphans = conn.execute(
            "SELECT COUNT(*) AS n FROM star_systems ss LEFT JOIN sectors s ON s.id = ss.sector_id"
            " WHERE ss.sector_id IS NOT NULL AND s.id IS NULL").fetchone()["n"]
        conn.rollback()
        return sectors, short, orphans
    finally:
        conn.close()


_LAUNCH_GENERATE = ("from planetgen.generation import run_galaxy; run_galaxy.require_datasets = lambda: None; "
                    "from planetgen.cli import generate; generate.main()")
"""str: `python -m planetgen.cli.generate` with the data gate (GEN.590) passed: a test machine downloads no datasets."""


@pytest.mark.slow
@pytest.mark.parametrize("sig, group", [
    (signal.SIGINT, True),    # Ctrl+C in a terminal reaches the run and its workers
    (signal.SIGTERM, True),   # Cancel on the Generate page stops the step's whole process group
    (signal.SIGTERM, False),  # `kill <pid>` of the run alone
])
def test_interrupting_a_parallel_galaxy_run_leaves_no_half_written_sector(control_config, redis_server, sig, group):
    _plan_wide_galaxy(control_config)
    env = dict(os.environ, PLANETGEN_CONTROL_DATABASE=control_config.database, PLANETGEN_WORKERS="2")
    run = subprocess.Popen(
        [sys.executable, "-c", _LAUNCH_GENERATE, "sectors", "--shell", str(GALAXY_SHELL),
         "--num-systems", str(NUM_SYSTEMS), "--workers", "2"] + _mysql_argv(control_config),
        cwd=REPO, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, start_new_session=True,
        # A process inherits an ignored Ctrl+C, and the run then leaves it
        # alone ("someone else handles it") and finishes with status 0.
        # Another test in this pytest worker can leave SIGINT ignored
        # (TEST.101).
        preexec_fn=lambda: signal.signal(signal.SIGINT, signal.SIG_DFL))
    try:
        deadline = time.time() + START_WAIT_SECONDS
        while _sector_counts(control_config)[0] < 1:
            assert run.poll() is None, run.stderr.read().decode(errors="replace")[-2000:]
            assert time.time() < deadline, "no sector was saved"
            time.sleep(0.05)
        if group:
            os.killpg(run.pid, sig)
        else:
            os.kill(run.pid, sig)
        run.wait(timeout=300)
    finally:
        if run.poll() is None:
            os.killpg(run.pid, signal.SIGKILL)
            run.wait()
    assert run.returncode != 0

    sectors, short, orphans = _sector_counts(control_config)
    assert sectors >= 1 and short == [] and orphans == 0
    # The run's record on Redis (PERF.98): the run is another process, so it is under the real key prefix.
    assert _run_states(redis_server, control_config.database, prefix="planetgen") == ["canceled"]
