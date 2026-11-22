"""OPS.292: a run can name its start. The `sectors` entry of a job that gives `--shell S --band B --wedge W --radius-pc R` starts at
that sector and charts the sectors within R of it, as the random start it names did, so two full runs of one seed can be
compared digest to digest. A random start says the arguments that name it again, and a pinned run stopped partway and carried on
by `finish-run` ends with the galaxy of the same pinned run done whole."""

import pytest

from planetgen.db import fingerprint, store
from planetgen.generation import run_galaxy
from planetgen.queue import redisqueue, runrecord, units

from tests.bughunt_support import mysql_argv, run_cli
from tests.test_finish_run_generation import _last_run, _stops
from tests.test_galaxy_gen import _plan_wide_galaxy
from tests.test_gen_resume import second_mysql_config  # noqa: F401 -- a fixture

RADIUS_PC = 4.0
SEED = bytes.fromhex("0123456789ABCDEF0123456789ABCDEF")
SAID = "To start another run at this same sector, give its sectors entry: "


@pytest.fixture
def conn(redis_server):
    return redisqueue.connect(redis_server)


def _content(config):
    """The fingerprint of the galaxy a database holds, the sectors, systems, stars and bright stars in it."""
    conn = store.get_connection(config)
    try:
        return fingerprint.region_fingerprint(conn).lines()
    finally:
        conn.close()


def test_the_arguments_that_name_a_start_say_its_address_and_radius():
    assert run_galaxy.pinned_start_arguments((319, 59, 1516), 12.0) == "--shell 319 --band 59 --wedge 1516 --radius-pc 12"
    assert run_galaxy.pinned_start_arguments((2, -1, 0), 4.5) == "--shell 2 --band -1 --wedge 0 --radius-pc 4.5"


def test_a_random_start_says_the_arguments_that_name_the_same_start_again(capsys, monkeypatch):
    monkeypatch.setattr(run_galaxy.run_sector, "_log_summary", lambda summary: None)
    saved = {"name": "Alxipromae", "sector_id": 1, "summary": None}
    run_galaxy._start_announcer(12.0)(saved, (319, 59, 1516), "")
    said = [line for line in capsys.readouterr().out.splitlines() if SAID in line]
    assert said and said[0].split(SAID)[1] == "--shell 319 --band 59 --wedge 1516 --radius-pc 12"


def test_a_pinned_run_stopped_and_carried_on_gives_the_galaxy_of_the_whole_run(
        mysql_config, second_mysql_config, tmp_path, conn):
    pinned = ["--shell", "900", "--band", "0", "--wedge", "0", "--radius-pc", str(RADIUS_PC)]
    _plan_wide_galaxy(second_mysql_config, galaxy_seed=SEED)
    run_cli("sectors", [*pinned, *mysql_argv(second_mysql_config)])
    whole = _content(second_mysql_config)
    assert len(whole) > 3

    _plan_wide_galaxy(mysql_config, galaxy_seed=SEED)
    with pytest.MonkeyPatch.context() as failing:
        _stops(failing, tmp_path, run_galaxy, "generate_and_save_sector_at", 3)
        with pytest.raises(units.RunStopped, match="Interrupted"):
            run_cli("sectors", [*pinned, *mysql_argv(mysql_config)])
    stopped = _last_run(conn)
    assert runrecord.read_run(conn, stopped)["state"] == "held"
    assert _content(mysql_config) != whole

    run_cli("finish-run", mysql_argv(mysql_config))
    finished = runrecord.read_run(conn, _last_run(conn))
    assert finished["state"] == "done" and finished["finishes"] == stopped
    assert _content(mysql_config) == whole
