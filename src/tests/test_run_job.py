# tests/test_run_job.py

"""
PERF.203 step 3: `planetgen run-job` puts the commands of a Generate job in one run. The spec is checked as the commands'
own parsers check theirs, the stages are the ones the commands give (the block light totals once, at the last place), what a
command checks before it starts is checked when its own stage opens, and the stored output is what the commands give one after
another.
"""

import argparse
import json

import pytest

from planetgen.cli import generate as generate_cli
from planetgen.db import store
from planetgen.generation import followed, job, run_plan, stages
from planetgen.util import log
from tests.test_bright_star_scatter import THRESHOLD, _scatter_args, _seed_galaxy


@pytest.fixture(autouse=True)
def _the_console_log_writes_to_stdout():
    log.configure(log.NORMAL)


def _entries(mysql_config, spec):
    """`job.parse_entries` for `spec` (a list of `(command, [options])`) against `mysql_config`, as `process_args` does."""
    _parser, parsers = generate_cli.build_parser()
    shared = parsers["run-job"].parse_args(["x", "--mysql-host", mysql_config.host, "--mysql-port", str(mysql_config.port),
                                            "--mysql-user", mysql_config.user, "--mysql-password", mysql_config.password,
                                            "--mysql-database", mysql_config.database])
    text = json.dumps([{"command": command, "argv": argv} for command, argv in spec])
    return shared, job.parse_entries(text, shared, parsers, generate_cli.validate_command)



def _job_args(mysql_config, spec):
    shared, entries = _entries(mysql_config, spec)
    shared.command, shared.job = "run-job", entries
    return shared


FLOOR = ["--bright-star-min-luminosity", str(THRESHOLD)]
NEW_GALAXY = [("plan", []), ("sectors", ["--radius-pc", "4"]), ("link", []),
              ("scatter", ["--mass", *FLOOR]), ("scatter", ["--luminosity", *FLOOR]), ("scatter", ["--phenomena"]),
              ("backfill", []), ("settle", [])]


def test_the_spec_is_a_list_of_commands_and_their_options():
    with pytest.raises(ValueError, match="not JSON"):
        job.parse_spec("{")
    with pytest.raises(ValueError, match="list of"):
        job.parse_spec("{}")
    with pytest.raises(ValueError, match="list of"):
        job.parse_spec("[]")
    with pytest.raises(ValueError, match="entry 1 is not"):
        job.parse_spec('[{"command": "plan", "argv": "x"}]')
    with pytest.raises(ValueError, match="a job holds"):
        job.parse_spec('[{"command": "reset"}]')
    assert job.parse_spec('[{"command": "link"}, {"command": "scatter", "argv": ["--mass"]}]') == [
        ("link", []), ("scatter", ["--mass"])]


def test_an_entry_is_validated_by_its_own_command(mysql_config):
    with pytest.raises(ValueError, match="entry 1 .*not valid"):
        _entries(mysql_config, [("scatter", ["--mass", "--phenomenon-min-mass", "3"])])
    with pytest.raises(ValueError, match="entry 1 .*not valid"):
        _entries(mysql_config, [("scatter", [])])


def test_what_is_not_a_stage_of_a_run_is_refused(mysql_config):
    with pytest.raises(ValueError, match="band and the nebula field"):
        _entries(mysql_config, [("scatter", ["--band", "--down-to", "100"])])
    with pytest.raises(ValueError, match="band and the nebula field"):
        _entries(mysql_config, [("scatter", ["--nebula-field"])])
    with pytest.raises(ValueError, match="--sector-id"):
        _entries(mysql_config, [("settle", ["--sector-id", "4"])])


def test_a_job_has_the_stages_of_its_commands_with_the_block_lights_once_at_the_end_of_the_last_scatter(mysql_config):
    args = _job_args(mysql_config, NEW_GALAXY)
    parts = job.build_parts(args.job, args)           # nothing is read from a database that has no galaxy yet
    keys = [stage.key for stage in job.stage_list(parts)]
    assert keys == ["skeleton", "start", "neighborhood", "link", "mass", "luminosity", "phenomena", "phenomena-index",
                    "nebula-draw", "block-lights", "backfill", "settle"]
    followers = followed.dedupe([f for part in parts for f in part.followers])
    assert [f.key for f in followers] == keys
    for part in parts:
        assert [f.key for f in part.followers] == [stage.key for stage in part.stages]


def test_a_stage_that_is_not_repeatable_comes_once():
    class Stage(followed.Follower):
        key = "mass"

    with pytest.raises(ValueError, match="comes twice"):
        followed.dedupe([Stage(None), Stage(None)])


def test_the_settings_of_a_stage_are_the_ones_its_own_command_was_given(mysql_config):
    args = _job_args(mysql_config, [("scatter", ["--mass", "--phenomenon-min-mass", "8", *FLOOR]),
                                    ("scatter", ["--phenomena", "--compact-min-mass", "6"])])
    parts = job.build_parts(args.job, args)
    mass, phenomena = (part.args for part in parts)
    assert stages.stage_settings("mass", mass)["mass_limit_sol"] == 8.0
    assert stages.stage_settings("phenomena", phenomena)["compact_limit_sol"] == 6.0
    assert stages.stage_settings("phenomena", mass)["compact_limit_sol"] != 6.0


def test_the_clock_of_the_later_stages_is_the_one_the_sectors_stage_opens_at(mysql_config):
    args = _job_args(mysql_config, [("sectors", ["--radius-pc", "4"]), ("link", []), ("settle", [])])
    parts = job.build_parts(args.job, args)
    assert parts[1].followers[0].params()["since"] == {"deferred": "stage-clock", "stage": "start"}
    assert parts[2].followers[0].params()["since"] == {"deferred": "stage-clock", "stage": "start"}


def _stars(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        return [{key: value for key, value in dict(row).items() if key != "created_at"}
                for row in conn.execute("SELECT * FROM bright_stars ORDER BY id").fetchall()]
    finally:
        conn.close()


def _phenomena(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        return [{key: value for key, value in dict(row).items() if key != "created_at"}
                for row in conn.execute("SELECT * FROM phenomenon_scatter ORDER BY id").fetchall()]
    finally:
        conn.close()


@pytest.mark.db
def test_a_job_of_scatters_stores_what_the_verbs_store_one_after_another(mysql_config):
    _seed_galaxy(mysql_config)
    args = _job_args(mysql_config, [("scatter", ["--mass", *FLOOR]), ("scatter", ["--luminosity", *FLOOR]),
                                    ("scatter", ["--phenomena"])])
    job.run_job(args)
    by_job = (_stars(mysql_config), _phenomena(mysql_config))
    assert by_job[0] and by_job[1]

    conn = store.get_connection(mysql_config)
    try:
        store.clear_bright_stars(conn)
        store.clear_phenomenon_scatter(conn)
    finally:
        conn.close()
    run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--mass"), ("mass",), refuse="mass")
    run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--luminosity"), ("luminosity",), refuse="luminosity")
    from tests.test_phenomenon_scatter import _scatter_args as phenomena_args
    run_plan.scatter_phenomena(phenomena_args(mysql_config, "--phenomena"), refuse="phenomena")
    assert (_stars(mysql_config), _phenomena(mysql_config)) == by_job


@pytest.mark.db
def test_a_job_that_clears_a_scatter_draws_it_again_instead_of_refusing(mysql_config):
    _seed_galaxy(mysql_config)
    job.run_job(_job_args(mysql_config, [("scatter", ["--mass", *FLOOR]), ("scatter", ["--luminosity", *FLOOR])]))
    first = _stars(mysql_config)
    assert first
    job.run_job(_job_args(mysql_config, [("clear-scatter", ["--luminosity"]), ("scatter", ["--luminosity", *FLOOR])]))
    again = _stars(mysql_config)
    assert [star["id"] for star in again[:1]] == [star["id"] for star in first[:1]]      # the massive stars were left alone
    without_ids = lambda stars: [{key: value for key, value in star.items() if key != "id"} for star in stars]
    assert without_ids(again) == without_ids(first)         # the bright stars are drawn again, the same, under new ids


@pytest.mark.db
def test_a_job_scatter_of_a_kind_that_is_stored_refuses_in_words_when_its_stage_opens(mysql_config, capsys):
    _seed_galaxy(mysql_config)
    job.run_job(_job_args(mysql_config, [("scatter", ["--mass", *FLOOR])]))
    stored = _stars(mysql_config)
    with pytest.raises(SystemExit):
        job.run_job(_job_args(mysql_config, [("scatter", ["--mass", *FLOOR])]))
    assert "scatter is already stored -- run 'planetgen clear-scatter --mass' first" in capsys.readouterr().out
    assert _stars(mysql_config) == stored


@pytest.mark.db
def test_a_job_that_plans_the_galaxy_can_scatter_it_in_the_same_run(mysql_config):
    """Nothing is planned when the job's followers are built: the skeleton is what its first stage makes, and the scatters take
    their seed and mass limit from it when their stages open."""
    args = _job_args(mysql_config, [("plan", ["--seed", "0123456789ABCDEF0123456789ABCDEF", "--edge", "300pc"]), ("scatter", ["--mass", *FLOOR])])
    job.run_job(args)
    conn = store.get_connection(mysql_config)
    try:
        assert store.get_galaxy_shape(conn) is not None
        assert store.scatter_stored(conn, "mass")
    finally:
        conn.close()


@pytest.mark.db
def test_a_job_whose_sectors_have_no_plan_to_work_from_stops_in_words(mysql_config, capsys):
    with pytest.raises(SystemExit):
        job.run_job(_job_args(mysql_config, [("sectors", ["--radius-pc", "4"])]))
    assert "The galaxy has never been planned" in capsys.readouterr().out


@pytest.mark.db
def test_the_run_job_command_numbers_its_stages_across_the_job_and_records_one_run(mysql_config, monkeypatch, tmp_path):
    """What the Generate page's runner gives the one `run-job` step: the stages numbered from the step's offset, the page's progress
    file naming the stage the run is in, one entry in the run history, and the stored galaxy the commands make."""
    import sys

    from planetgen.queue import progress_file
    from planetgen.util import logrecord  # noqa: F401 -- the stage marks are log records

    progress = tmp_path / "progress.json"
    monkeypatch.setenv(progress_file.ENV_VAR, str(progress))
    monkeypatch.setenv(stages.STAGE_OFFSET_ENV, "2")
    monkeypatch.setenv(stages.STAGE_TOTAL_ENV, "6")
    spec = json.dumps([{"command": "plan", "argv": ["--seed", "0123456789ABCDEF0123456789ABCDEF", "--edge", "300pc"]},
                       {"command": "scatter", "argv": ["--mass", *FLOOR]},
                       {"command": "scatter", "argv": ["--luminosity", *FLOOR]}])
    monkeypatch.setattr(sys, "argv", [
        "planetgen", "run-job", spec, "--disregard-missing-data", "--mysql-host", mysql_config.host, "--mysql-port",
        str(mysql_config.port), "--mysql-user", mysql_config.user, "--mysql-password", mysql_config.password,
        "--mysql-database", mysql_config.database])
    generate_cli.main()

    reported = json.loads(progress.read_text())["stage"]
    assert (reported["index"], reported["total"], reported["label"]) == (4, 4, "Make the block light totals")
    conn = store.get_connection(mysql_config)
    try:
        assert store.get_galaxy_shape(conn) is not None
        assert store.scatter_stored(conn, "mass") and store.scatter_stored(conn, "luminosity")
        runs = conn.execute("SELECT command, outcome FROM generation_runs").fetchall()
    finally:
        conn.close()
    assert [(row["command"], row["outcome"]) for row in runs] == [("run-job", "ok")]
