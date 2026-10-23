"""PERF.119 (F4): `planetgen sector` and `planetgen system` are runs of one command through the wrapper: the entry point builds
the step list from what was asked (`wrapper.steps_for`) and hands it to `wrapper.run`, a worker generates and saves, and the
run is on the run record like any other."""

import json

import pytest

from planetgen.cli import generate
from planetgen.generation import run_sector
from planetgen.queue import commands, redisqueue, runrecord, wrapper

SECTOR = ["sector", "--num-sectors", "3", "--num-systems", "1"]
SYSTEM = ["system", "--name", "Wrapped Star"]


def _args(argv, run_seed=7):
    args = generate.build_parser()[0].parse_args(argv)
    args.database = {"host": "h", "port": 3306, "user": "u", "database": "d"}
    args.run_seed = run_seed
    return args


def _argv(config, *argv):
    return ["generate", *argv[:1], "--mysql-host", config.host, "--mysql-port", str(config.port),
            "--mysql-user", config.user, "--mysql-password", config.password, "--mysql-database", config.database,
            *argv[1:]]


@pytest.mark.parametrize("argv, kind", [(SECTOR, "generate-unplaced-sector"), (SYSTEM, "generate-system")])
def test_the_request_gives_a_list_of_the_one_registered_command(argv, kind):
    args = _args(argv)
    steps = wrapper.steps_for(args.command, args)
    assert [step.kind for step in steps] == [kind]
    assert kind in commands.REGISTRY and commands.REGISTRY[kind].unit and commands.REGISTRY[kind].expand
    assert steps == wrapper.steps_for(args.command, _args(argv))
    assert [commands.Step.from_dict(json.loads(json.dumps(step.as_dict()))) for step in steps] == steps


def test_a_sector_run_is_labelled_with_its_count_and_a_system_run_has_one_unit():
    assert wrapper.steps_for("sector", _args(SECTOR))[0].label == "Sectors (3 unplaced)"
    assert wrapper.steps_for("sector", _args(["sector"]))[0].params["count"] == 1
    [step] = wrapper.steps_for("system", _args(SYSTEM))
    expansion = commands.expand_function(step.kind)(step.params)
    assert len(expansion.units) == 1 and expansion.units[0][1] == 1.0


def test_the_sectors_are_cut_into_chunks_of_consecutive_numbers_the_same_at_any_worker_count():
    for count in (1, 5, 256, 257, 1000):
        step = wrapper.steps_for("sector", _args(["sector", "--num-sectors", str(count)]))[0]
        units = commands.expand_function(step.kind)(step.params).units
        numbers = [n for unit, _cost in units for n in unit["numbers"]]
        assert numbers == list(range(count)) and len(units) <= run_sector.UNPLACED_UNITS_MAX
        assert sum(cost for _unit, cost in units) == count


def test_each_sector_draws_from_a_stream_of_its_own_that_the_run_seed_and_its_number_fix():
    seeds = [run_sector._task_seed(5, n) for n in range(50)]
    assert len(set(seeds)) == 50
    assert seeds == [run_sector._task_seed(5, n) for n in range(50)]
    assert run_sector._task_seed(6, 0) != seeds[0]


def test_no_other_command_is_a_one_command_run():
    with pytest.raises(ValueError):
        commands.build_steps("phenomenon", _args(["galaxy"]))


def test_a_sector_run_starts_on_the_run_record_and_saves_its_sectors(mysql_config, monkeypatch, redis_server):
    seen = []
    real = wrapper.run

    def spy(title, steps, *args, **kwargs):
        run_id = real(title, steps, *args, **kwargs)
        seen.append((title, [step.kind for step in steps], kwargs.get("kind", runrecord.GENERATION), run_id))
        return run_id

    monkeypatch.setattr(wrapper, "run", spy)
    monkeypatch.setattr("sys.argv", _argv(mysql_config, *SECTOR))
    generate.main()
    [(title, kinds, kind, run_id)] = seen
    assert title == "Sectors (3 unplaced)" and kinds == ["generate-unplaced-sector"] and kind == runrecord.GENERATION
    run = runrecord.read_run(redisqueue.connect(redis_server), run_id)
    assert run["state"] == "done" and run["stages"][0]["done"] == run["stages"][0]["units"] >= 1
    from planetgen.db import store
    conn = store.get_connection(mysql_config)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM sectors").fetchone()["n"] >= 3
    finally:
        conn.close()


def test_a_system_run_with_output_writes_nothing_and_is_a_maintenance_run(mysql_config, monkeypatch, redis_server, capsys):
    seen = []
    real = wrapper.run

    def spy(title, steps, *args, **kwargs):
        seen.append(kwargs.get("kind"))
        return real(title, steps, *args, **kwargs)

    monkeypatch.setattr(wrapper, "run", spy)
    monkeypatch.setattr("sys.argv", _argv(mysql_config, *SYSTEM, "--output", "-"))
    generate.main()
    assert seen == [runrecord.MAINTENANCE]
    assert "Wrapped Star" in capsys.readouterr().out


def test_every_stage_of_the_command_line_is_started_by_the_wrapper():
    """Nothing in `planetgen.generation` or `planetgen.cli` queues a unit or starts a work queue of its own."""
    import pathlib
    import planetgen
    root = pathlib.Path(planetgen.__file__).parent
    offenders = []
    for folder in ("generation", "cli"):
        for path in sorted((root / folder).glob("*.py")):
            text = path.read_text(encoding="utf-8")
            if "WorkQueue(" in text:
                offenders.append(path.name)
    assert offenders == [], f"stages started outside the wrapper: {offenders}"
