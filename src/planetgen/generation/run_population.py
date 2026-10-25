# planetgen/generation/run_population.py

"""
Run Population
==============

The `population` and `territories` commands (PERF.125): the population pass (species, civilizations and polities) and
the territories pass (which polity owns which system), each from what is already stored (`planetgen.population.model`)
and each its own registered command, `run-population` and `run-territories`, of one unit. A `galaxy` run asks for the
two as consecutive stages (`run_population_stage`, `run_territories_stage`); neither knows the other.
"""

from planetgen.generation import run_common, stages, steps
from planetgen.db import store
from planetgen.generation import run_plan
from planetgen.population import model
from planetgen.queue import commands, units, work as workQueue, wrapper
from planetgen.util import log


def _population_summary(counts):
    return (f"{counts['new_species']:,} new species; {counts['species']:,} species in all, "
            f"{counts['spacefaring']:,} spacefaring; {counts['polities']:,} polities.")


def _territories_summary(counts):
    return f"{counts['polities']:,} polities holding {counts['owned_systems']:,} systems."


def expand_population(params):
    """The `run-population` command's stage: one unit, the whole pass."""
    return commands.Expansion([(params, 1.0)], {})


def population_unit(params, report):
    """
    The `run-population` command's unit: scan the life worlds, drop the species without a civilization, date the
    civilizations and found or dissolve the polities (`model.run_population`).

    Returns:
        dict: `new_species`, `species`, `spacefaring`, `polities` and `lines`, what the run says about it (this unit may
            run in a worker, so the process that started the run logs them).
    """
    conn = store.get_connection(run_plan._mysql_config_of(params))
    try:
        counts = model.run_population(conn, rescan=bool(params.get("rescan")))
    finally:
        conn.close()
    return {**counts, "lines": [f"Population: {_population_summary(counts)}"]}


def expand_territories(params):
    """The `run-territories` command's stage: one unit, the whole pass."""
    return commands.Expansion([(params, 1.0)], {})


def territories_unit(params, report):
    """
    The `run-territories` command's unit: recompute which polity owns which system from the polities stored
    (`model.run_territories`).

    Returns:
        dict: `polities`, `owned_systems` and `lines`.
    """
    conn = store.get_connection(run_plan._mysql_config_of(params))
    try:
        counts = model.run_territories(conn)
    finally:
        conn.close()
    return {**counts, "lines": [f"Territories: {_territories_summary(counts)}"]}


def _run_pass(args, key, label, step_kind, extra):
    """One pass as the one stage of a run on the queue: the step for stage `key` (`wrapper.step_for`), followed from here
    with a bar, each line the unit returns logged."""
    mysql_config = store.mysql_config_from_args(args)
    workers = run_common._worker_count(args)
    step, listing = wrapper.step_for(args, key, label, {"workers": workers, "database": wrapper.database_args(mysql_config),
                                                         **extra})
    row = stages.hand_over_open_stage()

    def on_result(_index, _unit, result):
        for line in result["lines"]:
            log.normal(line)

    with run_common._generation_progress() as progress:
        log.set_console(progress.console)
        try:
            with steps.Step(label, step_kind, 1, args=args, progress=progress):
                wrapper.run(label, [step], lambda _index, plan: None, on_result, workers_count=workers,
                            who=getattr(args, "who", None), database=mysql_config.database,
                            argv=getattr(args, "argv", None), seed=None, listing=listing,
                            rows=None if row is None else {key: row})
        finally:
            log.reset_console()


def run_population(args):
    """
    `planetgen population` (PERF.125): the population pass, and nothing else.
    See docs/design/habitability/population-and-politics.md.

    Args:
        args (argparse.Namespace): Validated arguments (`command == "population"`).
    """
    stages.begin(stages.for_command("population", args), args, "population")
    run_population_stage(args)
    stages.finish()


def run_territories(args):
    """`planetgen territories` (PERF.125): the territories pass, and nothing else."""
    stages.begin(stages.for_command("territories", args), args, "territories")
    run_territories_stage(args)
    stages.finish()


def run_population_stage(args):
    """The population pass as a stage of a run (a `galaxy` run's, or `planetgen population`'s own)."""
    stages.enter("population")
    with workQueue.job_node("population", "Population pass"):
        _run_pass(args, "population", "Population pass", "population", {"rescan": bool(getattr(args, "rescan", False))})


def run_territories_stage(args):
    """The territories pass as a stage of a run."""
    stages.enter("territories")
    with workQueue.job_node("territories", "Territories pass"):
        _run_pass(args, "territories", "Territories pass", "territories", {})
