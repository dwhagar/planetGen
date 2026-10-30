# planetgen/generation/run_population.py

"""
Run Population
==============

The `population` and `territories` commands (PERF.125): the population pass (species, civilizations and polities) and
the territories pass (which polity owns which system), each from what is already stored (`planetgen.population.model`)
and each its own registered command, `run-population` and `run-territories`, of one unit. A `galaxy` run asks for the
two as consecutive stages (`run_population_stage`, `run_territories_stage`); neither knows the other.
"""

from planetgen.generation import followed, stages, steps
from planetgen.db import store
from planetgen.generation import run_plan
from planetgen.population import model
from planetgen.queue import commands
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


class PassFollower(followed.Follower):
    """One pass as a stage of a run on the queue: followed from here with a bar, each line the unit returns logged."""

    def __init__(self, args, key, label, step_kind, extra):
        super().__init__(args)
        self.key, self.label, self.step_kind, self.extra = key, label, step_kind, extra

    def params(self):
        return dict(self.extra)

    def on_plan(self, plan):
        self.hold(steps.Step(self.label, self.step_kind, 1, args=self.args, progress=self.progress))

    def on_result(self, unit, result):
        for line in result["lines"]:
            log.normal(line)


def population_follower(args):
    """The follower of the population stage of a run."""
    return PassFollower(args, "population", "Population pass", "population", {"rescan": bool(getattr(args, "rescan", False))})


def territories_follower(args):
    """The follower of the territories stage of a run."""
    return PassFollower(args, "territories", "Territories pass", "territories", {})


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
    """The population pass as a run of its own (`planetgen population`)."""
    stages.enter("population")
    followed.run_planned(args, planned_population(args))


def planned_population(args):
    """The population pass as a `followed.Planned`."""
    return followed.Planned([population_follower(args)], "Population pass")


def run_territories_stage(args):
    """The territories pass as a run of its own (`planetgen territories`)."""
    stages.enter("territories")
    followed.run_planned(args, planned_territories(args))


def planned_territories(args):
    """The territories pass as a `followed.Planned`."""
    return followed.Planned([territories_follower(args)], "Territories pass")
