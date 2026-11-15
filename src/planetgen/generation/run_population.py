# planetgen/generation/run_population.py

"""
Run Population
==============

The `population` and `territories` commands (PERF.125): the population pass (species, civilizations and polities) and
the territories pass (which polity owns which system), each from what is already stored (`planetgen.population.model`)
and each registered commands of one unit: the population pass is four (`scan-species`, `date-civilizations`, `refresh-technology`, `count-population`, OPS.239) and the territories pass is `run-territories`. A `galaxy` run asks for
the two as consecutive stages (`run_population_stage`, `run_territories_stage`); neither knows the other.
"""

from planetgen.generation import followed, stages, steps
from planetgen.population import model
from planetgen.queue import commands, engine, unitbeat
from planetgen.util import log


def _territories_summary(counts):
    return f"{counts['polities']:,} polities holding {counts['owned_systems']:,} systems."


def expand_population(params):
    """The stage of each of the population pass's commands (`scan-species`, `date-civilizations`, `refresh-technology`, `count-population`): one unit, whose pieces the engine runs (`run-territories` cuts its stage the same way)."""
    return commands.Expansion([(params, 1.0)], {})


expand_territories = expand_population


def _reader(params):
    """A connection to read with while the pieces are listed, closed by the caller (the unit's own opens on its first piece)."""
    return engine.Context("scan-species", params, None)


def scan_species_pieces(params, cursor):
    """
    The `scan-species` command's pieces (OPS.203, OPS.239; the first of the four commands the `run-population` command was): forget every species for a rescan, then scan the planets one batch of ids a
    piece (`model.population_pieces`, the pieces of kinds `reset` and `scan`).
    """
    reader = _reader(params)
    try:
        listed = model.population_pieces(reader.db, rescan=bool(params.get("rescan")))
    finally:
        reader.close()
    return [piece for piece in listed if piece[0] in ("reset", "scan")][cursor:]


def date_civilizations_pieces(params, cursor):
    """The `date-civilizations` command's pieces (OPS.239): drop the species without a civilization, then date the civilizations and found or dissolve the polities."""
    return [["drop"], ["civilizations"]][cursor:]


def refresh_technology_pieces(params, cursor):
    """The `refresh-technology` command's pieces (OPS.239): refresh the map's technology rows (MAP.344), one piece."""
    return [["technology"]][cursor:]


def count_population_pieces(params, cursor):
    """The `count-population` command's pieces (OPS.239): count the species, the spacefaring ones and the polities, one piece."""
    return [["counts"]][cursor:]


def _under_the_lock(ctx, work):
    """Runs `work()` holding the pass lock (`model.take_pass_lock`; the beat says the unit waits on `mysql` while it waits), commits through `ctx.commit()` before the
    lock is let go, and returns what `work` returns."""
    with unitbeat.waiting("mysql"):
        model.take_pass_lock(ctx.db)
    try:
        ctx.commit()        # a fresh transaction, so the piece reads what the last holder of the lock committed
        result = work()
        ctx.commit()
        return result
    finally:
        model.release_pass_lock(ctx.db)


def population_do(params, piece, ctx):
    """One piece of the population pass (`model.do_population_piece`), under the pass lock and committed through `ctx.commit()`: the `do` of all four commands of the pass."""
    return _under_the_lock(ctx, lambda: model.do_population_piece(ctx.db, piece))


def scan_species_gather(params, results, ctx):
    """The `scan-species` unit's result: `added`, the species the scan found (and `lines`, what the run says about it)."""
    added = sum(result.get("added", 0) for result in results)
    return {"added": added, "lines": [f"Population: {added:,} new species scanned."]}


def date_civilizations_gather(params, results, ctx):
    """The `date-civilizations` unit's result: `dated`, always true."""
    return {"dated": True}


def refresh_technology_gather(params, results, ctx):
    """The `refresh-technology` unit's result: `refreshed`, always true."""
    return {"refreshed": True}


def count_population_gather(params, results, ctx):
    """
    What the `count-population` unit returns: the counts of the pass (`species`, `spacefaring`, `polities`) and `lines`, what the run says about it (this unit may run in a worker, so the process that
    started the run logs them).
    """
    counts = next(result for result in reversed(results) if "species" in result)
    return {**counts, "lines": [f"Population: {counts['species']:,} species in all, {counts['spacefaring']:,} spacefaring; {counts['polities']:,} polities."]}


def territories_pieces(params, cursor):
    """The `run-territories` command's pieces (OPS.203): recompute which polity owns which system (`model.territories_pieces`)."""
    return model.territories_pieces()[cursor:]


def territories_do(params, piece, ctx):
    """One piece of the territories pass (`model.do_territories_piece`), counted through `ctx.step` and `ctx.advance`, under the pass lock and committed through
    `ctx.commit()`."""
    return _under_the_lock(ctx, lambda: model.do_territories_piece(ctx.db, piece, ctx))


def territories_gather(params, results, ctx):
    """
    What the unit returns: `polities`, `owned_systems` and `lines`.
    """
    counts = next(result for result in reversed(results) if "polities" in result)
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


def population_followers(args):
    """The followers of the population stages of a run: the scan, the dating of the civilizations, the technology rows and the count (OPS.239)."""
    rescan = {"rescan": bool(getattr(args, "rescan", False))}
    return [PassFollower(args, "species", "Scan the planets for species", "population", rescan),
            PassFollower(args, "civilizations", "Date the civilizations", "population", {}),
            PassFollower(args, "technology", "Refresh the technology rows", "population", {}),
            PassFollower(args, "population", "Count the population", "population", {})]


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
    followed.run_planned(args, planned_population(args))


def planned_population(args):
    """The population pass as a `followed.Planned`."""
    return followed.Planned(population_followers(args), "Population pass")


def run_territories_stage(args):
    """The territories pass as a run of its own (`planetgen territories`)."""
    stages.enter("territories")
    followed.run_planned(args, planned_territories(args))


def planned_territories(args):
    """The territories pass as a `followed.Planned`."""
    return followed.Planned([territories_follower(args)], "Territories pass")
