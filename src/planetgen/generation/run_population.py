# planetgen/generation/run_population.py

"""
Run Population
==============

The `population` command: species, civilizations and territories from
what is already stored (`planetgen.population.model`), also run after a
`sector` or `galaxy` run given `--population`.
"""

from planetgen.cli.stage_progress import StageProgress
from planetgen.queue import work as workQueue
from planetgen.db import store
from planetgen.population import model
from planetgen.util import log



def _population_summary(counts):
    return (f"{counts['new_species']:,} new species; {counts['species']:,} species in all, "
            f"{counts['spacefaring']:,} spacefaring; {counts['polities']:,} polities holding "
            f"{counts['owned_systems']:,} systems.")


def run_population(args):
    """
    Runs the population pass (`model.run_pass`): names the dominant
    species of every new life world, dates civilizations, founds polities
    and recomputes territories. See docs/design/population-and-politics.md.

    Args:
        args (argparse.Namespace): Validated arguments (`command ==
            "population"`).
    """
    conn = store.get_connection(store.mysql_config_from_args(args))
    try:
        with StageProgress(model.pass_stage_count(args.territories_only)) as bar:
            counts = model.run_pass(conn, rescan=args.rescan, territories_only=args.territories_only,
                                    on_stage=bar.stage)
    finally:
        conn.close()
    log.normal(f"Population: {_population_summary(counts)}")


def run_population_after(args):
    """The population pass after a `sector` or `galaxy` run, only with
    `--population` (off by default, Boss 2026-10-01)."""
    if not getattr(args, "population", False):
        return
    conn = store.get_connection(store.mysql_config_from_args(args))
    try:
        with workQueue.job_node("population", "Population pass"):
            with StageProgress(model.pass_stage_count()) as bar:
                counts = model.run_pass(conn, on_stage=bar.stage)
    finally:
        conn.close()
    log.normal(f"Population: {_population_summary(counts)}")
