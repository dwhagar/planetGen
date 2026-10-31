"""PERF.125: the population pass and the territories pass are two commands; the tests that want a populated galaxy
run both, in order, the way a `galaxy` run's two stages do."""

import contextlib

from planetgen import tuning
from planetgen.population import model


@contextlib.contextmanager
def every_window_has_a_civilization():
    """Gives a world whose timeline reached the civilization stage a civilization for the length of the block: a request for intelligent life only raises
    the odds (GEN.1130) and a civilization is 1 in 1,000, so a test that needs species sets the chance to 1."""
    saved = tuning.CIVILIZATION_CHANCE
    tuning.CIVILIZATION_CHANCE = 1.0
    try:
        yield
    finally:
        tuning.CIVILIZATION_CHANCE = saved


def populate(conn, rescan=False):
    """`model.run_population` then `model.run_territories` over `conn`; the counts of both in one dict."""
    with every_window_has_a_civilization():
        counts = model.run_population(conn, rescan=rescan)
    return {**counts, **model.run_territories(conn)}
