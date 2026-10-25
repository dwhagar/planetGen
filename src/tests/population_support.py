"""PERF.125: the population pass and the territories pass are two commands; the tests that want a populated galaxy
run both, in order, the way a `galaxy` run's two stages do."""

from planetgen.population import model


def populate(conn, rescan=False):
    """`model.run_population` then `model.run_territories` over `conn`; the counts of both in one dict."""
    counts = model.run_population(conn, rescan=rescan)
    return {**counts, **model.run_territories(conn)}
