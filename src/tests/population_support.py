"""PERF.125: the population pass and the territories pass are two commands; the tests that want a populated galaxy
run both, in order, the way a `galaxy` run's two stages do."""

import contextlib
import threading

from planetgen import tuning
from planetgen.population import model


_CHANCE_LOCK = threading.Lock()
_CHANCE_USERS = {"count": 0, "saved": None}


@contextlib.contextmanager
def every_window_has_a_civilization():
    """Gives a world whose timeline reached the civilization stage a civilization for the length of the block: a request for intelligent life only raises
    the odds (GEN.1130) and a civilization is 1 in 1,000, so a test that needs species sets the chance to 1.

    The setting is one value for the whole process, so blocks that overlap (two passes on two threads, TEST.622) share it: the first in saves the
    chance and the last out puts it back. A block that saved and restored for itself let a pass that finished first restore the chance under the one
    still running, and the one that finished last left the chance at 1 for every test after it."""
    with _CHANCE_LOCK:
        if _CHANCE_USERS["count"] == 0:
            _CHANCE_USERS["saved"] = tuning.CIVILIZATION_CHANCE
            tuning.CIVILIZATION_CHANCE = 1.0
        _CHANCE_USERS["count"] += 1
    try:
        yield
    finally:
        with _CHANCE_LOCK:
            _CHANCE_USERS["count"] -= 1
            if _CHANCE_USERS["count"] == 0:
                tuning.CIVILIZATION_CHANCE = _CHANCE_USERS["saved"]


def populate(conn, rescan=False):
    """`model.run_population` then `model.run_territories` over `conn`; the counts of both in one dict."""
    with every_window_has_a_civilization():
        counts = model.run_population(conn, rescan=rescan)
    return {**counts, **model.run_territories(conn)}
