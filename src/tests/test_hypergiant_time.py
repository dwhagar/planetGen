# tests/test_hypergiant_time.py

"""
PERF.249: a hypergiant system with every orbit-producing flag forced on is made in minutes of a worker's time at most, not in the 5 to 11 minutes and 8.4 GB the first bug
hunt measured (`band_climate` split the record of a fast-spinning world into one window per solar day: part 1, PR 5634; the fixed point of the seasonal swing formed its
whole Gram matrix every pass: part 2).
"""

import time

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem

from tests.fuzz_support import deterministic_entropy

HYPERGIANT_SECONDS = 240.0
"""float: What the K00 system may take. It takes about 30 s alone on a four-core worker; it took 311 s before PERF.249 and the limit of a test is 900 s."""


def test_a_hypergiant_system_with_every_orbit_flag_forced_is_made_in_minutes():
    config = SystemConfig()
    config.STAR_TYPE, config.MOONS, config.COMETS, config.MAX_PLANETS = "K00", True, True, True
    started = time.perf_counter()
    with deterministic_entropy(0):
        system = StarSystem(system_config=config)
    seconds = time.perf_counter() - started
    assert system.all_planets
    assert seconds < HYPERGIANT_SECONDS, f"a K00 system took {seconds:.0f} s"
