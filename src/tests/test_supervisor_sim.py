"""OPS.212: the supervisor under a fake clock with faults (`tests/supervisor_sim.py`): killed supervisors, starters that stop, stages whose end was killed, refused Redis commands.
A failing seed replays exactly: `pytest -k seed_17`."""

import pytest

from planetgen.queue import redisqueue
from tests.supervisor_sim import Sim


@pytest.fixture
def real(redis_server, key_prefix, monkeypatch):
    from planetgen.queue import runrecord
    monkeypatch.setattr(runrecord, "PREFIX", key_prefix)
    return redisqueue.connect(redis_server)


@pytest.mark.parametrize("seed", range(12), ids=lambda seed: f"seed_{seed}")
def test_the_supervisor_keeps_its_promises_under_every_fault_the_seed_draws(real, seed):
    sim = Sim(real, seed).run()
    sim.check()
    assert sim.arrivals > 0 and not sim.pending(), f"seed {seed}: still waiting after the faults stopped"


def test_the_simulator_draws_the_faults_it_claims_to(real):
    sim = Sim(real, 3).run()
    assert sim.killed > 0 and sim.conn.refused > 0 and len(sim.ends) > 0
