"""OPS.266: a sector's generation tells the unit that runs it what it is making (`generation/bodybeat.py`): the system under way as the second-level count's label, each body made as one more on it, a finished
phase of generation as progress, and the save of the sector as its own count."""

from planetgen.generation import bodybeat
from planetgen.queue import redisqueue, unitbeat
from planetgen.util import log


def _beat(redis_server, now):
    return unitbeat.UnitBeat(redisqueue.connect(redis_server), "bodies", 0, "generate-sectors", clock=lambda: now[0])


def test_each_body_made_counts_under_the_label_of_its_system(redis_server):
    now = [0.0]
    with _beat(redis_server, now) as beat:
        bodybeat.system("System 2 of 40: bodies made")
        for _ in range(3):
            now[0] += 4.0
            bodybeat.made()
        assert (beat.record()["step"], beat.record()["step_done"]) == ("System 2 of 40: bodies made", 3.0)
        assert beat.longest_gap() == 4.0


def test_a_phase_of_generation_that_ends_is_progress_and_a_save_is_a_count_of_its_own(redis_server):
    now = [0.0]
    with _beat(redis_server, now) as beat:
        with log.timed_phase("settle pass"):
            now[0] += 3.0
        now[0] += 3.0
        bodybeat.settled()
        saved = bodybeat.SaveStep(10)
        saved.update(completed=0)
        for _ in range(4):
            saved.advance()
        record = beat.record()
        assert (record["step"], record["step_done"], record["step_total"]) == ("Saving the sector", 4.0, 10.0)
        assert beat.longest_gap() == 3.0


def test_outside_a_unit_every_call_does_nothing():
    bodybeat.system("x")
    bodybeat.made()
    bodybeat.settled()
    saved = bodybeat.SaveStep(3)
    saved.update(completed=1, total=3)
    saved.advance()
