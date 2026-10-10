# tests/test_steps.py

"""PERF.51: the shared progress step: a bar drawn only for work predicted to pass 15 seconds."""

import queue
import threading
import time

import pytest

from planetgen.generation import run_common, steps
from planetgen.generation.stats import Bucket, GenerationStats


def _stats_with(kind, workers, seconds_per_unit):
    stats = GenerationStats()
    stats.buckets[(kind, workers, 0)] = Bucket(kind, 0, workers=workers, samples=10,
                                               seconds_per_task=seconds_per_unit * 10,
                                               seconds_per_system=seconds_per_unit)
    return stats


@pytest.fixture
def progress():
    with run_common._generation_progress(disable=True) as display:
        yield display


def test_a_step_with_no_recorded_speed_draws_its_bar_at_once(progress):
    with steps.Step("Neighbours", "link", 100, stats=GenerationStats(), progress=progress) as bar:
        assert bar.shown and bar.predicted is None
        assert progress.tasks[0].description == "Neighbours" and progress.main_task == bar.task


def test_a_step_predicted_over_the_limit_draws_at_once(progress):
    stats = _stats_with("link", 1, 1.0)           # 1 s a unit: 100 units = 100 s
    with steps.Step("Neighbours", "link", 100, stats=stats, progress=progress) as bar:
        assert bar.predicted == pytest.approx(100.0) and bar.shown


def test_a_step_predicted_short_draws_nothing(progress):
    stats = _stats_with("link", 1, 0.01)          # 100 units = 1 s
    with steps.Step("Neighbours", "link", 100, stats=stats, progress=progress) as bar:
        assert not bar.shown and not progress.tasks
        bar.update(advance=40)
    assert not progress.tasks


def test_a_short_step_that_runs_past_the_limit_gets_its_bar_when_it_does(progress):
    stats = _stats_with("link", 1, 0.0001)
    with steps.Step("Neighbours", "link", 100, stats=stats, progress=progress, threshold=0.05) as bar:
        bar.update(advance=40)
        assert not bar.shown
        time.sleep(0.2)
        assert bar.shown
        assert progress.tasks[0].completed == 40        # it starts where the work is


def test_the_second_step_shown_is_the_bar_under_the_first(progress):
    with steps.Step("Sectors", "sector", 10, stats=GenerationStats(), progress=progress) as main:
        with steps.Step("Containment", "containment", 5, stats=GenerationStats(), progress=progress) as sub:
            assert progress.main_task == main.task and progress.detail_task == sub.task
        assert progress.detail_task is None and len(progress.tasks) == 1


def test_a_finished_step_records_its_speed_under_its_kind_and_workers():
    stats = GenerationStats()
    with steps.Step("Paths", "paths", 50, stats=stats, workers=3, progress=None) as bar:
        time.sleep(0.02)
        bar.update(advance=50)
    (bucket,) = stats.buckets.values()
    assert (bucket.kind, bucket.workers, bucket.samples) == ("paths", 3, 1)
    assert bucket.seconds_per_system > 0 and stats.pool_rate("paths", 3) > 0


def test_a_failed_or_unrecorded_step_records_nothing():
    stats = GenerationStats()
    with pytest.raises(RuntimeError):
        with steps.Step("Paths", "paths", 50, stats=stats, progress=None):
            raise RuntimeError("boom")
    with steps.Step("Sectors", "sector", 50, stats=stats, progress=None, record=False):
        time.sleep(0.01)
    assert not stats.buckets


def test_a_step_without_a_display_still_works_and_draws_nothing():
    with steps.Step("Paths", "paths", 5, stats=GenerationStats(), progress=None) as bar:
        bar.update(advance=2)
        assert not bar.shown and bar.done == 2


def test_an_unmeasured_step_leaves_no_bar_behind(progress):
    with steps.Step("Backfill (finding)", "backfill", None, stats=GenerationStats(), progress=progress) as bar:
        assert bar.shown
    assert not [task for task in progress.tasks if task.total is None]


def test_the_step_bar_ends_full(progress):
    with steps.Step("Paths", "paths", 8, stats=GenerationStats(), progress=progress):
        pass
    assert progress.tasks[0].completed == 8


def test_a_worker_step_reports_down_a_channel_and_the_relay_draws_it(progress):
    channel = queue.Queue()

    class Pipe:
        def put(self, item):
            channel.put(item)

    with steps.worker_step(Pipe(), "sector 4", "Saving sector 4", "save", 10, workers=2) as bar:
        bar.update(advance=3, total=10)
    relay = steps.Relay(progress=progress)
    items = []
    while not channel.empty():
        items.append(channel.get())
    assert [item[0] for item in items] == ["start", "progress", "end"]
    relay.put(items[0])
    relay.steps["sector 4"].stats = GenerationStats()
    assert relay.steps["sector 4"].shown and progress.tasks[0].description == "Saving sector 4"
    relay.put(items[1])
    assert progress.tasks[0].completed == 3
    relay.put(items[2])
    assert not relay.steps


def test_a_worker_step_without_a_channel_does_nothing():
    with steps.worker_step(None, "k", "name") as bar:
        bar.update(advance=1)
        bar.advance()


def test_a_relay_ends_the_steps_a_dead_worker_left_open(progress):
    relay = steps.Relay(progress=progress)
    relay.put(("start", "k", "Saving", "save", 10, 1))
    assert relay.steps
    relay.close()
    assert not relay.steps


def test_the_relay_drains_a_channel_until_stopped(progress):
    channel = queue.Queue()
    stop = threading.Event()
    relay = steps.Relay(progress=progress)
    channel.put(("start", "k", "Saving", "save", 10, 1))
    channel.put(("end", "k", True))
    thread = threading.Thread(target=relay.drain, args=(channel, stop))
    thread.start()
    time.sleep(0.4)
    stop.set()
    thread.join(timeout=3)
    assert not thread.is_alive() and not relay.steps
