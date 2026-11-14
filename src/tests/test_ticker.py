"""OPS.169: the one progress ticker (`util.ticker`) and the places that use it instead of keeping their own clock."""

import threading

from planetgen import tuning
from planetgen.cli import warm_progress
from planetgen.datasets import fetch
from planetgen.db import block_lights, uncharted_index
from planetgen.util import ticker


class _Clock:
    def __init__(self):
        self.now = 100.0

    def __call__(self):
        return self.now


def test_a_tick_says_nothing_until_the_cadence_has_passed_and_then_once():
    clock, lines = _Clock(), []
    count = ticker.Ticker("Things", total=4000, say=lines.append, clock=clock)
    clock.now += tuning.PROGRESS_LOG_SECONDS - 1
    assert count.tick(10) is False and lines == []
    clock.now += 1
    assert count.tick(1240) is True
    assert lines == ["Things: 1,240 of 4,000 done (31%), 41.3 per s, 30.0 s elapsed."]
    assert count.tick(1300) is False and len(lines) == 1
    clock.now += tuning.PROGRESS_LOG_SECONDS
    count.tick(2000, "3 failed")
    assert lines[-1].startswith("Things: 2,000 of 4,000 done (50%)") and ", 3 failed, " in lines[-1] and lines[-1].endswith("1:00 elapsed.")


def test_the_cadence_is_the_tuning_constant_unless_one_is_given():
    assert ticker.Ticker("A").every == tuning.PROGRESS_LOG_SECONDS
    assert ticker.Ticker("A", every=2.0).every == 2.0


def test_a_ticker_on_a_bars_scale_says_the_share_and_not_a_rate():
    clock = _Clock()
    count = ticker.Ticker("Block light totals", total=1.0, percent=True, clock=clock)
    clock.now += 247
    assert count.text(0.37) == "Block light totals: 37% done, 4:07 elapsed."
    assert ticker.Ticker("Block light totals", total=1.0, percent=True, clock=clock).text(0.5, "still working").endswith("50% done, still working, 0 ms elapsed.")


def test_a_ticker_with_no_total_counts_and_finish_speaks_whatever_the_clock():
    clock, lines = _Clock(), []
    count = ticker.Ticker("Rows", say=lines.append, clock=clock)
    count.finish(12)
    assert lines == ["Rows: 12 done, 0 ms elapsed."]


def test_the_stamp_puts_the_time_since_the_start_after_the_text():
    clock = _Clock()
    began = clock()
    clock.now += 247
    assert ticker.stamp("Adding up the levels...", began, clock) == "Adding up the levels (4:07 elapsed)."
    assert ticker.stamp("Done.", began, clock) == "Done (4:07 elapsed)."


def test_the_pulse_says_every_cadence_while_it_runs_and_survives_a_text_that_raises():
    said, rounds = [], iter([None, "Still working."])
    done = threading.Event()

    def text():
        value = next(rounds)
        if value is None:
            raise RuntimeError("a round that fails")
        return value

    def say(line):
        said.append(line)
        done.set()

    with ticker.Pulse(say, text, every=0.01):
        assert done.wait(5)
    assert said == ["Still working."]


def test_no_place_keeps_its_own_progress_clock():
    """The block lights, the index, the downloads and the warm-up say where they are at the one cadence."""
    assert not hasattr(block_lights, "PROGRESS_SAY_SECONDS") and not hasattr(uncharted_index, "SAY_SECONDS")
    assert fetch.PROGRESS_SECONDS == tuning.PROGRESS_LOG_SECONDS == warm_progress.PROGRESS_SECONDS

