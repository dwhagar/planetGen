# tests/test_step_timer.py

"""`util/step_timer.py` (PERF.163): the steps of a piece of work with the seconds each took, as a table."""

import pytest

from planetgen.util.step_timer import StepTimer


def test_a_step_is_recorded_after_its_block_and_the_total_is_their_sum():
    timer = StepTimer()
    with timer.step("first"):
        pass
    with timer.step("second"):
        sum(range(10000))
    assert [label for label, _seconds in timer.steps] == ["first", "second"]
    assert timer.total() == pytest.approx(sum(seconds for _label, seconds in timer.steps))


def test_a_block_that_raises_leaves_no_step_behind():
    timer = StepTimer()
    with pytest.raises(ValueError):
        with timer.step("broken"):
            raise ValueError("no")
    assert timer.steps == []


def test_the_table_lists_the_costliest_step_first_with_its_share():
    timer = StepTimer()
    timer.steps = [("small", 1.0), ("big", 3.0)]
    assert timer.lines("Work") == ["Work (4.00s):", "      3.00s  75.0%  big", "      1.00s  25.0%  small"]
    assert StepTimer().lines("Nothing") == ["Nothing (0.00s):"]
