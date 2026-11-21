"""A command that ends on SIGTERM says who sent it (stop_reason), as the warm-up does."""

import os
import signal

import pytest

from planetgen.admin import stop_reason
from planetgen.util import process


@pytest.fixture
def term_handler(tmp_path, monkeypatch):
    monkeypatch.setattr(stop_reason, "folder", lambda: str(tmp_path / "stop-reasons"))
    previous = process.interrupt_on_sigterm()
    yield
    signal.signal(signal.SIGTERM, previous)


def test_a_reason_a_canceller_left_is_carried_by_the_interrupt(term_handler):
    stop_reason.leave(os.getpid(), "update.sh pulled a new build")
    with pytest.raises(process.Stopped) as caught:
        os.kill(os.getpid(), signal.SIGTERM)
    assert caught.value.reason == "update.sh pulled a new build"
    assert isinstance(caught.value, KeyboardInterrupt)


def test_a_signal_with_no_reason_says_so(term_handler):
    with pytest.raises(process.Stopped) as caught:
        os.kill(os.getpid(), signal.SIGTERM)
    assert caught.value.reason == process.NO_REASON
