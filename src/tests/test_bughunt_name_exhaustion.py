# tests/test_bughunt_name_exhaustion.py

"""
Tier 1 bug-hunt coverage: `utils.generate_phoneme_salad_name`'s retry loop
under name-validity exhaustion.

Confirmed (via a real 5-second `timeout`-wrapped subprocess, before any
fix) that this function hung forever -- an unconditional `while True` with
`is_name_valid` mocked to always return False never terminated. Fixed by
capping the loop at `utils.MAX_NAME_GENERATION_ATTEMPTS` and raising a
clean `RuntimeError` on exhaustion (see that constant's own docstring).
Every star/planet/moon/sector name in this generator funnels through this
one function, so an unbounded hang here is an unbounded hang for the
entire generator -- this test (now passing) is the permanent regression
guard, using an in-process hard wall-clock timeout (`signal.alarm`, POSIX-
only, matching this project's own Linux/Apache deployment target -- see
`README.md`'s "Web Interface" section) rather than a real subprocess, so
it runs as fast as any other test while still proving termination, not
just "eventually raises" with no time bound.
"""

import signal
import time
from unittest import mock

import pytest

from stellarObjects import utils


class _AlarmTimeout(Exception):
    pass


def _run_with_wall_clock_timeout(seconds, fn):
    """Runs `fn()` under a hard SIGALRM timeout, raising `_AlarmTimeout`
    if it doesn't return in time -- proves an upper bound on how long a
    call can run, not just that it eventually raises *something*."""
    def _handler(signum, frame):
        raise _AlarmTimeout(f"did not return within {seconds}s")

    old_handler = signal.signal(signal.SIGALRM, _handler)
    signal.alarm(seconds)
    try:
        return fn()
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old_handler)


def test_phoneme_salad_name_raises_cleanly_when_every_candidate_invalid():
    """Regression: previously hung forever (confirmed via a real 5s
    subprocess timeout before the fix); now bounded by
    MAX_NAME_GENERATION_ATTEMPTS and raises RuntimeError."""
    with mock.patch.object(utils, "is_name_valid", return_value=False):
        with pytest.raises(RuntimeError):
            _run_with_wall_clock_timeout(
                10,
                lambda: utils.generate_phoneme_salad_name(
                    ["test", "testing", "example"], ["pre"], ["suf"]
                ),
            )


def test_phoneme_salad_name_exhaustion_is_fast_not_just_bounded():
    """The cap isn't just "eventually terminates" -- 10,000 attempts of
    cheap string validation should finish in well under a second, so a
    real misconfiguration fails fast during generation, not after a
    multi-second stall."""
    with mock.patch.object(utils, "is_name_valid", return_value=False):
        start = time.monotonic()
        with pytest.raises(RuntimeError):
            utils.generate_phoneme_salad_name(["test", "testing", "example"], ["pre"], ["suf"])
        elapsed = time.monotonic() - start
    assert elapsed < 5.0, f"exhaustion took {elapsed}s, expected well under 5s"


def test_phoneme_salad_name_still_succeeds_normally():
    """Sanity check the fix didn't change normal (non-exhausted)
    behavior -- a real name pool with is_name_valid unmocked must still
    return quickly, well within the cap."""
    name = utils.generate_phoneme_salad_name(
        ["aurora", "nebula", "zenith", "corvus"], ["ka", "el", "or"], ["us", "ia", "on"]
    )
    assert isinstance(name, str) and len(name) > 0


def test_sector_name_generator_also_bounded_under_exhaustion():
    """generate_sector_name joins two generate_phoneme_salad_name calls --
    confirms the fix propagates to it too, not just the lower-level
    function directly."""
    with mock.patch.object(utils, "is_name_valid", return_value=False):
        with pytest.raises(RuntimeError):
            _run_with_wall_clock_timeout(10, utils.generate_sector_name)
