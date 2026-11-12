# tests/dead_fixtures.py

"""
A fixture that hung once is not tried again for the next module (TEST.190).

A session fixture that fails is raised again for every test that needs it at once (pytest keeps the outcome), but a fixture of a smaller scope is
set up again by each module (or class) that needs it. The browser tests' server is one: when it does not come up, and the setup waits instead of
failing, every module waits the whole `timeout` (900 s) in turn, and one dead server costs a run hours.

This plugin remembers a module-, class- or package-scoped fixture whose setup was stopped by `pytest-timeout`, and fails every later setup of the
same fixture in the same process at once, naming the test that waited first. A function-scoped fixture is left alone (one slow test is not a dead
server). `PLANETGEN_TEST_RETRY_DEAD_FIXTURES=1` turns it off. Loaded with `-p tests.dead_fixtures` (`pytest.ini`), like `tests.hang_guard`.
"""

import os

import pytest

OFF_ENV = "PLANETGEN_TEST_RETRY_DEAD_FIXTURES"
"""str: Set to 1 to set a hung fixture up again for every module, as pytest does."""

TIMED_OUT = "Timeout (>"
"""str: What the `Failed` that `pytest-timeout` raises starts with in its message."""

_dead = {}
"""dict: `{fixture function: the id of the test whose setup it hung in}`, for the fixtures of a scope wider than one test."""


def _key(fixturedef):
    """The function behind a fixture, so a fixture a helper module defines and several test modules import (each gets a definition of its own) is one."""
    func = getattr(fixturedef.func, "__wrapped__", fixturedef.func)
    return f"{func.__module__}.{func.__qualname__}"


def dead():
    """The fixtures remembered as hung in this process."""
    return dict(_dead)


def _wide(fixturedef):
    return fixturedef.scope != "function" and not os.environ.get(OFF_ENV)


@pytest.hookimpl(tryfirst=True)
def pytest_runtest_setup(item):
    """Before any fixture of the test is set up (a failure raised inside a fixture's own setup would leave pytest's bookkeeping of it half done)."""
    if not _dead or os.environ.get(OFF_ENV):
        return
    for name in item.fixturenames:
        definitions = item._fixtureinfo.name2fixturedefs.get(name)
        key = _key(definitions[-1]) if definitions else None
        if key in _dead and _wide(definitions[-1]):
            pytest.fail(f"TEST.190: fixture {name!r} hung in the setup for {_dead[key]} and is not set up again "
                        f"(set {OFF_ENV}=1 to try it for every module)", pytrace=False)


@pytest.hookimpl(hookwrapper=True)
def pytest_fixture_setup(fixturedef, request):
    outcome = yield
    if _wide(fixturedef) and outcome.excinfo is not None and TIMED_OUT in str(outcome.excinfo[1]):
        _dead.setdefault(_key(fixturedef), request.node.nodeid)
