# tests/test_budget_queries.py

"""
The big-galaxy query budget (TEST.133, opt-in): every page and list query, on a galaxy of millions of systems, must
read about as many rows as the page shows, not as many as the table holds.

Run it with `scripts/test.py budget` (it is skipped everywhere else: building 2,000,000 systems takes about two
minutes and 3 GB, kept between runs in a database named for its size and the schema version). Set
`PLANETGEN_BUDGET_SYSTEMS` and `PLANETGEN_BUDGET_SECTORS` for a smaller galaxy while writing a case, and
`PLANETGEN_BUDGET_REBUILD=1` to build again.

Each case runs the project's own query function for one page and counts the rows the server read for it
(`budget_support.rows_read`). A case fails when that is over its budget, or when the statement limit stops it.
A case that is known to read too much carries `known_over` naming the TODO item that will fix it: it is an
expected failure that turns into a failure of its own when the fix lands, so the case is then given its real
budget. Add a case here with every new page or list query (TEST.135 adds the nebula reads).
"""

import pytest

from planetgen.db import query
from tests import budget_support
from tests.conftest import _test_server_kwargs

pytestmark = [pytest.mark.budget, pytest.mark.fresh_database]

PAGE = 50
STATEMENT_LIMIT_S = 30
PAGE_BUDGET = 5_000
"""int: Rows a page of 50 may read: its own rows, the rows of its joins, and some index slack."""


def _systems():
    return budget_support.sizes()[0]


def _cases():
    """(name, call, budget callable of the system count, known_over or None)."""
    everything = lambda n: n * 2  # noqa: E731 -- a table once, with room for the join to the stars
    return [
        ("systems: page 1 by name", lambda c: query.list_systems(c, limit=PAGE, offset=0, sort="name"), lambda n: PAGE_BUDGET, None),
        ("systems: page 1 by name, star type G", lambda c: query.list_systems(c, star_type_prefix="G", limit=PAGE, offset=0, sort="name"),
         lambda n: PAGE_BUDGET, "a star-type filter walks the name order (PERF.71 section 4)"),
        ("systems: rows 20,000 to 20,050 by name", lambda c: query.list_systems(c, limit=PAGE, offset=20_000, sort="name"),
         lambda n: PAGE_BUDGET, "paging by key instead of OFFSET (PERF.71 section 4)"),
        ("systems: offset half the table by name", lambda c: query.list_systems(c, limit=PAGE, offset=_systems() // 4, sort="name"),
         lambda n: PAGE_BUDGET, "paging by key instead of OFFSET (PERF.71 section 4)"),
        ("systems: page 1 sorted by sector", lambda c: query.list_systems(c, limit=PAGE, offset=0, sort="sector"),
         lambda n: PAGE_BUDGET, None),
        ("systems: page 1 sorted by octant", lambda c: query.list_systems(c, limit=PAGE, offset=0, sort="octant"),
         lambda n: PAGE_BUDGET, None),
        ("systems: count, no filter", lambda c: query.count_systems(c), everything, None),
        ("systems: count, octant I", lambda c: query.count_systems(c, octants=("I",)), everything, None),
        ("systems: count, star type G", lambda c: query.count_systems(c, star_type_prefix="G"), lambda n: n * 4, None),
        ("sectors: count", lambda c: query.count_sectors(c), everything, None),
        ("sectors: page 1 by name", lambda c: query.list_sectors(c, limit=PAGE, offset=0, sort="name"), lambda n: PAGE_BUDGET,
         "keyset paging for the other data tables (PERF.75)"),
        ("sectors: page 1 by systems", lambda c: query.list_sectors(c, limit=PAGE, offset=0, sort="systems"),
         lambda n: PAGE_BUDGET, "the sector sorts read the whole table (PERF.75)"),
        ("sectors: page 1 by density", lambda c: query.list_sectors(c, limit=PAGE, offset=0, sort="density"),
         lambda n: PAGE_BUDGET, "the sector sorts read the whole table (PERF.75)"),
    ]


@pytest.fixture(scope="session")
def big_galaxy():
    server = _test_server_kwargs()
    config = budget_support.ensure_big_galaxy(server)
    yield config


def _case_params():
    params = []
    for name, call, budget, known_over in _cases():
        marks = [pytest.mark.xfail(strict=True, reason=f"reads too many rows until {known_over}")] if known_over else []
        params.append(pytest.param(call, budget, id=name, marks=marks))
    return params


@pytest.mark.parametrize("call, budget", _case_params())
def test_a_page_reads_about_the_rows_it_shows(big_galaxy, call, budget):
    conn = query.open_readonly(big_galaxy, statement_timeout_s=STATEMENT_LIMIT_S)
    try:
        read, seconds = budget_support.rows_read(conn, call)
    finally:
        conn.close()
    allowed = budget(_systems())
    assert read <= allowed, f"read {read:,} rows in {seconds:.1f} s; the budget is {allowed:,}"
