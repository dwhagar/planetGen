"""GEN.1499 slice 1: the forced planets census (`tests/audit/planet_forced_census.py`)."""

from tests.audit import planet_forced_census as census

ROWS = [(1, 3, 2, 0.05, 0.09, 1.0, 2.0, "exchanged", "exchanged"), (2, 1, 0, 0.1, 0.1, 1.0, 1.0, "unchanged", "unchanged"), (3, 2, 2, 0.02, 0.04, 0.5, 1.5, "exchanged", "outside the theory")]


def test_the_summary_counts_the_moved_stars_and_planets():
    s = census.summarize(ROWS)
    assert (s["stars"], s["planets"], s["moved_stars"], s["moved_planets"], s["turned_outside"]) == (3, 6, 2, 4, 1)


def test_the_report_names_the_moved_stars():
    assert "stars whose planets the companion terms move | 2" in census.format_report(ROWS, "h")
