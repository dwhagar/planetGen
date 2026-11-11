"""GEN.1492 slice 0: the belt forcing census (`tests/audit/belt_forcing_census.py`)."""

from tests.audit import belt_forcing_census as census


def test_the_summary_counts_the_belts_outside_the_theory_the_forced_and_the_diluted():
    rows = [(1, 0.1, 0.1, 0.2, 0.05, False), (2, 0.1, 0.1, 0.05, 0.05, False), (3, 0.5, 0.1, 0.0, 0.0, True)]
    s = census.summarize(rows)
    assert (s["belts"], s["held"], s["reached"], s["over_rms"], s["diluted"]) == (3, 1, 2, 1, 1)
    assert s["median"] == 1.25


def test_the_report_names_the_forced_share():
    rows = [(1, 0.1, 0.1, 0.2, 0.2, False)]
    assert "forcing larger than the drawn rms | 1 (100 percent)" in census.format_report(rows, "h")
