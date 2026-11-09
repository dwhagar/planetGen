"""GEN.1493 slice 1: the late accretion census (`tests/audit/late_accretion_census.py`)."""

from tests.audit import late_accretion_census as census


def test_the_summary_counts_the_floor_and_the_ceiling_and_the_time():
    rows = [(1, 1.0, 0.05, 0.5, True, False), (2, 1.0, 3.0, 0.5, False, True), (3, 1.0, 0.5, 0.5, False, False), ("time", 1.0, 10.0)]
    s = census.summarize(rows)
    assert (s["planets"], s["floor"], s["cap"], s["median"]) == (3, 1, 1, 1.0) and s["pass_s"] == 1.0 and s["total_s"] == 10.0


def test_the_report_names_the_pass_cost():
    rows = [(1, 1.0, 0.5, 0.5, False, False), ("time", 2.0, 10.0)]
    assert "20 percent is the pass" in census.format_report(rows, "h")
