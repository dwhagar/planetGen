"""GEN.1493 slice 0: the late supply census (`tests/audit/late_supply_census.py`)."""

from tests.audit import late_supply_census as census


def test_the_summary_leaves_out_planets_that_take_nothing_from_the_median():
    rows = [(1, 1.0, 1.0, True, False, 0.0), (2, 1.0, 1.0, True, False, 2.0), (3, 1.0, 1.0, True, True, 4.0)]
    table = census.summarize(rows)
    assert table["all"][:3] == (3, 1, 3.0) and table["with a giant"][:3] == (1, 0, 4.0) and table["without a giant"][:3] == (2, 1, 2.0)


def test_the_report_names_the_share_that_matches_the_present_law():
    rows = [(i, 1.0, 1.0, True, False, 10.0) for i in range(12)]
    assert "0.1" in census.format_report(rows, "h")
