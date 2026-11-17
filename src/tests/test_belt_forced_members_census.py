"""GEN.1492 slice 1: the forced members census (`tests/audit/belt_forced_members_census.py`)."""

from tests.audit import belt_forced_members_census as census


def test_the_summary_counts_the_moved_belts_and_members():
    rows = [(1, 10, 8, 0.1, 0.2, True, 6, 5.0, 9.0), (2, 5, 0, 0.1, 0.1, False, 0, 5.0, 5.0), (3, 20, 20, 0.1, 0.12, False, 20, 4.0, 6.0)]
    s = census.summarize(rows)
    assert (s["belts"], s["moved_belts"], s["members"], s["moved_members"], s["over_half"]) == (3, 2, 35, 28, 1)
    assert (s["tilted_belts"], s["tilted_members"]) == (2, 26)


def test_the_report_names_the_moved_share():
    assert "belts whose members the forcing moves | 1" in census.format_report([(1, 10, 8, 0.1, 0.2, True, 6, 5.0, 9.0)], "h")


def test_the_report_names_the_tilted_belts():
    assert "belts whose members the tilt moves | 1" in census.format_report([(1, 10, 8, 0.1, 0.2, True, 6, 5.0, 9.0)], "h")
