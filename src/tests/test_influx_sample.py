# tests/test_influx_sample.py

"""GEN.1440, GEN.1453: the influx sample (`tests/audit/influx_sample.py`) reads whole systems and reports by the kind of belt they hold."""

from tests.audit import influx_sample as sample


def test_the_report_has_a_row_for_each_kind_of_system():
    text = sample.format_report([("none", 0.0, 0.0, 0, 0), ("inner", 5.0, 10.0, 1, 1), ("inner", 8.0, 8.0, 1, 0), ("outer", 3.0, 3.0, 2, 0)], "h")
    assert "| none belts | 1 | 0 |" in text and "| inner belts | 2 | 2 |" in text and "| outer belts | 1 | 1 |" in text
    assert text.count("| 1 |") >= 1


def test_one_system_gives_rows_for_its_planets_only():
    rows = sample.system_rows(1)
    assert all(row[0] in ("none", "inner", "outer") and row[1] >= 0.0 and row[1] <= row[2] + 1e-9 for row in rows)
