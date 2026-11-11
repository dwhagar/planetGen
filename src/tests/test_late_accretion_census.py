"""GEN.1493 slice 1: the late accretion census (`tests/audit/late_accretion_census.py`)."""

from tests.audit import late_accretion_census as census


def _rows():
    return [(1, 1.0, 0.05, 0.5, True, False, 0.0, 0.0), (2, 1.0, 3.0, 0.5, False, True, 0.02, 0.01), (3, 1.0, 0.5, 0.5, False, False, 0.0, 0.0), ("time", 1.0, 2.0, 10.0, 1, 2)]


def test_the_summary_counts_the_floor_the_ceiling_the_changes_and_the_time():
    s = census.summarize(_rows())
    assert (s["planets"], s["floor"], s["cap"], s["median"], s["same"], s["over_one_percent"]) == (3, 1, 1, 1.0, 2, 1)
    assert s["dp_max"] == 0.02 and s["dt_max"] == 0.01 and (s["pass_s"], s["full_s"], s["total_s"], s["made"], s["applies"]) == (1.0, 2.0, 10.0, 1, 2)


def test_the_report_names_the_pass_cost_and_the_difference_from_the_full_replay():
    text = census.format_report(_rows(), "h")
    assert "10.0 percent of the whole systems" in text and "1 of 2" in text and "0.02" in text
