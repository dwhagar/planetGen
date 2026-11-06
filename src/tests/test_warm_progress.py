"""The progress `planetgen warm --progress` shows for a stage cut into chunks (cli/warm_progress.py)."""

from planetgen.cli import warm_progress


def _chunk(warmed=64, built=64, seconds=38.0):
    return {"warmed": warmed, "built": built, "cached": warmed - built, "left_out": 0, "noun": "nodes",
            "seconds": seconds, "line": "unused"}


def _run(chunks, show=True, each=40.0, own=38.0, cost=None, expected=None):
    lines, now = [], [0.0]

    def clock():
        return now[0]

    progress = warm_progress.WarmProgress(["lists", "galaxy-deep-nodes"], lines.append, show=show, clock=clock,
                                          expected=expected)
    progress.plan(1, {"units": chunks, "cost": cost})
    for _ in range(chunks):
        now[0] += each
        progress.result(_chunk(seconds=own))
    return lines


def test_a_progress_line_comes_when_a_minute_and_ten_percent_have_passed():
    """OPS.99: a chunk every 40 s of a stage of 8 is a line every two chunks (80 s, 25 points), not one per chunk."""
    lines = _run(8)
    progress = [line for line in lines if "%" in line]
    assert progress[0] == "galaxy-deep-nodes 25% (2 of 8 chunks), 1:20 elapsed"
    assert len(progress) == 3 and progress[2].startswith("galaxy-deep-nodes 75% (6 of 8 chunks)")
    assert lines[-1] == ("galaxy-deep-nodes: 512 nodes warmed (512 built, 0 already stored) in 5:20 over 8 chunks of up to 64 nodes; "
                         "1 chunk at once, 38.0 s a chunk on its worker; step 2 of 2")


def test_no_time_left_is_shown_from_fewer_than_a_few_chunks():
    """PERF.237: the first chunks are the first of the workers that run side by side; they say nothing of the rest."""
    lines = _run(8)
    first = next(line for line in lines if "%" in line)
    assert "left" not in first and "2 of 8" in first
    later = [line for line in lines if "%" in line][-1]
    assert "4:00 elapsed, about 1:20 left" in later


def test_the_closing_line_says_how_many_chunks_run_at_once_and_how_long_a_chunk_takes_on_its_worker():
    """PERF.237: ten workers each take 37 s a chunk, so a chunk is done every 3.7 s overall. The progress line leaves that out; the stage's last line has it."""
    lines = _run(40, each=3.7, own=37.0)
    assert lines[-1].endswith("; 10 chunks at once, 37.0 s a chunk on its worker; step 2 of 2")
    assert not any("at once" in line for line in lines[:-1])


def test_a_long_stage_leaves_a_bounded_number_of_lines():
    """A stage of 180 chunks that takes ten minutes is a line a minute at most, and a stage that takes hours one at each tenth."""
    quick = _run(180, each=3.0)
    assert 8 <= len([line for line in quick if "%" in line]) <= 10
    slow = _run(180, each=60.0)
    progress = [line for line in slow if "%" in line]
    assert 8 <= len(progress) <= 10 and progress[0].startswith("galaxy-deep-nodes 10% (18 of 180 chunks)")
    assert "180 chunks" in slow[-1]


def test_a_stage_of_one_unit_prints_its_own_line_once():
    lines = []
    progress = warm_progress.WarmProgress(["lists", "search"], lines.append, show=True, expected={"lists": 1.0, "search": 1.0})
    progress.plan(0, {"units": 1, "cost": 1.0})
    progress.result({"line": "lists: 20 pages warmed"})
    assert lines == ["lists: 20 pages warmed"]


def test_without_progress_only_the_summary_is_left():
    lines = _run(4, show=False)
    assert len(lines) == 1 and "over 4 chunks" in lines[0]


def test_the_overall_bar_weighs_a_step_by_its_expected_cost():
    """PERF.234: a long step moves the bar as far as it takes. The first of three steps is a tenth of the cost, then the rest."""
    progress = warm_progress.WarmProgress(["lists", "systems", "search"], lambda line: None, show=True,
                                          expected={"lists": 10.0, "systems": 980.0, "search": 10.0})
    progress.plan(0, {"units": 1, "cost": 10.0})
    assert progress._weight_done() == 0.0
    progress.plan(1, {"units": 10, "cost": 980.0})
    assert progress._weight_done() == 10.0
    progress.done = 5
    assert progress._weight_done() == 10.0 + 490.0


def test_a_queued_steps_cost_replaces_its_expected_cost():
    progress = warm_progress.WarmProgress(["lists", "systems"], lambda line: None, show=True,
                                          expected={"lists": 10.0, "systems": 100.0})
    assert progress.weights == [10.0, 100.0]
    progress.plan(0, {"units": 1, "cost": 4.0})
    assert progress.weights == [4.0, 100.0]
    progress.plan(1, {"units": 20, "cost": 575.0})
    assert progress.weights == [4.0, 575.0]


def test_the_whole_warm_ups_time_left_is_a_rough_guess_until_every_step_is_queued():
    now = [0.0]
    lines = []
    progress = warm_progress.WarmProgress(["systems", "search"], lines.append, show=True, clock=lambda: now[0],
                                          expected={"systems": 100.0, "search": 100.0})
    progress.plan(0, {"units": 20, "cost": 100.0})
    for _ in range(warm_progress.MIN_CHUNKS):
        now[0] += 10.0
        progress.result(_chunk(warmed=250, seconds=10.0))
    assert progress.rough()
    assert progress.overall_left() == 50.0 / 25.0 * (200.0 - 25.0)     # elapsed over cost done, times the cost to go
    assert "(rough guess)" in progress._overall_left_text()
    progress.plan(1, {"units": 1, "cost": 100.0})
    assert not progress.rough() and "(rough guess)" not in progress._overall_left_text()


def test_the_whole_warm_up_shows_no_time_left_from_fewer_than_a_few_chunks():
    now = [0.0]
    progress = warm_progress.WarmProgress(["systems"], lambda line: None, show=True, clock=lambda: now[0],
                                          expected={"systems": 100.0})
    progress.plan(0, {"units": 20, "cost": 100.0})
    for _ in range(warm_progress.MIN_CHUNKS - 1):
        now[0] += 10.0
        progress.result(_chunk(seconds=10.0))
    assert progress.overall_left() is None and progress._overall_left_text() == ""


def test_clock_text():
    assert warm_progress.clock_text(65) == "1:05"
    assert warm_progress.clock_text(3912) == "1:05:12"
    assert warm_progress.clock_text(None) == "?"


def test_a_unit_that_reports_part_of_itself_moves_the_stage_and_the_whole_warm_up():
    now = [0.0]
    progress = warm_progress.WarmProgress(["systems", "search"], lambda line: None, show=True, clock=lambda: now[0],
                                          expected={"systems": 100.0, "search": 100.0})
    progress.plan(0, {"units": 4, "cost": 100.0})
    now[0] = 10.0
    for unit in range(4):
        progress.unit_progress(0, unit, 25, 100)        # four big chunks, each a quarter done: one chunk's worth
    assert progress._position() == 1.0
    assert progress._weight_done() == 25.0
    assert progress.overall_left() is None               # one chunk's worth is too little to say how long the rest takes
    for unit in range(4):
        progress.unit_progress(0, unit, 75, 100)        # three chunks' worth
        progress.unit_progress(0, unit, 50, 100)        # a report that goes backwards changes nothing
    now[0] = 30.0
    for unit in range(4):
        progress.unit_progress(0, unit, 100, 100)
    assert progress._weight_done() == 100.0


def test_a_finished_unit_counts_whole_and_its_partial_report_is_dropped():
    now = [0.0]
    progress = warm_progress.WarmProgress(["systems"], lambda line: None, show=True, clock=lambda: now[0],
                                          expected={"systems": 100.0})
    progress.plan(0, {"units": 2, "cost": 100.0})
    progress.unit_progress(0, 0, 90, 100)
    now[0] = 5.0
    progress.result(_chunk(), 0)
    assert progress._position() == 1.0 and progress._weight_done() == 50.0
    progress.unit_progress(1, 1, 50, 100)                # another stage's report is not ours
    progress.unit_progress(0, 1, 50, 0)                  # a unit that has no estimate says nothing
    assert progress._position() == 1.0


def test_the_time_left_uses_the_parts_of_the_chunks_that_have_reported():
    now = [0.0]
    progress = warm_progress.WarmProgress(["systems"], lambda line: None, show=True, clock=lambda: now[0],
                                          expected={"systems": 100.0})
    progress.plan(0, {"units": 20, "cost": 100.0})
    now[0] = 100.0
    for unit in range(10):
        progress.unit_progress(0, unit, 50, 100)         # ten chunks half done: five chunks' worth in 100 s
    assert progress._stage_left(100.0) == 100.0 / 5.0 * 15.0
    assert progress.overall_left() == 100.0 * 75.0 / 25.0
