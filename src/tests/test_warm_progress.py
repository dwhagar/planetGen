"""The progress `planetgen warm --progress` shows for a stage cut into chunks (cli/warm_progress.py)."""

from planetgen.cli import warm_progress


def _chunk(warmed=64, built=64, seconds=38.0):
    return {"warmed": warmed, "built": built, "cached": warmed - built, "left_out": 0, "noun": "nodes",
            "seconds": seconds, "line": "unused"}


def _run(chunks, show=True):
    lines, now = [], [0.0]

    def clock():
        return now[0]

    progress = warm_progress.WarmProgress(["lists", "galaxy-deep-nodes"], lines.append, show=show, clock=clock)
    progress.plan(1, {"units": chunks})
    for _ in range(chunks):
        now[0] += 40.0
        progress.result(_chunk())
    return lines


def test_every_chunk_says_its_number_and_the_total():
    lines = _run(4)
    assert lines[0].startswith("galaxy-deep-nodes chunk 1 of 4: this chunk 64 nodes in 38.0 s; so far 64 nodes (64 built, 0 already stored)")
    assert lines[0].endswith("0:40 elapsed, about 2:00 left")
    assert lines[2].startswith("galaxy-deep-nodes chunk 3 of 4")
    assert lines[-2] == ("galaxy-deep-nodes: 256 nodes warmed (256 built, 0 already stored) in 2:40 over 4 chunks of up "
                         "to 64 nodes")
    assert lines[-1] == "[##############################] 2/2 galaxy-deep-nodes"


def test_a_long_stage_leaves_a_bounded_number_of_lines():
    lines = _run(180)
    chunk_lines = [line for line in lines if " chunk " in line]
    assert 10 <= len(chunk_lines) <= warm_progress.LOG_CHUNK_LINES
    assert chunk_lines[0].startswith("galaxy-deep-nodes chunk 9 of 180")
    assert "180 chunks" in lines[-2]


def test_a_stage_of_one_unit_prints_its_own_line_once():
    lines = []
    progress = warm_progress.WarmProgress(["lists", "search"], lines.append, show=True)
    progress.plan(0, {"units": 1})
    progress.result({"line": "lists: 20 pages warmed"})
    assert lines == ["lists: 20 pages warmed", "[###############---------------] 1/2 lists"]


def test_without_progress_only_the_summary_is_left():
    lines = _run(4, show=False)
    assert len(lines) == 1 and "over 4 chunks" in lines[0]


def test_clock_text():
    assert warm_progress.clock_text(65) == "1:05"
    assert warm_progress.clock_text(3912) == "1:05:12"
    assert warm_progress.clock_text(None) == "?"
