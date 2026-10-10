# tests/test_benchmark.py

"""
`planetgen benchmark` (PERF.31): the report's layout and the worker-count
argument. The run itself (child processes over scratch databases) is
exercised by hand and its numbers live in docs/design.
"""

import argparse

import pytest

from planetgen.cli import benchmark, generate


def _args(**kw):
    values = dict(sectors=2, ring=12, max_ring=12, seed="ab" * 16)
    values.update(kw)
    return argparse.Namespace(**values)


def _run(workers=1):
    return {"workers": workers, "plan_seconds": 10.0, "fill_seconds": 30.0, "total_seconds": 40.0,
            "stages": [{"command": "galaxy", "key": "sectors", "label": "Generate the sectors", "skipped": False,
                        "seconds": 20.0, "metrics": {"sectors": 2}},
                       {"command": "galaxy", "key": "scatter", "label": "Scatter", "skipped": True,
                        "seconds": 0.0, "metrics": {}}],
            "database": {"Questions": 400, "Com_delete": 0}}


def test_workers_parse_to_whole_numbers():
    assert benchmark._workers("1,2, 4") == [1, 2, 4]


@pytest.mark.parametrize("text", ["", "0", "a", "1,-2"])
def test_bad_worker_counts_are_refused(text):
    with pytest.raises(SystemExit):
        benchmark._workers(text)


def test_the_report_lists_each_run_and_its_stage_shares():
    text = benchmark.format_report(_args(), [_run(1), _run(2)], None)
    assert "1 worker(s)" in text and "2 worker(s)" in text
    assert "Generate the sectors" in text and "50.0%" in text and "sectors 2" in text
    assert "Scatter " not in text  # a skipped stage is left out
    assert "Questions 400 (10/s)" in text and "Com_delete" not in text  # zero counters are left out


def test_the_report_lists_the_profile():
    text = benchmark.format_report(_args(), [_run()], (12.5, [(3.25, 7, "db/store.py:10 save_sector")]))
    assert "Profile of one single-worker fill" in text and "db/store.py:10 save_sector" in text


def test_the_command_is_known_to_the_cli():
    assert "benchmark" in generate.READ_ONLY_COMMANDS
    parser = argparse.ArgumentParser()
    benchmark.add_benchmark_arguments(parser)
    got = parser.parse_args([])
    assert got.sectors == benchmark.DEFAULT_SECTORS and got.ring == benchmark.DEFAULT_RING
