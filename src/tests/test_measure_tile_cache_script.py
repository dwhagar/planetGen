"""scripts/measure_tile_cache.py (PERF.233): the report it makes from a sample of the Galaxy Map's nodes."""

import importlib.util
import pathlib

_PATH = pathlib.Path(__file__).resolve().parents[2] / "scripts" / "measure_tile_cache.py"
_SPEC = importlib.util.spec_from_file_location("measure_tile_cache", _PATH)
measure = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(measure)


def test_the_report_projects_each_levels_sample_to_its_node_count():
    counts = {177147: 2, 3: 1000, 1: 10}
    measured = {177147: [(0.2, 20000), (0.4, 30000)], 3: [(0.1, 40000)], 1: [(0.1, 80000), (0.3, 100000)]}
    report = measure.summarize(counts, measured, nodes_parallel=100, seconds_parallel=5.0, workers=4, budget_mb=1024.0)
    by_level = {level["m"]: level for level in report["levels"]}
    assert by_level[3]["total_bytes"] == 40000 * 1000
    assert abs(by_level[3]["total_seconds"] - 100.0) < 1e-9
    assert by_level[177147]["mean_bytes"] == 25000 and by_level[1]["mean_bytes"] == 90000
    assert report["nodes"] == 1012 and report["parallel_nodes_per_second"] == 20.0
    assert report["bytes"] == 25000 * 2 + 40000 * 1000 + 90000 * 10


def test_a_level_with_no_sample_costs_nothing_and_the_text_names_every_level():
    report = measure.summarize({243: 5, 1: 0}, {}, nodes_parallel=10, seconds_parallel=1.0, workers=2, budget_mb=512.0)
    text = measure.render(report)
    assert "243" in text and "512" in text and report["bytes"] == 0.0


def test_the_arguments_default_to_a_sample_of_forty_on_four_workers():
    args = measure._parse(["--database", "planetgen"])
    assert (args.sample, args.workers) == (40, 4)
