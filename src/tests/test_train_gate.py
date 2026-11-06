# tests/test_train_gate.py

"""the merge train's gate script splits a file list by time, resumes from its results file and tells a red the train caused from a known one."""

import importlib.util
import pathlib

SCRIPT = pathlib.Path(__file__).resolve().parents[2] / "scripts" / "train_gate.py"
spec = importlib.util.spec_from_file_location("train_gate", SCRIPT)
gate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate)


def test_shards_hold_every_file_once_and_balance_the_time():
    durations = {f"tests/test_{n}.py": float(t) for n, t in enumerate([900, 600, 500, 400, 300, 200, 100, 50, 50, 10])}
    files = sorted(durations) + ["tests/test_unknown.py"]
    lists, totals = gate.shards(files, 3, durations)
    assert sorted(f for lst in lists for f in lst) == sorted(files)
    assert max(totals) - min(totals) <= max(durations.values())


def test_a_results_file_names_what_passed_and_what_failed(tmp_path):
    path = tmp_path / "results.txt"
    path.write_text("PASS tests/test_a.py 12\nFAIL tests/test_b.py 90 rc=1 1 failed\ndone 1/1\n")
    found = gate.read_results([str(path)])
    assert found["tests/test_a.py"] == ("PASS", 12.0, "")
    assert found["tests/test_b.py"][0] == "FAIL" and "1 failed" in found["tests/test_b.py"][2]


def test_a_known_red_is_not_new(tmp_path, monkeypatch):
    path = tmp_path / "results.txt"
    path.write_text("FAIL tests/test_known.py 5 rc=1 x\nFAIL tests/test_new.py 5 rc=1 y\nPASS tests/test_stale.py 5\n")
    known = tmp_path / "known.json"
    known.write_text('{"tests/test_known.py": "red on main", "tests/test_stale.py": "was red"}')
    monkeypatch.setattr(gate, "KNOWN", str(known))
    new, explained, stale = gate.reds([str(path)])
    assert new == ["tests/test_new.py"] and explained == ["tests/test_known.py"] and stale == ["tests/test_stale.py"]
