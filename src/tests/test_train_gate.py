# tests/test_train_gate.py

"""The merge train's gate script splits a file list by time, resumes from its results file and tells a red the train caused from a known one."""

import importlib.util
import os
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


def test_the_sampling_tests_are_real_files_and_join_when_generation_changes(monkeypatch):
    """TEST.681: tests that count what a few seeded bodies do name no changed module, so the gate lists them by name."""
    tests_dir = pathlib.Path(__file__).resolve().parent
    for name in gate.GUARD_TESTS + gate.SAMPLE_TESTS:
        assert (tests_dir / f"{name}.py").exists(), name
    monkeypatch.setattr(gate, "_git", lambda *args: "src/planetgen/physics/surface_air.py\n")
    assert {f"tests/{name}.py" for name in gate.SAMPLE_TESTS} <= set(gate.gate_files())
    monkeypatch.setattr(gate, "_git", lambda *args: "docs/TODO.md\n")
    assert "tests/test_sheet_fill.py" not in gate.gate_files()


def test_lanes_share_the_files_by_claiming_them_and_a_restarted_lane_takes_back_its_own(tmp_path):
    claims = str(tmp_path / "claims")
    assert gate.claim(claims, "tests/test_a.py", "lane-05")
    assert not gate.claim(claims, "tests/test_a.py", "lane-06")
    assert gate.claim(claims, "tests/test_a.py", "lane-05")
    assert gate.claim(claims, "tests/test_b.py", "lane-06")


def test_the_status_counts_done_failed_running_and_waiting_and_a_dead_lanes_claim_is_released(tmp_path):
    claims = str(tmp_path / "claims")
    results = tmp_path / "r.txt"
    results.write_text("PASS tests/test_a.py 5\nFAIL tests/test_b.py 5 rc=1 x\n")
    files = ["tests/test_a.py", "tests/test_b.py", "tests/test_c.py", "tests/test_d.py"]
    gate.claim(claims, "tests/test_c.py", "lane-07")
    found = gate.status(files, claims, [str(results)])
    assert found == {"done": ["tests/test_a.py"], "failed": ["tests/test_b.py"], "running": ["tests/test_c.py"], "waiting": ["tests/test_d.py"]}
    assert gate.release_stale(claims, files, {"tests/test_a.py", "tests/test_b.py"}, 3600) == []
    old = os.path.getmtime(gate._claim_path(claims, "tests/test_c.py")) - 7200
    os.utime(gate._claim_path(claims, "tests/test_c.py"), (old, old))
    assert gate.release_stale(claims, files, {"tests/test_a.py", "tests/test_b.py"}, 3600) == ["tests/test_c.py"]
    assert gate.status(files, claims, [str(results)])["waiting"] == ["tests/test_c.py", "tests/test_d.py"]


def test_two_lanes_running_on_one_claims_folder_run_each_file_once(tmp_path, monkeypatch):
    ran = []

    class Done:
        returncode = 0

    monkeypatch.setattr(gate.subprocess, "run", lambda command, **kw: ran.append(command[3]) or Done())
    monkeypatch.setattr(gate, "DURATIONS", str(tmp_path / "d.json"))
    files = ["tests/test_a.py", "tests/test_b.py", "tests/test_c.py"]
    claims = str(tmp_path / "claims")
    gate.run(files, str(tmp_path / "r1.txt"), claims=claims, owner="one", say=lambda line: None)
    gate.run(files, str(tmp_path / "r2.txt"), claims=claims, owner="two", say=lambda line: None)
    assert sorted(ran) == files
    assert gate.status(files, claims, [str(tmp_path / "r1.txt"), str(tmp_path / "r2.txt")])["done"] == files


def test_learning_with_replace_keeps_the_longest_pass_so_a_skipped_database_run_does_not_win(tmp_path, monkeypatch):
    durations = tmp_path / "d.json"
    durations.write_text('{"tests/test_db.py": 6.0}')
    monkeypatch.setattr(gate, "DURATIONS", str(durations))
    first, second = tmp_path / "a.txt", tmp_path / "b.txt"
    first.write_text("PASS tests/test_db.py 2400\nFAIL tests/test_x.py 9 rc=1 y\n")
    second.write_text("PASS tests/test_db.py 30\nPASS tests/test_y.py 12\n")
    gate.learn([str(first), str(second)], replace=True)
    kept = gate._load(str(durations), {})
    assert kept == {"tests/test_db.py": 2400.0, "tests/test_y.py": 12.0}
