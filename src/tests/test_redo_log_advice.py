"""PERF.222: `scripts/redo_log_advice.py` sizes MySQL's redo log from the measured write rate and prints the line to set;
it never changes anything."""

import importlib.util
import os

ROOT = os.path.realpath(os.path.join(os.path.dirname(__file__), "..", ".."))
_spec = importlib.util.spec_from_file_location("redo_log_advice", os.path.join(ROOT, "scripts", "redo_log_advice.py"))
advice = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(advice)

GB, MB = 1024 ** 3, 1024 ** 2


def test_the_log_holds_thirty_minutes_of_writes_between_one_and_sixteen_gb():
    assert advice.recommend(0) == 1 * GB                                  # an idle server still gets a log that is not 100 MB
    assert advice.recommend(1 * MB) == 2 * GB                             # 1 MB/s for 30 minutes is 1.76 GB, rounded up
    assert advice.recommend(3 * MB) == 6 * GB
    assert advice.recommend(500 * MB) == 16 * GB                          # capped: recovery time and disk


def test_capacity_comes_from_the_new_variable_or_the_old_pair():
    assert advice.capacity_of({"innodb_redo_log_capacity": str(100 * MB)}) == (100 * MB, True)
    assert advice.capacity_of({"innodb_log_file_size": str(48 * MB), "innodb_log_files_in_group": "2"}) == (96 * MB, False)


def test_a_log_that_is_too_small_gets_the_set_global_line():
    lines = advice.advice(2 * MB, 100 * MB, True)
    text = "\n".join(lines)
    assert "innodb_redo_log_capacity = 4G" in text and f"SET GLOBAL innodb_redo_log_capacity = {4 * GB}" in text
    assert "100 MB" in text and "2 MB/s written" in text


def test_a_log_that_is_big_enough_changes_nothing_unless_the_checkpointer_lags():
    assert "nothing to change" in "\n".join(advice.advice(1 * MB, 4 * GB, True))
    lagging = "\n".join(advice.advice(1 * MB, 4 * GB, True, warnings=84))
    assert "84 checkpointer warning" in lagging and "innodb_redo_log_capacity = 8G" in lagging     # doubled


def test_an_old_server_is_told_to_restart_and_never_to_set_global():
    text = "\n".join(advice.advice(2 * MB, 96 * MB, False))
    assert "innodb_log_file_size = 2G" in text and "restart MySQL" in text and "SET GLOBAL" not in text


def test_the_error_log_warnings_are_counted_or_unread(tmp_path):
    log = tmp_path / "error.log"
    log.write_text("a\n[Warning] [MY-014084] [InnoDB] log_checkpointer lagging\nb\n[Warning] [MY-014084] again\n")
    assert advice.count_warnings(str(log)) == 2
    assert advice.count_warnings(str(tmp_path / "missing.log")) is None


def test_a_server_that_cannot_be_asked_is_skipped_with_a_message(monkeypatch, capsys):
    from planetgen.db import store
    monkeypatch.setattr(store, "get_connection", lambda *a, **k: (_ for _ in ()).throw(OSError("no server")))
    assert advice.main([]) == 0
    assert "not checked" in capsys.readouterr().err
