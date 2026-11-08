# tests/test_battery_digest_cache.py

"""
PERF.209, PERF.250: the battery digest is kept under a key of every file it reads, so the workers and the next run do not work it out again
and a changed file is never given an old answer. A test run keeps it where `ENV` says; an install keeps it with its kept star tables.
"""

import os
import sys
import time

import pytest

from planetgen.generation import battery_digest_cache as cache
from planetgen.generation.generator_battery import BatteryDigest, CaseDigest


def _digest(tag):
    return BatteryDigest(tag, {"sector-a": CaseDigest(tag + "f", tag + "s", tag + "x")})


class Counter:
    def __init__(self):
        self.calls = 0

    def __call__(self):
        self.calls += 1
        return _digest(f"d{self.calls}")


def test_without_a_folder_the_digest_is_worked_out_every_time(monkeypatch):
    monkeypatch.setenv(cache.ENV, "")
    compute = Counter()
    cache.digest(compute, BatteryDigest, CaseDigest)
    cache.digest(compute, BatteryDigest, CaseDigest)
    assert compute.calls == 2


def test_a_kept_digest_is_read_back_for_the_same_files(monkeypatch, tmp_path):
    monkeypatch.setenv(cache.ENV, str(tmp_path / "kept"))
    compute = Counter()
    first = cache.digest(compute, BatteryDigest, CaseDigest)
    second = cache.digest(compute, BatteryDigest, CaseDigest)
    assert compute.calls == 1 and second == first


def test_a_changed_file_is_not_given_the_old_digest(monkeypatch, tmp_path):
    package = tmp_path / "package"
    package.mkdir()
    (package / "module.py").write_text("x = 1\n")
    monkeypatch.setattr(cache, "_PACKAGE", str(package))
    monkeypatch.setenv(cache.ENV, str(tmp_path / "kept"))
    compute = Counter()
    cache.digest(compute, BatteryDigest, CaseDigest)
    (package / "module.py").write_text("x = 2\n")
    os.utime(package / "module.py", ns=(1, 1))
    second = cache.digest(compute, BatteryDigest, CaseDigest)
    assert compute.calls == 2 and second.overall == "d2"
    (package / "other.py").write_text("")
    cache.digest(compute, BatteryDigest, CaseDigest)
    assert compute.calls == 3


def test_compiled_files_do_not_change_the_key(monkeypatch, tmp_path):
    package = tmp_path / "package"
    (package / "__pycache__").mkdir(parents=True)
    (package / "module.py").write_text("x = 1\n")
    monkeypatch.setattr(cache, "_PACKAGE", str(package))
    before = cache.key()
    (package / "__pycache__" / "module.cpython-313.pyc").write_bytes(b"\0")
    assert cache.key() == before


def test_a_damaged_copy_is_worked_out_again(monkeypatch, tmp_path):
    folder = tmp_path / "kept"
    monkeypatch.setenv(cache.ENV, str(folder))
    compute = Counter()
    cache.digest(compute, BatteryDigest, CaseDigest)
    for path in folder.iterdir():
        path.write_text("{not json")
    cache.digest(compute, BatteryDigest, CaseDigest)
    assert compute.calls == 2


def test_an_install_keeps_it_with_its_star_tables_when_nothing_names_a_folder(monkeypatch, tmp_path):
    """PERF.250: on the server the variable is not set; the digest goes where the install keeps its kept star tables."""
    from planetgen.generation import band_table_cache
    monkeypatch.delenv(cache.ENV, raising=False)
    monkeypatch.setattr(band_table_cache, "root", lambda: str(tmp_path / "band-tables"))
    compute = Counter()
    first = cache.digest(compute, BatteryDigest, CaseDigest)
    second = cache.digest(compute, BatteryDigest, CaseDigest)
    assert compute.calls == 1 and second == first
    assert [path.name for path in (tmp_path / "band-tables").iterdir()] == [f"battery-{cache.key()}.json"]


def test_a_host_without_the_install_folder_works_it_out_every_time(monkeypatch):
    from planetgen.generation import band_table_cache
    monkeypatch.delenv(cache.ENV, raising=False)
    monkeypatch.setattr(band_table_cache, "root", lambda: None)
    compute = Counter()
    cache.digest(compute, BatteryDigest, CaseDigest)
    cache.digest(compute, BatteryDigest, CaseDigest)
    assert compute.calls == 2


def test_a_kept_digest_is_readable_by_the_other_user_who_runs_the_jobs(monkeypatch, tmp_path):
    monkeypatch.setenv(cache.ENV, str(tmp_path / "kept"))
    old_umask = os.umask(0o077)
    try:
        cache.digest(Counter(), BatteryDigest, CaseDigest)
    finally:
        os.umask(old_umask)
    (path,) = (tmp_path / "kept").iterdir()
    assert path.stat().st_mode & 0o044 == 0o044


@pytest.mark.parametrize("what", ["python", "numpy", "scipy", "machine"])
def test_a_changed_library_or_platform_is_not_given_the_old_digest(monkeypatch, what):
    before = cache.key()
    if what == "python":
        monkeypatch.setattr(sys, "version", sys.version + " other")
    elif what == "machine":
        import platform
        monkeypatch.setattr(platform, "machine", lambda: "other-machine")
    else:
        import importlib
        monkeypatch.setattr(importlib.import_module(what), "__version__", "0.0.1")
    assert cache.key() != before


def test_a_changed_data_file_is_not_given_the_old_digest(monkeypatch, tmp_path):
    from planetgen.datasets import fetch
    data = tmp_path / "data"
    data.mkdir()
    (data / "table.dat").write_text("1")
    monkeypatch.setattr(fetch, "cache_dir", lambda: str(data))
    before = cache.key()
    (data / "table.dat").write_text("22")
    assert cache.key() != before


def test_kept_digests_of_older_code_are_removed_after_a_week_and_recent_ones_stay(monkeypatch, tmp_path):
    folder = tmp_path / "kept"
    folder.mkdir()
    old, recent, other = (folder / "battery-old.json"), (folder / "battery-recent.json"), (folder / "tables.json")
    for path in (old, recent, other):
        path.write_text("{}")
    long_ago = time.time() - (cache.STALE_DAYS + 1) * 86400
    os.utime(old, (long_ago, long_ago))
    os.utime(other, (long_ago, long_ago))
    monkeypatch.setenv(cache.ENV, str(folder))
    cache.digest(Counter(), BatteryDigest, CaseDigest)
    assert not old.exists() and recent.exists() and other.exists()
    assert (folder / f"battery-{cache.key()}.json").exists()


def test_a_folder_that_cannot_be_made_is_ignored(monkeypatch, tmp_path):
    blocker = tmp_path / "blocker"
    blocker.write_text("a file, not a folder")
    monkeypatch.setenv(cache.ENV, str(blocker / "below"))
    compute = Counter()
    assert cache.digest(compute, BatteryDigest, CaseDigest).overall == "d1"
