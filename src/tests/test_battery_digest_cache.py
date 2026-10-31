# tests/test_battery_digest_cache.py

"""
PERF.209: a test run keeps the battery digest it works out, under a key of every file it reads, so the workers and the next run do not work it
out again and a changed file is never given an old answer.
"""

import os

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
    monkeypatch.delenv(cache.ENV, raising=False)
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
