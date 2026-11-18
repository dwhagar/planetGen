# tests/test_band_table_cache.py

"""
PERF.218: a test run keeps each bright-star sampling table under a key of every file, dataset and tuning constant it reads, so the workers and the
next run do not build it again, a table read back is the table that was made, and a changed input never reads another's table.
"""

import os
import sys
import threading
import time

from planetgen import tuning
from planetgen.generation import band_table_cache as cache
from planetgen.generation import star_population


class Counter:
    def __init__(self):
        self.calls = 0

    def __call__(self):
        self.calls += 1
        return {"cells": [(1.0, 2.0, 2.3, 0.5)], "n": self.calls}


def _kept(compute, arguments=("a", 1.5)):
    return cache.kept("test", arguments, compute, lambda table: table["cells"], lambda form: {"cells": [tuple(c) for c in form], "n": 0})


def test_without_a_folder_the_table_is_built_every_time(monkeypatch, tmp_path):
    monkeypatch.delenv(cache.ENV, raising=False)
    monkeypatch.setattr(cache, "DEFAULT_ROOT", str(tmp_path / "not-made"))
    compute = Counter()
    _kept(compute)
    _kept(compute)
    assert compute.calls == 2
    assert not (tmp_path / "not-made").exists(), "the program never makes the install's folder, update.sh does (PERF.188)"


def test_an_install_keeps_the_tables_in_its_folder_and_only_the_current_code_keys_tables(monkeypatch, tmp_path):
    root = tmp_path / "band-tables"
    root.mkdir()
    older = root / "0123456789abcdef0123"
    older.mkdir()
    (older / "bright-old.json").write_text("[]")
    monkeypatch.delenv(cache.ENV, raising=False)
    monkeypatch.setattr(cache, "DEFAULT_ROOT", str(root))
    compute = Counter()
    first = _kept(compute)
    again = _kept(compute)
    assert compute.calls == 1 and again["cells"] == first["cells"]
    assert [path.name for path in root.iterdir()] == [cache._code_key()[:20]], "the tables of an older version of the code can never be read again"


def test_a_variable_set_but_empty_keeps_nothing_even_on_an_install(monkeypatch, tmp_path):
    root = tmp_path / "band-tables"
    root.mkdir()
    monkeypatch.setenv(cache.ENV, "")
    monkeypatch.setattr(cache, "DEFAULT_ROOT", str(root))
    compute = Counter()
    _kept(compute)
    _kept(compute)
    assert compute.calls == 2 and not list(root.iterdir())


def test_a_kept_table_is_read_back_for_the_same_arguments_and_not_for_others(monkeypatch, tmp_path):
    monkeypatch.setenv(cache.ENV, str(tmp_path / "kept"))
    compute = Counter()
    first = _kept(compute)
    again = _kept(compute)
    other = _kept(compute, ("a", 2.5))
    assert compute.calls == 2 and again["cells"] == first["cells"] and again["n"] == 0 and other["n"] == 2


def test_a_changed_tuning_constant_reads_no_old_table(monkeypatch, tmp_path):
    monkeypatch.setenv(cache.ENV, str(tmp_path / "kept"))
    compute = Counter()
    _kept(compute)
    monkeypatch.setattr(tuning, "BRIGHT_TABLE_BOUND_MARGIN", tuning.BRIGHT_TABLE_BOUND_MARGIN * 1.01)
    _kept(compute)
    assert compute.calls == 2


def test_an_unreadable_kept_file_is_built_again(monkeypatch, tmp_path):
    monkeypatch.setenv(cache.ENV, str(tmp_path / "kept"))
    compute = Counter()
    _kept(compute)
    for folder, _names, files in os.walk(tmp_path):
        for name in files:
            with open(os.path.join(folder, name), "w", encoding="utf-8") as handle:
                handle.write("not json")
    _kept(compute)
    assert compute.calls == 2


def test_a_folder_that_cannot_be_made_builds_the_table(monkeypatch, tmp_path):
    blocker = tmp_path / "file"
    blocker.write_text("x")
    monkeypatch.setenv(cache.ENV, str(blocker / "below"))
    compute = Counter()
    assert _kept(compute)["n"] == 1


def test_a_sampling_table_read_back_is_the_table_that_was_made(monkeypatch, tmp_path):
    monkeypatch.setenv(cache.ENV, str(tmp_path / "kept"))
    population = star_population.POPULATIONS[0]
    args = (1000.0, population, None, None, -0.5)
    star_population._bright_table.cache_clear()
    made = star_population._bright_table(*args)
    star_population._bright_table.cache_clear()
    read = star_population._bright_table(*args)
    star_population._bright_table.cache_clear()
    assert read == made and read is not made
    assert any(name.startswith("bright-") for _folder, _names, files in os.walk(tmp_path) for name in files)
    band = (1.0, None)
    star_population._mass_band_table.cache_clear()
    made = star_population._mass_band_table(band, population, None, -0.5, True)
    star_population._mass_band_table.cache_clear()
    assert star_population._mass_band_table(band, population, None, -0.5, True) == made
    star_population._mass_band_table.cache_clear()


# --- PERF.281: the key follows what a table reads, and one process builds a table ---------------------------------------------------------------


def _tracing(build):
    """The files of the `planetgen` package that `build()` runs a function of."""
    ran = set()

    def profile(frame, event, _arg):
        if event == "call":
            ran.add(frame.f_code.co_filename)

    sys.setprofile(profile)
    try:
        build()
    finally:
        sys.setprofile(None)
    return {os.path.relpath(name, cache._PACKAGE) for name in ran if name.startswith(cache._PACKAGE + os.sep)}


def test_a_table_build_runs_only_files_the_key_covers(monkeypatch):
    """The key reads `INPUT_FILES`, so a build that ran a file outside it would keep a stale table after that file changed: every population, both kinds of table, a bin of each composition."""
    monkeypatch.setenv(cache.ENV, "")                                  # build, keep nothing
    ran = set()
    for population in star_population.SCATTERED_POPULATIONS if hasattr(star_population, "SCATTERED_POPULATIONS") else star_population.POPULATIONS:
        for fe_bin in (-1.0, 0.3):
            star_population._bright_table.cache_clear()
            star_population._mass_band_table.cache_clear()
            ran |= _tracing(lambda: (star_population._bright_table(5000.0, population, None, None, fe_bin),
                                     star_population._bright_table(60.0, population, (1.0, 2.0), None, fe_bin),
                                     star_population._mass_band_table((1.0, 2.0), population, None, fe_bin, True),
                                     star_population._mass_band_table((5.0, None), population, None, fe_bin)))
    star_population._bright_table.cache_clear()
    star_population._mass_band_table.cache_clear()
    assert ran - set(cache.INPUT_FILES) == set(), "add these to band_table_cache.INPUT_FILES (the key must follow every file a table build runs)"
    assert "physics/post_main_sequence.py" in ran and "generation/star_population.py" in ran


def test_the_key_follows_the_contents_of_the_inputs_and_not_the_rest_of_the_package(monkeypatch, tmp_path):
    (tmp_path / "read.py").write_text("x = 1\n")
    (tmp_path / "other.py").write_text("y = 1\n")
    monkeypatch.setattr(cache, "_PACKAGE", str(tmp_path))
    monkeypatch.setattr(cache, "INPUT_FILES", ("read.py",))
    cache._code_key.cache_clear()
    cache._tuning_names.cache_clear()
    try:
        first = cache._code_key()
        (tmp_path / "other.py").write_text("y = 2\n")
        (tmp_path / "read.py").touch()
        cache._code_key.cache_clear()
        assert cache._code_key() == first, "a file the tables do not read, or the same code with new times, keeps the key"
        (tmp_path / "read.py").write_text("x = 2\n")
        cache._code_key.cache_clear()
        assert cache._code_key() != first
    finally:
        cache._code_key.cache_clear()
        cache._tuning_names.cache_clear()


def test_the_tuning_constants_in_the_key_are_the_ones_the_table_files_read():
    names = cache._tuning_names()
    assert "BRIGHT_TABLE_BOUND_MARGIN" in names and "IMF_SLOPES" in names
    assert "ASTEROID_FIELD_RADIUS_RANGE_LY" not in names and len(names) < len([n for n in vars(tuning) if n.isupper()]) / 2


def test_two_processes_that_ask_for_one_table_build_it_once(monkeypatch, tmp_path):
    monkeypatch.setenv(cache.ENV, str(tmp_path / "kept"))
    calls, found = [], []

    def slow():
        calls.append(1)
        time.sleep(0.5)
        return {"cells": [(1.0, 2.0, 2.3, 0.5)], "n": 1}

    threads = [threading.Thread(target=lambda: found.append(_kept(slow))) for _ in range(3)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    assert len(calls) == 1 and len(found) == 3 and all(table["cells"] == [(1.0, 2.0, 2.3, 0.5)] for table in found)
