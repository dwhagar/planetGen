# tests/test_band_table_cache.py

"""
PERF.218: a test run keeps each bright-star sampling table under a key of every file, dataset and tuning constant it reads, so the workers and the
next run do not build it again, a table read back is the table that was made, and a changed input never reads another's table.
"""

import os

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


def test_without_a_folder_the_table_is_built_every_time(monkeypatch):
    monkeypatch.delenv(cache.ENV, raising=False)
    compute = Counter()
    _kept(compute)
    _kept(compute)
    assert compute.calls == 2


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
