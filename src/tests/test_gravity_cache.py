# tests/test_gravity_cache.py

"""MAP.387: the gravity map's grid files (`galaxy.gravity_cache`): the round trip, the stamp, torn and foreign files, pruning, and a cache that cannot be written."""

import os
import stat

import numpy as np
import pytest

from planetgen.galaxy import gravity_cache as cache


def _values(zones=4, seed=0):
    return np.random.default_rng(seed).normal(size=(zones ** 3, cache.COLUMNS)).astype(np.float32)


def test_a_grid_reads_back_what_was_written(tmp_path):
    values = _values()
    assert cache.write(str(tmp_path), 7, 4, "0123456789abcdef", values)
    got = cache.read(str(tmp_path), 7, 4, "0123456789abcdef")
    assert got.dtype == np.float32 and got.shape == (64, 10) and got.tobytes() == values.tobytes()
    assert os.listdir(tmp_path) == ["g7_4_0123456789abcdef.bin"]


def test_a_missing_short_torn_or_foreign_file_reads_as_nothing(tmp_path):
    directory = str(tmp_path)
    assert cache.read(directory, 1, 4, "0" * 16) is None
    assert cache.read(os.path.join(directory, "nowhere"), 1, 4, "0" * 16) is None
    cache.write(directory, 1, 4, "0" * 16, _values())
    path = os.path.join(directory, cache.file_name(1, 4, "0" * 16))
    raw = open(path, "rb").read()
    for broken in (raw[:-1], raw[:5], b"", b"XXXX" + raw[4:], raw + b"\0"):
        open(path, "wb").write(broken)
        assert cache.read(directory, 1, 4, "0" * 16) is None
    open(path, "wb").write(raw)
    assert cache.read(directory, 1, 4, "0" * 16) is not None
    assert cache.read(directory, 1, 8, "0" * 16) is None        # another zone count is another file


def test_writing_a_sectors_grid_removes_its_older_files_and_only_its_own(tmp_path):
    directory = str(tmp_path)
    for sector, stamp in ((1, "a" * 16), (1, "b" * 16), (2, "a" * 16)):
        cache.write(directory, sector, 4, stamp, _values())
    assert sorted(os.listdir(directory)) == ["g1_4_" + "b" * 16 + ".bin", "g2_4_" + "a" * 16 + ".bin"]
    (tmp_path / "notes.txt").write_text("keep")
    assert cache.clear(directory) == 2 and os.listdir(directory) == ["notes.txt"]


def test_the_stamp_changes_with_any_part_and_nothing_else():
    base = cache.stamp((1, 16, (3, "2026-11-05 10:00:00.000")))
    assert base == cache.stamp((1, 16, (3, "2026-11-05 10:00:00.000"))) and len(base) == 16
    assert len({base, cache.stamp((2, 16, (3, "2026-11-05 10:00:00.000"))), cache.stamp((1, 16, (3, "2026-11-05 10:00:00.001"))), cache.stamp((1, 8, (3, "x")))}) == 4


def test_a_wrong_shape_is_refused_and_an_unwritable_directory_is_not_an_error(tmp_path):
    with pytest.raises(ValueError):
        cache.write(str(tmp_path), 1, 4, "0" * 16, _values(zones=3))
    if os.geteuid() == 0:
        pytest.skip("root can write anywhere")
    locked = tmp_path / "locked"
    locked.mkdir()
    locked.chmod(stat.S_IRUSR | stat.S_IXUSR)
    try:
        assert cache.write(str(locked), 1, 4, "0" * 16, _values()) is False
    finally:
        locked.chmod(stat.S_IRWXU)
