"""PERF.274: looking at how full the disk cache is must not cost a long walk before every chunk of nodes."""

import pytest

from planetgen.web.lib import tilecache


@pytest.fixture
def root(tmp_path, monkeypatch):
    """A cache of two folders and 3,000 bytes; the memory of earlier walks is empty."""
    for folder in ("a", "b"):
        (tmp_path / folder).mkdir()
        (tmp_path / folder / "node.bin").write_bytes(b"x" * 1500)
    (tmp_path / "stamp.json").write_text("{}")      # a kept file never counts
    monkeypatch.setattr(tilecache, "_walked", {})
    return str(tmp_path)


class _Clock:
    """A clock the test moves: every walk of the cache takes `cost` of it."""

    def __init__(self, monkeypatch, cost):
        self.now, self.cost, self.walks = 1000.0, cost, 0
        real = tilecache.cache_files

        def walk(root):
            self.walks += 1
            self.now += self.cost
            return real(root)

        monkeypatch.setattr(tilecache.time, "monotonic", lambda: self.now)
        monkeypatch.setattr(tilecache, "cache_files", walk)


def test_a_quick_walk_is_made_again_every_time(root, monkeypatch):
    clock = _Clock(monkeypatch, cost=0.01)
    assert not tilecache.cache_is_full(root, max_bytes=3001)
    assert tilecache.cache_is_full(root, max_bytes=3000)
    assert tilecache.cache_is_full(root, max_bytes=2000)
    assert clock.walks == 3, "a small cache is read as it is each time"


def test_a_long_walk_is_not_made_again_for_ten_times_as_long(root, monkeypatch):
    clock = _Clock(monkeypatch, cost=20.0)
    assert not tilecache.cache_is_full(root, max_bytes=5000)
    clock.now += 100.0
    assert tilecache.cache_is_full(root, max_bytes=2000), "the bytes the last walk found answer for another budget too"
    assert clock.walks == 1
    clock.now += 101.0      # 20 s walk, ten times that is 200 s after it finished
    assert not tilecache.cache_is_full(root, max_bytes=5000)
    assert clock.walks == 2, "after ten times as long it looks again"


def test_a_walk_tells_the_hook_of_each_folder_it_passes(root, monkeypatch):
    seen = []
    monkeypatch.setattr(tilecache, "ON_DIRECTORY", lambda: seen.append(1))
    files = tilecache.cache_files(root)
    assert len(files) == 2 and len(seen) == 3, "the root and its two folders"
