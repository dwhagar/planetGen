# tests/test_private_dirs.py

"""
The tile cache's and the Generate jobs' fallback directory in the system
temp directory (`planetgen/web/lib/privatedir.py`, used by `planetgen/web/lib/tilecache.py` and
`web/jobs.py`): created private (0700), and refused when another user
could have planted it -- a symlink, a directory owned by someone else, or
one other users can write to.
"""

from planetgen.util import settings as settings_model
import os
import stat


import pytest  # noqa: E402

from planetgen.web.lib import tilecache  # noqa: E402
from planetgen.web.lib.privatedir import ensure_private_dir  # noqa: E402
from planetgen.web import jobs  # noqa: E402

UNCREATABLE = "/proc/planetgen-test-cannot-create/dir"


def test_a_new_directory_is_created_mode_0700(tmp_path):
    path = tmp_path / "planetgen-jobs"
    assert ensure_private_dir(str(path)) == str(path)
    assert stat.S_IMODE(os.lstat(path).st_mode) == 0o700


def test_an_existing_private_directory_is_reused(tmp_path):
    path = tmp_path / "planetgen-jobs"
    path.mkdir(mode=0o700)
    (path / "keep").write_text("x")
    assert ensure_private_dir(str(path)) == str(path)
    assert (path / "keep").exists()


@pytest.mark.parametrize("mode", [0o777, 0o770, 0o722, 0o1777])
def test_a_directory_others_can_write_is_refused(tmp_path, mode):
    path = tmp_path / "planetgen-jobs"
    path.mkdir()
    os.chmod(path, mode)
    with pytest.raises(OSError, match="writable by other users"):
        ensure_private_dir(str(path))


def test_a_symlink_is_refused_even_to_a_private_directory(tmp_path):
    target = tmp_path / "elsewhere"
    target.mkdir(mode=0o700)
    link = tmp_path / "planetgen-jobs"
    link.symlink_to(target)
    with pytest.raises(OSError, match="not a directory"):
        ensure_private_dir(str(link))


def test_a_file_is_refused(tmp_path):
    path = tmp_path / "planetgen-jobs"
    path.write_text("x")
    with pytest.raises(OSError, match="not a directory"):
        ensure_private_dir(str(path))


@pytest.mark.skipif(os.geteuid() != 0, reason="needs root to give the directory to another user")
def test_a_directory_owned_by_another_user_is_refused(tmp_path):
    path = tmp_path / "planetgen-jobs"
    path.mkdir(mode=0o700)
    os.chown(path, 54321, -1)
    with pytest.raises(OSError, match="owned by uid 54321"):
        ensure_private_dir(str(path))


@pytest.fixture
def fallback_tmp(tmp_path, monkeypatch):
    """The system temp directory is `tmp_path`, and the default
    directories can't be created, so the fallback is what's used."""
    monkeypatch.setattr(jobs.tempfile, "gettempdir", lambda: str(tmp_path))
    monkeypatch.setattr(tilecache.tempfile, "gettempdir", lambda: str(tmp_path))
    monkeypatch.setattr(jobs, "DEFAULT_JOBS_DIR", UNCREATABLE)
    monkeypatch.setattr(tilecache, "DEFAULT_CACHE_DIR", UNCREATABLE)
    monkeypatch.setattr(jobs, "_config", lambda: settings_model.Jobs())
    monkeypatch.setattr(tilecache, "_config", lambda: settings_model.TileCache())
    for var in ("PLANETGEN_JOBS_DIR", "PLANETGEN_TILE_CACHE_DIR", "PLANETGEN_TILE_CACHE_MAX_MB"):
        monkeypatch.delenv(var, raising=False)
    return tmp_path


def test_jobs_fallback_is_created_private(fallback_tmp):
    path = jobs.jobs_dir()
    assert path == str(fallback_tmp / "planetgen-jobs")
    assert stat.S_IMODE(os.lstat(path).st_mode) == 0o700


def test_tile_cache_fallback_is_created_private(fallback_tmp):
    path = tilecache.cache_dir()
    assert path == str(fallback_tmp / "planetgen-tiles")
    assert stat.S_IMODE(os.lstat(path).st_mode) == 0o700


def test_a_planted_jobs_fallback_is_refused(fallback_tmp):
    planted = fallback_tmp / "planetgen-jobs"
    planted.mkdir()
    os.chmod(planted, 0o777)
    with pytest.raises(OSError, match="no writable jobs directory"):
        jobs.jobs_dir()


def test_a_symlinked_tile_cache_fallback_is_refused(fallback_tmp):
    target = fallback_tmp / "attacker"
    target.mkdir(mode=0o700)
    (fallback_tmp / "planetgen-tiles").symlink_to(target)
    assert tilecache.cache_dir() is None


def test_a_configured_directory_is_used_as_is(tmp_path, monkeypatch):
    # Only the shared-temp fallback is checked; a directory the admin
    # configured keeps working as before (create-cache-dir.sh makes it).
    target = tmp_path / "jobs"
    monkeypatch.setenv("PLANETGEN_JOBS_DIR", str(target))
    assert jobs.jobs_dir() == str(target)
