"""scripts/band_table_share.py (PERF.232): the shared copy of the sampling tables is read only by the same code."""

import importlib.util
import pathlib

import pytest

from planetgen.generation import band_table_cache

_PATH = pathlib.Path(__file__).resolve().parents[2] / "scripts" / "band_table_share.py"
_SPEC = importlib.util.spec_from_file_location("band_table_share", _PATH)
share = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(share)


@pytest.fixture
def folders(tmp_path, monkeypatch):
    local = tmp_path / "local"
    monkeypatch.setenv(band_table_cache.ENV, str(local))
    monkeypatch.setenv(share.SHARE_ENV, str(tmp_path / "shared"))
    return local, tmp_path / "shared"


def _make_table(local, name):
    folder = pathlib.Path(share._local_folder())
    (folder / name).write_text("[1.0]")
    return folder


def test_the_key_follows_the_contents_of_the_code_and_not_the_file_times(tmp_path):
    package = tmp_path / "pkg"
    package.mkdir()
    (package / "a.py").write_text("x = 1")
    datasets = tmp_path / "data"
    datasets.mkdir()
    first = share.content_key(str(package), str(datasets))
    (package / "a.py").touch()
    assert share.content_key(str(package), str(datasets)) == first
    (package / "a.py").write_text("x = 2")
    assert share.content_key(str(package), str(datasets)) != first


def test_saved_tables_come_back_into_a_fresh_local_folder_and_only_for_the_same_code(folders, monkeypatch):
    local, shared = folders
    _make_table(local, "bright-aaa.json")
    assert share.save() == 1 and share.save() == 0
    for path in local.iterdir():
        for item in path.iterdir():
            item.unlink()
    assert share.restore() == 1
    assert [item.name for path in local.iterdir() for item in path.iterdir()] == ["bright-aaa.json"]
    monkeypatch.setattr(share, "content_key", lambda *args: "otherkey")
    assert share.restore() == 0
