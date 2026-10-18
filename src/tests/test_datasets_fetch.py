# tests/test_datasets_fetch.py

"""The shared setup-time download step (GEN.548): sha256 check, atomic keep, stand-in, User-Agent, never at run time."""

import hashlib
import io
import os
import pathlib
import re
import urllib.error

import pytest

from planetgen.cli import datasets as datasets_cli
from planetgen.datasets import catalogue, fetch
from planetgen.util import user_agent as ua
from planetgen._version import REPO_URL, __version__
from planetgen.api.version import API_VERSION

SRC = pathlib.Path(__file__).resolve().parents[1] / "planetgen"


def _piece(name, payload, url=None):
    return fetch.Piece(name, url or f"https://example.org/{name}", hashlib.sha256(payload).hexdigest(), len(payload))


def _dataset(*payloads, key="fixture"):
    pieces = tuple(_piece(f"part{i}.dat", payload, f"https://example.org/{key}/part{i}.dat") for i, payload in enumerate(payloads))
    return fetch.Dataset(key, "A fixture", "Nobody 2026", "open", pieces, "a flat table", "XX1")


class _Response(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


class _Host:
    """An opener that serves bytes the test built, remembers the requests and can fail."""

    def __init__(self, files, error=None):
        self.files, self.error, self.requests = files, error, []

    def __call__(self, request, timeout=None):
        self.requests.append((request.full_url, dict(request.header_items()), timeout))
        if self.error:
            raise self.error
        return _Response(self.files[request.full_url])


@pytest.fixture(autouse=True)
def cache(tmp_path, monkeypatch):
    monkeypatch.setenv(fetch.CACHE_DIR_ENV, str(tmp_path / "cache"))
    fetch._VERIFIED.clear()
    return tmp_path / "cache"


def _serve(dataset, *payloads):
    return _Host({piece.url: payload for piece, payload in zip(dataset.pieces, payloads)})


def test_a_fetch_downloads_checks_and_keeps_every_file(cache):
    payloads = (b"alpha" * 100, b"beta" * 50)
    dataset = _dataset(*payloads)
    result = fetch.fetch(dataset, opener=_serve(dataset, *payloads))
    assert (result.status, result.reasons) == ("fetched", ())
    for piece, payload in zip(dataset.pieces, payloads):
        assert (cache / "fixture" / piece.name).read_bytes() == payload
    assert not [name for name in os.listdir(cache / "fixture") if name.startswith(".")]


def test_a_second_fetch_asks_for_nothing(cache):
    dataset = _dataset(b"x" * 10)
    fetch.fetch(dataset, opener=_serve(dataset, b"x" * 10))
    host = _serve(dataset, b"x" * 10)
    assert fetch.fetch(dataset, opener=host).status == "present"
    assert host.requests == []


def test_every_request_sends_the_standard_user_agent(cache):
    dataset = _dataset(b"abc")
    host = _serve(dataset, b"abc")
    fetch.fetch(dataset, opener=host)
    (url, headers, timeout), = host.requests
    assert headers["User-agent"] == ua.user_agent()
    assert timeout == fetch.TIMEOUT_SECONDS


def test_the_user_agent_is_product_tokens_with_the_release_url_and_api_version():
    text = ua.user_agent()
    assert text == f"PlanetGen/{__version__} (+{REPO_URL}) PlanetGen-API/{API_VERSION}"
    assert re.fullmatch(r"PlanetGen/\d+\.\d+\.\d+ \(\+https://\S+\) PlanetGen-API/\d+", text)


def test_a_wrong_hash_is_not_kept_and_names_both_hashes(cache):
    dataset = _dataset(b"good bytes")
    result = fetch.fetch(dataset, opener=_serve(dataset, b"evil bytes"))
    assert result.status == "stand-in"
    assert "sha256" in result.reasons[0] and dataset.pieces[0].sha256 in result.reasons[0]
    assert not (cache / "fixture" / "part0.dat").exists()


def test_a_file_longer_or_shorter_than_the_catalogue_says_is_not_kept(cache):
    dataset = _dataset(b"12345678")
    for wrong in (b"1234", b"12345678" + b"9" * 100000):
        result = fetch.fetch(dataset, opener=_serve(dataset, wrong))
        assert result.status == "stand-in" and "bytes" in result.reasons[0]
    assert not (cache / "fixture" / "part0.dat").exists()


@pytest.mark.parametrize("error", [
    urllib.error.URLError("name resolution failed"),
    urllib.error.HTTPError("https://example.org/fixture/part0.dat", 503, "Unavailable", {}, None),
    TimeoutError("timed out"),
    ConnectionResetError("reset"),
])
def test_a_host_that_is_down_gives_a_stand_in_and_a_reason_never_an_exception(cache, error):
    dataset = _dataset(b"abc")
    result = fetch.fetch(dataset, opener=_Host({}, error))
    assert result.status == "stand-in" and "example.org/fixture/part0.dat" in result.reasons[0]


def test_a_failure_keeps_the_good_files_already_there(cache):
    dataset = _dataset(b"first", b"second")
    fetch.fetch(dataset, opener=_serve(dataset, b"first", b"second"))
    (cache / "fixture" / "part1.dat").write_bytes(b"damaged")
    host = _Host({dataset.pieces[0].url: b"first"}, urllib.error.URLError("down"))
    result = fetch.fetch(dataset, opener=host)
    assert result.status == "stand-in" and len(host.requests) == 1
    assert (cache / "fixture" / "part0.dat").read_bytes() == b"first"


def test_a_changed_source_replaces_the_old_file(cache):
    old = _dataset(b"old values")
    fetch.fetch(old, opener=_serve(old, b"old values"))
    new = _dataset(b"new values!")
    assert fetch.fetch(new, opener=_serve(new, b"new values!")).status == "fetched"
    assert (cache / "fixture" / "part0.dat").read_bytes() == b"new values!"


def test_locate_gives_the_folder_when_every_file_is_there(cache):
    dataset = _dataset(b"abc", b"defg")
    fetch.fetch(dataset, opener=_serve(dataset, b"abc", b"defg"))
    located = fetch.locate(dataset)
    assert not located.standin and located.label.startswith("downloaded: ")
    assert pathlib.Path(located.path("part1.dat")).read_bytes() == b"defg"
    with pytest.raises(KeyError):
        located.path("other.dat")


def test_locate_on_an_empty_cache_is_the_labelled_stand_in(cache):
    dataset = _dataset(b"abc")
    located = fetch.locate(dataset)
    assert located.standin and "has not been downloaded" in located.reason
    assert located.label == "stand-in (XX1): a flat table"
    with pytest.raises(FileNotFoundError):
        located.path("part0.dat")


def test_locate_refuses_a_file_that_was_altered_after_it_was_kept(cache):
    dataset = _dataset(b"abcd")
    fetch.fetch(dataset, opener=_serve(dataset, b"abcd"))
    (cache / "fixture" / "part0.dat").write_bytes(b"abce")
    assert fetch.locate(dataset).standin
    (cache / "fixture" / "part0.dat").write_bytes(b"ab")
    assert "size or sha256" in fetch.locate(dataset).reason


def test_locate_never_opens_a_connection(cache, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("a run-time read opened a connection")
    monkeypatch.setattr("urllib.request.urlopen", forbidden)
    monkeypatch.setattr("socket.socket.connect", forbidden)
    assert fetch.locate(_dataset(b"abc")).standin


def test_locate_hashes_a_file_once_per_process(cache, monkeypatch):
    dataset = _dataset(b"abcd")
    fetch.fetch(dataset, opener=_serve(dataset, b"abcd"))
    fetch._VERIFIED.clear()
    calls = []
    real = fetch._sha256_of
    monkeypatch.setattr(fetch, "_sha256_of", lambda path: calls.append(path) or real(path))
    fetch.locate(dataset)
    fetch.locate(dataset)
    assert len(calls) == 1


@pytest.mark.parametrize("make", [
    lambda: fetch.Piece("../x", "https://e.org/x", "0" * 64, 1),
    lambda: fetch.Piece("x", "http://e.org/x", "0" * 64, 1),
    lambda: fetch.Piece("x", "https://e.org/x", "ABC", 1),
    lambda: fetch.Piece("x", "https://e.org/x", "0" * 64, 0),
    lambda: fetch.Dataset("a/b", "t", "s", "l", (_piece("x", b"1"),), "s", "f"),
    lambda: fetch.Dataset("a", "t", "s", "l", (), "s", "f"),
    lambda: fetch.Dataset("a", "t", "s", "l", (_piece("x", b"1"), _piece("x", b"2")), "s", "f"),
])
def test_a_malformed_catalogue_entry_is_refused(make):
    with pytest.raises(ValueError):
        make()


def test_every_catalogue_entry_is_well_formed_and_keys_are_unique():
    keys = [dataset.key for dataset in catalogue.ALL]
    assert len(keys) == len(set(keys))
    for dataset in catalogue.ALL:
        assert dataset.source and dataset.licence and dataset.standin and dataset.fill


def test_the_cache_folder_defaults_outside_the_checkout_and_can_move(monkeypatch):
    monkeypatch.delenv(fetch.CACHE_DIR_ENV)
    assert fetch.cache_dir() == "/var/cache/planetgen/datasets"
    monkeypatch.setenv(fetch.CACHE_DIR_ENV, "/elsewhere")
    assert fetch.cache_dir() == "/elsewhere"


def test_the_command_fetches_reports_and_exits_1_on_a_stand_in(cache, monkeypatch, capsys):
    good, bad = _dataset(b"abc", key="good"), _dataset(b"xyz", key="bad")
    monkeypatch.setattr(catalogue, "ALL", (good, bad))
    monkeypatch.setattr(fetch.urllib.request, "urlopen", _Host(
        {good.pieces[0].url: b"abc", bad.pieces[0].url: b"wrong"}))
    assert datasets_cli.main(["fetch"]) == 1
    out = capsys.readouterr().out
    assert "good: fetched" in out and "bad: stand-in (XX1)" in out and "uses: a flat table" in out
    assert datasets_cli.main(["status", "good"]) == 0
    assert datasets_cli.main(["status", "bad"]) == 1


def test_the_command_with_nothing_registered_is_a_clean_no_op(monkeypatch, capsys):
    monkeypatch.setattr(catalogue, "ALL", ())
    assert datasets_cli.main(["fetch"]) == 0
    assert "No datasets" in capsys.readouterr().out


def test_the_command_refuses_an_unknown_dataset(monkeypatch):
    monkeypatch.setattr(catalogue, "ALL", (_dataset(b"a"),))
    with pytest.raises(SystemExit) as exit_info:
        datasets_cli.main(["fetch", "nope"])
    assert "unknown dataset nope" in str(exit_info.value)


def test_nothing_that_runs_with_the_site_downloads():
    """Downloads happen at setup only (Boss, 2026-10-18): the generator's packages never open a URL themselves."""
    pattern = re.compile(r"^\s*(import|from)\s+(urllib\.request|http\.client|requests|httpx|aiohttp)\b", re.M)
    offenders = [str(path.relative_to(SRC)) for folder in ("physics", "generation", "galaxy", "population", "describe", "names")
                 for path in (SRC / folder).rglob("*.py") if pattern.search(path.read_text())]
    assert not offenders, offenders


def test_the_scripts_run_the_download_step_at_setup():
    root = SRC.parents[1]
    common, update, install = [(root / name).read_text() for name in ("scripts/deploy-common.sh", "update.sh", "install.sh")]
    assert "run_cli datasets fetch" in common
    assert "ensure_datasets" in update and "ensure_datasets" in install
    assert update.index("ensure_nltk_words") < update.index("ensure_datasets") < update.index("== 4/9")
