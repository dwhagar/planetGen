# tests/test_datasets_fetch.py

"""The shared setup-time download step (GEN.548): sha256 check, atomic keep, missing, User-Agent, never at run time."""

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
    monkeypatch.setattr(fetch.time, "sleep", lambda seconds: None)
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
    assert result.status == "missing"
    assert "sha256" in result.reasons[0] and dataset.pieces[0].sha256 in result.reasons[0]
    assert not (cache / "fixture" / "part0.dat").exists()


def test_a_file_longer_or_shorter_than_the_catalogue_says_is_not_kept(cache):
    dataset = _dataset(b"12345678")
    for wrong in (b"1234", b"12345678" + b"9" * 100000):
        result = fetch.fetch(dataset, opener=_serve(dataset, wrong))
        assert result.status == "missing" and "bytes" in result.reasons[0]
    assert not (cache / "fixture" / "part0.dat").exists()


@pytest.mark.parametrize("error", [
    urllib.error.URLError("name resolution failed"),
    urllib.error.HTTPError("https://example.org/fixture/part0.dat", 503, "Unavailable", {}, None),
    TimeoutError("timed out"),
    ConnectionResetError("reset"),
])
def test_a_host_that_is_down_gives_missing_and_a_reason_never_an_exception(cache, error):
    dataset = _dataset(b"abc")
    result = fetch.fetch(dataset, opener=_Host({}, error))
    assert result.status == "missing" and "example.org/fixture/part0.dat" in result.reasons[0]


def test_a_failure_keeps_the_good_files_already_there(cache):
    dataset = _dataset(b"first", b"second")
    fetch.fetch(dataset, opener=_serve(dataset, b"first", b"second"))
    (cache / "fixture" / "part1.dat").write_bytes(b"damaged")
    host = _Host({dataset.pieces[0].url: b"first"}, urllib.error.URLError("down"))
    result = fetch.fetch(dataset, opener=host)
    assert result.status == "missing" and len(host.requests) == fetch.ATTEMPTS
    assert (cache / "fixture" / "part0.dat").read_bytes() == b"first"


class _Flaky(_Host):
    """Fails the first `failures` requests, then serves."""

    def __init__(self, files, failures):
        super().__init__(files)
        self.failures = failures

    def __call__(self, request, timeout=None):
        if len(self.requests) < self.failures:
            self.requests.append((request.full_url, {}, timeout))
            raise urllib.error.URLError("reset")
        return super().__call__(request, timeout)


def test_a_host_that_answers_on_a_second_attempt_is_fine(cache):
    dataset = _dataset(b"abcdef")
    host = _Flaky({dataset.pieces[0].url: b"abcdef"}, failures=2)
    assert fetch.fetch(dataset, opener=host).status == "fetched" and len(host.requests) == 3


def test_a_wrong_answer_is_not_retried(cache):
    dataset = _dataset(b"abcdef")
    host = _serve(dataset, b"zzzzzz")
    assert fetch.fetch(dataset, opener=host).status == "missing" and len(host.requests) == 1


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
    assert not located.missing and located.label.startswith("downloaded: ")
    assert pathlib.Path(located.path("part1.dat")).read_bytes() == b"defg"
    with pytest.raises(KeyError):
        located.path("other.dat")


def test_locate_on_an_empty_cache_is_missing_with_its_label(cache):
    dataset = _dataset(b"abc")
    located = fetch.locate(dataset)
    assert located.missing and "has not been downloaded" in located.reason
    assert located.label == "missing (XX1, replaces a flat table)"
    with pytest.raises(FileNotFoundError):
        located.path("part0.dat")


def test_locate_refuses_a_file_that_was_altered_after_it_was_kept(cache):
    dataset = _dataset(b"abcd")
    fetch.fetch(dataset, opener=_serve(dataset, b"abcd"))
    (cache / "fixture" / "part0.dat").write_bytes(b"abce")
    assert fetch.locate(dataset).missing
    (cache / "fixture" / "part0.dat").write_bytes(b"ab")
    assert "size or sha256" in fetch.locate(dataset).reason


def test_locate_never_opens_a_connection(cache, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("a run-time read opened a connection")
    monkeypatch.setattr("urllib.request.urlopen", forbidden)
    monkeypatch.setattr("socket.socket.connect", forbidden)
    assert fetch.locate(_dataset(b"abc")).missing


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
        assert dataset.source and dataset.licence and dataset.replaces and dataset.fill


def test_the_cache_folder_defaults_outside_the_checkout_and_can_move(monkeypatch):
    monkeypatch.delenv(fetch.CACHE_DIR_ENV)
    assert fetch.cache_dir() == "/var/cache/planetgen/datasets"
    monkeypatch.setenv(fetch.CACHE_DIR_ENV, "/elsewhere")
    assert fetch.cache_dir() == "/elsewhere"


def test_the_command_fetches_reports_and_exits_1_on_a_missing_dataset(cache, monkeypatch, capsys):
    good, bad = _dataset(b"abc", key="good"), _dataset(b"xyz", key="bad")
    monkeypatch.setattr(catalogue, "ALL", (good, bad))
    monkeypatch.setattr(fetch.urllib.request, "urlopen", _Host(
        {good.pieces[0].url: b"abc", bad.pieces[0].url: b"wrong"}))
    assert datasets_cli.main(["fetch"]) == 1
    out = capsys.readouterr().out
    assert "good: fetched" in out and "bad: missing (XX1)" in out and "generation will not start" in out and "replaces: a flat table" in out
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


# ---- datasets that keep compact files derived from a large archive ----

def _archive_dataset(derive=None, version=1, payload=b"1,2\n3,4\n"):
    piece = fetch.Piece("big.csv", "https://example.org/big.csv", hashlib.sha256(payload).hexdigest(), len(payload),
                        transient=True)

    def sums(folder):
        text = (pathlib.Path(folder) / "big.csv").read_text()
        return {"sum.txt": str(sum(int(v) for line in text.split() for v in line.split(","))).encode()}

    return fetch.Dataset("archive", "A fixture", "Nobody 2026", "open", (piece,), "a flat table", "XX2",
                         derived=("sum.txt",), derive=derive or sums, derive_version=version), payload


def test_an_archive_is_derived_once_and_then_removed(cache):
    dataset, payload = _archive_dataset()
    host = _Host({dataset.pieces[0].url: payload})
    assert fetch.fetch(dataset, opener=host).status == "fetched"
    assert (cache / "archive" / "sum.txt").read_text() == "10"
    assert not (cache / "archive" / "big.csv").exists()
    located = fetch.locate(dataset)
    assert not located.missing and pathlib.Path(located.path("sum.txt")).read_text() == "10"
    again = _Host({dataset.pieces[0].url: payload})
    assert fetch.fetch(dataset, opener=again).status == "present" and again.requests == []


def test_a_new_derive_version_or_catalogue_hash_derives_again(cache):
    dataset, payload = _archive_dataset()
    fetch.fetch(dataset, opener=_Host({dataset.pieces[0].url: payload}))
    newer, _ = _archive_dataset(version=2)
    assert fetch.locate(newer).missing
    assert fetch.fetch(newer, opener=_Host({newer.pieces[0].url: payload})).status == "fetched"
    assert not fetch.locate(newer).missing
    other, other_payload = _archive_dataset(payload=b"5,6\n")
    assert fetch.locate(other).missing


def test_a_derived_file_altered_afterwards_is_missing(cache):
    dataset, payload = _archive_dataset()
    fetch.fetch(dataset, opener=_Host({dataset.pieces[0].url: payload}))
    (cache / "archive" / "sum.txt").write_text("11")
    assert fetch.locate(dataset).missing


def test_a_derive_function_that_fails_is_missing_not_an_exception(cache):
    def broken(folder):
        raise ValueError("a column moved")
    dataset, payload = _archive_dataset(derive=broken)
    result = fetch.fetch(dataset, opener=_Host({dataset.pieces[0].url: payload}))
    assert result.status == "missing" and "a column moved" in result.reasons[0]
    assert fetch.locate(dataset).missing and not (cache / "archive" / ".derived.json").exists()


def test_derive_must_make_exactly_the_files_the_catalogue_names(cache):
    dataset, payload = _archive_dataset(derive=lambda folder: {"other.txt": b"x"})
    assert fetch.fetch(dataset, opener=_Host({dataset.pieces[0].url: payload})).status == "missing"


def test_a_transient_file_needs_a_derive_function():
    with pytest.raises(ValueError):
        fetch.Dataset("a", "t", "s", "l", (fetch.Piece("x", "https://e.org/x", "0" * 64, 1, transient=True),), "s", "f")


# ---- the URL of every file is printed (GEN.600) ----

def test_fetch_prints_the_url_of_a_fetched_a_kept_and_a_failed_file(cache, monkeypatch, capsys):
    first, second = _dataset(b"abc", b"def", key="mixed"), _dataset(b"xyz", key="bad")
    monkeypatch.setattr(catalogue, "ALL", (first, second))
    kept = cache / "mixed"
    kept.mkdir(parents=True)
    (kept / "part0.dat").write_bytes(b"abc")                       # already there with the right hash
    host = _Host({first.pieces[1].url: b"def", second.pieces[0].url: b"zzz"})
    monkeypatch.setattr(fetch.urllib.request, "urlopen", host)
    assert datasets_cli.main(["fetch"]) == 1
    out = capsys.readouterr().out
    assert f"kept     {first.pieces[0].url}" in out
    assert f"fetched  {first.pieces[1].url}" in out
    assert f"failed   {second.pieces[0].url} (" in out and "sha256" in out
    assert out.count(second.pieces[0].url) >= 2                    # in the dataset's summary and on its own line
    assert [url for url, _headers, _timeout in host.requests] == [first.pieces[1].url, second.pieces[0].url]


def test_fetch_gives_every_piece_a_result(cache):
    dataset = _dataset(b"abc", b"def")
    result = fetch.fetch(dataset, opener=_serve(dataset, b"abc", b"def"))
    assert [(r.piece.name, r.state) for r in result.pieces] == [("part0.dat", "fetched"), ("part1.dat", "fetched")]
    again = fetch.fetch(dataset, opener=_serve(dataset, b"abc", b"def"))
    assert [r.state for r in again.pieces] == ["kept", "kept"]


def test_status_prints_the_url_of_every_missing_file_and_none_of_a_present_one(cache, monkeypatch, capsys):
    dataset = _dataset(b"abc", b"def", key="half")
    monkeypatch.setattr(catalogue, "ALL", (dataset,))
    (cache / "half").mkdir(parents=True)
    (cache / "half" / "part0.dat").write_bytes(b"abc")
    assert datasets_cli.main(["status"]) == 1
    out = capsys.readouterr().out
    assert dataset.pieces[1].url in out and dataset.pieces[0].url not in out
    assert "part1.dat has not been downloaded" in out
    (cache / "half" / "part0.dat").write_bytes(b"changed")
    assert datasets_cli.main(["status"]) == 1
    out = capsys.readouterr().out
    assert dataset.pieces[0].url in out and dataset.pieces[1].url in out


def test_status_of_a_present_dataset_prints_no_url(cache, monkeypatch, capsys):
    dataset = _dataset(b"abc")
    monkeypatch.setattr(catalogue, "ALL", (dataset,))
    fetch.fetch(dataset, opener=_serve(dataset, b"abc"))
    assert datasets_cli.main(["status"]) == 0
    assert "https://" not in capsys.readouterr().out


def test_status_of_a_stale_derived_dataset_names_the_archive_to_download(cache, monkeypatch, capsys):
    dataset, payload = _archive_dataset()
    monkeypatch.setattr(catalogue, "ALL", (dataset,))
    assert datasets_cli.main(["status"]) == 1
    assert dataset.pieces[0].url in capsys.readouterr().out
