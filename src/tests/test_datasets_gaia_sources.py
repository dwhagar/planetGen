# tests/test_datasets_gaia_sources.py

"""OPS.183: the Gaia sample (`datasets/gaia_nearby.py`) is fetched from the archive, then its asynchronous query, the Heidelberg mirror and two VizieR queries in turn; setup fails only when all of them do."""

import dataclasses
import hashlib
import io
import urllib.error
import urllib.parse

import pytest

from planetgen.datasets import fetch, gaia_nearby

ROWS = 5323


def _file(rows=ROWS):
    lines = [",".join(gaia_nearby.COLUMNS)]
    for i in range(rows):
        lines.append(f"{1000 + i},{(i * 7) % 360}.25,{(i % 90)}.5,{50 + (i % 40) + 0.125},0.0{(i % 9) + 1},1.5,-2.5,{9 + (i % 5)}.5,{1 + (i % 3)}.25")
    return ("\n".join(lines) + "\n").encode("utf-8")


def _vizier_file(rows=ROWS):
    lines = ["Source,RA_ICRS,DE_ICRS,Plx,e_Plx,pmRA,pmDE,Gmag,BP-RP"]
    for i in range(rows):
        lines.append(f"{1000 + i},{(i * 7) % 360}.25,{(i % 90)}.5,{50 + (i % 40) + 0.125},0.0{(i % 9) + 1},1.5,-2.5,{9 + (i % 5)}.5,{1 + (i % 3)}.25")
    return ("\n".join(lines) + "\n").encode("utf-8")


class _Response(io.BytesIO):
    def __init__(self, data, url):
        super().__init__(data)
        self.url = url

    def geturl(self):
        return self.url

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


class _Net:
    """An opener serving what a test says each host answers: a bytes (served), an exception (raised) or a function of the request."""

    def __init__(self, answers):
        self.answers, self.requests = answers, []

    def __call__(self, request, timeout=None):
        url = request.full_url
        self.requests.append(url)
        for prefix, answer in self.answers.items():
            if url.startswith(prefix):
                if callable(answer):
                    answer = answer(request)
                if isinstance(answer, Exception):
                    raise answer
                return _Response(answer, url)
        raise urllib.error.URLError(f"nothing serves {url}")


ESA_SYNC = gaia_nearby.URL.split("?")[0]
ESA_ASYNC = gaia_nearby.ESA_TAP + "/async"
JOB = ESA_ASYNC + "/abc123"
ARI = gaia_nearby.ARI_TAP
VIZIER = gaia_nearby.VIZIER_TAP


def _http(code):
    return urllib.error.HTTPError("https://x", code, "no", {}, None)


def _dataset(payload):
    piece = fetch.Piece(gaia_nearby.FILE, gaia_nearby.URL, hashlib.sha256(payload).hexdigest(), len(payload), mirrors=gaia_nearby.MIRRORS, describe=gaia_nearby.describe_rows)
    return dataclasses.replace(gaia_nearby.DATASET, pieces=(piece,))


def _async_ok(payload):
    def post(request):
        return b"" if request.data is not None and b"PHASE=RUN" not in request.data else b"<job/>"
    return {ESA_ASYNC + "/abc123/phase": b"COMPLETED", ESA_ASYNC + "/abc123/results/result": payload, JOB: post, ESA_ASYNC: lambda r: b"<job/>"}


@pytest.fixture(autouse=True)
def cache(tmp_path, monkeypatch):
    monkeypatch.setenv(fetch.CACHE_DIR_ENV, str(tmp_path / "cache"))
    fetch._VERIFIED.clear()
    monkeypatch.setattr(fetch.time, "sleep", lambda seconds: None)
    monkeypatch.setattr(gaia_nearby.time, "sleep", lambda seconds: None)
    monkeypatch.delenv(gaia_nearby.ESA_USER_ENV, raising=False)
    monkeypatch.delenv(gaia_nearby.ESA_PASSWORD_ENV, raising=False)


class _AsyncNet(_Net):
    """The ESA archive's asynchronous query: the job address comes back from the POST, the phase and the result from the GETs."""

    def __init__(self, answers, payload, phase="COMPLETED"):
        super().__init__(answers)
        self.payload, self.phase = payload, phase

    def __call__(self, request, timeout=None):
        url = request.full_url
        if url == ESA_ASYNC:
            self.requests.append(url)
            return _Response(b"<job/>", JOB)
        if url == JOB + "/phase":
            self.requests.append(url)
            return _Response(self.phase.encode(), url)
        if url == JOB + "/results/result":
            self.requests.append(url)
            return _Response(self.payload, url)
        if url == JOB:
            self.requests.append(url + "#delete")
            return _Response(b"", url)
        return super().__call__(request, timeout)


def _fetch(payload, net):
    dataset = _dataset(payload)
    return dataset, fetch.fetch(dataset, opener=net)


def test_the_archives_own_answer_is_used_when_it_is_whole():
    payload = _file()
    net = _Net({ESA_SYNC: payload})
    _dataset_made, result = _fetch(payload, net)
    assert result.status == "fetched" and len(net.requests) == 1


def test_a_capped_answer_goes_to_the_asynchronous_query_which_serves_the_same_bytes():
    payload = _file()
    net = _AsyncNet({ESA_SYNC: _file(2000)}, payload)
    dataset, result = _fetch(payload, net)
    assert result.status == "fetched"
    assert fetch.locate(dataset).path(gaia_nearby.FILE)
    assert JOB + "#delete" in net.requests, "the finished job is deleted (the anonymous queue holds 50)"


def test_a_full_queue_and_a_failed_job_go_to_the_heidelberg_mirror():
    payload = _file()
    net = _AsyncNet({ESA_SYNC: _http(429), ARI: payload}, b"", phase="ERROR")
    _dataset_made, result = _fetch(payload, net)
    assert result.status == "fetched"
    assert any(url.startswith(ARI) for url in net.requests)


def test_the_vizier_file_is_taken_when_the_archive_and_the_mirror_fail_and_is_kept_by_the_next_setup(cache):
    payload = _file()
    net = _AsyncNet({ESA_SYNC: _http(429), ARI: _http(503), VIZIER: _vizier_file()}, b"", phase="ABORTED")
    dataset = _dataset(payload)
    assert fetch.fetch(dataset, opener=net).status == "fetched"
    path = fetch.locate(dataset).path(gaia_nearby.FILE)
    assert open(path, "rb").read().splitlines()[0] == b",".join(c.encode() for c in gaia_nearby.COLUMNS)
    asked = len(net.requests)
    assert fetch.fetch(dataset, opener=net).status == "present" and len(net.requests) == asked, "an accepted alternative is kept, not fetched again"
    assert not fetch.locate(dataset).missing


def test_the_nearby_stars_catalogue_is_the_last_resort():
    payload = _file()
    gcns = _vizier_file().replace(b"Source,", b"GaiaEDR3,", 1)
    answers = {ESA_SYNC: _http(429), ARI: _http(503), VIZIER: lambda request: gcns if b"gcnsmain" in request.data else _http(500)}
    net = _AsyncNet(answers, b"", phase="ERROR")
    dataset, result = _fetch(payload, net)
    assert result.status == "fetched"


def test_setup_fails_only_when_all_four_sources_fail_and_says_why_for_each():
    payload = _file()
    net = _AsyncNet({ESA_SYNC: _http(429), ARI: _http(503), VIZIER: _http(500)}, b"", phase="ERROR")
    dataset, result = _fetch(payload, net)
    assert result.status == "missing"
    reason = "; ".join(result.reasons)
    for label in ("its own address", "asynchronous query", "Heidelberg", "VizieR (Gaia DR3", "nearby-star catalog"):
        assert label in reason


def test_a_vizier_answer_with_too_few_rows_or_the_wrong_columns_is_refused():
    assert gaia_nearby.accept_rows(_file(2000)) and "capped" in gaia_nearby.accept_rows(_file(2000))
    assert gaia_nearby.accept_rows(b"a,b\n1,2\n")
    assert gaia_nearby.accept_rows(_file()) is None
    with pytest.raises(ValueError):
        gaia_nearby.vizier_to_file(b"Source,Plx\n1,50\n")


def test_the_vizier_file_keeps_only_well_measured_parallaxes_sorted_by_source():
    raw = (b"Source,RA_ICRS,DE_ICRS,Plx,e_Plx,pmRA,pmDE,Gmag,BP-RP\n"
           b"30,1,2,60,1,0,0,10,1\n20,1,2,60,10,0,0,10,1\n10,1,2,45,2,0,0,10,1\n")
    made = gaia_nearby.vizier_to_file(raw).decode().splitlines()
    assert made[0].split(",") == list(gaia_nearby.COLUMNS)
    assert [line.split(",")[0] for line in made[1:]] == ["10", "30"]


def test_the_account_password_is_registered_so_no_log_shows_it(monkeypatch):
    from planetgen.util import secrets_mask
    secrets_mask.forget_registered()
    payload = _file()
    password = "Esa-Account-Pw-9182736"
    net = _AsyncNet({"https://gea.esac.esa.int/tap-server/login": b""}, payload)
    monkeypatch.setattr(gaia_nearby.urllib.request, "build_opener", lambda *a, **k: types_opener(net))
    gaia_nearby.esa_async(net, 5, environ={gaia_nearby.ESA_USER_ENV: "someone", gaia_nearby.ESA_PASSWORD_ENV: password})
    assert password not in secrets_mask.mask_text(f"the account said {password}")
    secrets_mask.forget_registered()


def types_opener(net):
    class _Built:
        open = staticmethod(net)
    return _Built
