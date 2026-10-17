# tests/test_served.py

"""`tests.served.served` (TEST.193): the test web server waits for the requests still running when it stops."""

import threading
import time
import urllib.request
import warnings

from tests.served import InFlight, served


def _slow_app(started, release):
    def app(environ, start_response):
        started.set()
        release.wait(timeout=20)
        start_response("200 OK", [("Content-Type", "text/plain")])
        return [b"done"]
    return app


def _get(url, results):
    try:
        with urllib.request.urlopen(url, timeout=20) as response:
            results.append(response.read())
    except Exception as exc:  # the server stopped under it
        results.append(exc)


def test_it_serves_the_app_and_counts_nothing_afterwards():
    counted = []

    def app(environ, start_response):
        start_response("200 OK", [("Content-Type", "text/plain")])
        return [b"hello"]

    with served(app) as url:
        with urllib.request.urlopen(url, timeout=10) as response:
            counted.append(response.read())
    assert counted == [b"hello"]


def test_stopping_waits_for_a_request_that_is_still_running():
    started, release = threading.Event(), threading.Event()
    results = []
    began = time.monotonic()
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # a request that ends inside the wait raises no warning
        with served(_slow_app(started, release), drain_seconds=10) as url:
            client = threading.Thread(target=_get, args=(url, results), daemon=True)
            client.start()
            assert started.wait(timeout=10)
            threading.Timer(1.0, release.set).start()
    client.join(timeout=10)
    assert results == [b"done"]
    assert time.monotonic() - began >= 1.0, "the stop did not wait for the running request"


def test_stopping_gives_up_after_the_wait_and_names_the_request():
    started, release = threading.Event(), threading.Event()
    results = []
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            with served(_slow_app(started, release), drain_seconds=0.5) as url:
                client = threading.Thread(target=_get, args=(url, results), daemon=True)
                client.start()
                assert started.wait(timeout=10)
        assert any("1 request(s) were still running" in str(w.message) for w in caught), [str(w.message) for w in caught]
    finally:
        release.set()


def test_the_count_includes_a_request_that_raises():
    def app(environ, start_response):
        raise RuntimeError("boom")

    counted = InFlight(app)
    try:
        counted({}, lambda *a: None)
    except RuntimeError:
        pass
    assert counted.drain(0.1) == 0
