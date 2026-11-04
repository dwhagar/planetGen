# tests/test_served.py

"""`tests.served.served` (TEST.193): the test web server waits for the requests still running when it stops."""

import logging
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


def test_the_warning_names_the_route_still_running():
    started, release = threading.Event(), threading.Event()
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            with served(_slow_app(started, release), drain_seconds=0.5) as url:
                client = threading.Thread(target=_get, args=(url + "/admin/generate/shape-report?mass=1", []), daemon=True)
                client.start()
                assert started.wait(timeout=10)
        assert any("GET /admin/generate/shape-report x1" in str(w.message) for w in caught), [str(w.message) for w in caught]
    finally:
        release.set()


def test_a_static_file_still_being_read_is_not_waited_for_or_named():
    started, release = threading.Event(), threading.Event()
    began = time.monotonic()
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            with served(_slow_app(started, release), drain_seconds=10) as url:
                client = threading.Thread(target=_get, args=(url + "/static/vendor/shoelace/chunks/chunk.js", []), daemon=True)
                client.start()
                assert started.wait(timeout=10)
        assert time.monotonic() - began < 5
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


def test_the_request_log_stays_out_of_a_failing_tests_output(capfd, caplog):
    """TEST.187: werkzeug logs one line per request, which buried a failing browser test's assertion."""
    logger = logging.getLogger("werkzeug")
    saved = (logger.level, list(logger.handlers))
    logger.setLevel(logging.NOTSET)  # as werkzeug finds it before its first request
    logger.handlers[:] = []
    caplog.set_level(logging.DEBUG)
    try:
        def app(environ, start_response):
            start_response("200 OK", [("Content-Type", "text/plain")])
            return [b"hello"]

        with served(app) as url:
            for _ in range(3):
                urllib.request.urlopen(f"{url}/static/style.css", timeout=10).read()
        assert logger.level == logging.WARNING
        assert "GET /static/style.css" not in capfd.readouterr().err
        assert not [record for record in caplog.records if record.name == "werkzeug"]
    finally:
        logger.setLevel(saved[0])
        logger.handlers[:] = saved[1]
