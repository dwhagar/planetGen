# tests/served.py

"""
A test web server that is finished when its fixture ends (TEST.193).

`server.shutdown()` stops the accept loop only. A request a browser had already sent keeps running in its own
thread, holding a pooled database connection (and an open transaction), so the module fixture that dropped the
test database right after it waited on that connection's lock: 853 s in gate runs. `served` counts the requests
in flight and, when the server stops, waits a few seconds for them to end before the fixture goes on to drop its
database. A request still running after that is named in a warning, and `dbpool.drop_database` ends its
connection.
"""

import contextlib
import logging
import threading
import time
import warnings

from werkzeug.serving import make_server

def make_quiet_server(host, port, app, **kwargs):
    """`werkzeug.serving.make_server` with the `werkzeug` logger at WARNING (TEST.187): its one INFO line per
    request ("GET /static/... 200") otherwise fills the captured log, and a failing browser test's report ends
    in hundreds of them, below the assertion. Every test server is made here."""
    logging.getLogger("werkzeug").setLevel(logging.WARNING)
    return make_server(host, port, app, **kwargs)


DRAIN_SECONDS = 8
"""How long the stop waits for requests still running when the server is shut down."""


class _Counted:
    """The response iterable, counted as finished when the server closes it."""

    def __init__(self, result, done):
        self._result = result
        self._done = done

    def __iter__(self):
        return iter(self._result)

    def close(self):
        try:
            close = getattr(self._result, "close", None)
            if close:
                close()
        finally:
            self._done()


STATIC_PREFIX = "/static/"
"""A request under this path is a file read: it holds no database connection, so it is not counted (TEST.219: the leftovers
of a browser test page were a few of the Shoelace files the closing page had stopped reading, which no wait ends)."""


class InFlight:
    """WSGI middleware counting the requests that have started and not yet finished, those that can hold a database
    connection (not `STATIC_PREFIX`)."""

    def __init__(self, app):
        self.app = app
        self.active = 0
        self.routes = {}
        self._cond = threading.Condition()

    def _finished(self, key):
        with self._cond:
            self.active -= 1
            self.routes[key] -= 1
            self._cond.notify_all()

    def __call__(self, environ, start_response):
        if environ.get("PATH_INFO", "").startswith(STATIC_PREFIX):
            return self.app(environ, start_response)
        key = f"{environ.get('REQUEST_METHOD', '?')} {environ.get('PATH_INFO', '?')}"
        with self._cond:
            self.active += 1
            self.routes[key] = self.routes.get(key, 0) + 1
        try:
            result = self.app(environ, start_response)
        except BaseException:
            self._finished(key)
            raise
        return _Counted(result, lambda: self._finished(key))

    def running_routes(self):
        """`{"METHOD /path": how many}` of the requests started and not yet finished."""
        with self._cond:
            return {key: count for key, count in self.routes.items() if count > 0}

    def drain(self, seconds):
        """Waits up to `seconds` for every started request to finish; returns how many are still running."""
        deadline = time.monotonic() + seconds
        with self._cond:
            while self.active > 0 and time.monotonic() < deadline:
                self._cond.wait(timeout=max(deadline - time.monotonic(), 0.01))
            return self.active


@contextlib.contextmanager
def served(app, drain_seconds=DRAIN_SECONDS):
    """Runs `app` on a free local port in a thread and yields its base URL (`http://127.0.0.1:<port>`)."""
    counted = InFlight(app)
    server = make_quiet_server("127.0.0.1", 0, counted, threaded=True)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield f"http://127.0.0.1:{server.server_port}"
    finally:
        server.shutdown()
        thread.join(timeout=10)
        running = counted.drain(drain_seconds)
        if running:
            named = ", ".join(f"{key} x{count}" for key, count in sorted(counted.running_routes().items()))
            warnings.warn(f"TEST.193: {running} request(s) were still running {drain_seconds} s after the test "
                          f"server stopped ({named}); the connections they hold are ended when the database is dropped")
        server.server_close()
