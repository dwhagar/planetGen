# tests/webpage_support.py

"""
Subprocess-based test harness for the `src/html/*.py` CGI browser pages.

These pages have **no other test seam**: every one of them is plain
stdlib CGI (reads `os.environ`, writes headers+body to stdout) with
module-level code that runs a real request at import time (auth checks,
POST-body handling, then `run(handler)`) -- there's no factory function,
no `if __name__ == "__main__":` guard on several of them (e.g.
`system.py`), and they aren't WSGI apps, so no Flask/Werkzeug test client
applies. `docs/html-interface.md`'s own "Local testing without Apache"
section documents the only way that actually works: start the JSON API
separately, then invoke a page script directly with `QUERY_STRING`/
`CONTENT_LENGTH`/env set, exactly the way a real CGI server would. This
module automates that pattern for pytest:

- `live_api(mysql_config)`: starts `api.app.create_app()` on a real
  background thread (`werkzeug.serving.make_server`, not Flask's own
  `app.run()`, which blocks) bound to `127.0.0.1:0` (OS-assigned free
  port), pointed at the fixture's throwaway database, admin bootstrapped.
  Yields the base URL; shuts the server down on teardown.
- `run_page(base_url, page, query="", method="GET", body="", cookie=None,
  db=None)`: subprocess-invokes `python3 src/html/<page>.py` with the
  right CGI environment, returns a `PageResult(status_code, status_text,
  headers, body)` parsed from its stdout.
"""

from __future__ import annotations

import os
import subprocess
import sys
import threading
from dataclasses import dataclass, field
from urllib.parse import urlencode

import pytest
from werkzeug.serving import make_server

_HTML_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "html")
_PYTHON = sys.executable


@dataclass
class PageResult:
    status_code: int
    status_text: str
    headers: dict = field(default_factory=dict)
    body: str = ""
    stderr: str = ""

    @property
    def ok(self) -> bool:
        return 200 <= self.status_code < 300


def _parse_cgi_output(raw: bytes) -> PageResult:
    """Splits a CGI script's stdout into (status, headers, body) -- the
    header block is CRLF-terminated key: value lines up to the first
    blank line (page.py's own send_headers format)."""
    text = raw.decode("utf-8", errors="replace")
    if "\r\n\r\n" in text:
        header_block, body = text.split("\r\n\r\n", 1)
    elif "\n\n" in text:
        header_block, body = text.split("\n\n", 1)
    else:
        header_block, body = text, ""

    headers = {}
    status_code, status_text = 0, ""
    for line in header_block.splitlines():
        if ":" not in line:
            continue
        key, _, value = line.partition(":")
        key, value = key.strip(), value.strip()
        headers[key] = value
        if key.lower() == "status":
            parts = value.split(" ", 1)
            status_code = int(parts[0])
            status_text = parts[1] if len(parts) > 1 else ""
    return PageResult(status_code=status_code, status_text=status_text, headers=headers, body=body)


def run_page(base_url, page, query="", method="GET", body="", cookie=None, timeout=30, extra_env=None) -> PageResult:
    """
    Runs `src/html/<page>` as a real CGI subprocess against a live API.

    Args:
        base_url (str): The live_api fixture's base URL.
        page (str): Script filename, e.g. "system.py".
        query (str or dict): QUERY_STRING, or a dict to urlencode.
        method (str): "GET" or "POST".
        body (str or dict): POST body (application/x-www-form-urlencoded),
            or a dict to urlencode.
        cookie (str, optional): Raw Cookie header value (session cookie).
        extra_env (dict, optional): Additional/overriding env vars.

    Returns:
        PageResult
    """
    if isinstance(query, dict):
        query = urlencode(query)
    if isinstance(body, dict):
        body = urlencode(body)
    body_bytes = body.encode("utf-8")

    env = dict(os.environ)
    env["PLANETGEN_API_BASE_URL"] = base_url
    env["QUERY_STRING"] = query
    env["REQUEST_METHOD"] = method
    env["CONTENT_LENGTH"] = str(len(body_bytes)) if method == "POST" else ""
    env["CONTENT_TYPE"] = "application/x-www-form-urlencoded"
    if cookie:
        env["HTTP_COOKIE"] = cookie
    else:
        env.pop("HTTP_COOKIE", None)
    if extra_env:
        env.update(extra_env)

    script_path = os.path.join(_HTML_DIR, page)
    proc = subprocess.run(
        [_PYTHON, script_path],
        input=body_bytes,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=env,
        timeout=timeout,
    )
    result = _parse_cgi_output(proc.stdout)
    result.stderr = proc.stderr.decode("utf-8", errors="replace")
    return result


@pytest.fixture
def live_api(mysql_config):
    """Starts the real Flask API on a background thread against
    mysql_config's throwaway database; yields its base URL."""
    from api.app import create_app
    from api.config import Config

    class TestConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config

    app = create_app(TestConfig)
    server = make_server("127.0.0.1", 0, app)
    port = server.server_port
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield f"http://127.0.0.1:{port}/api"
    finally:
        server.shutdown()
        thread.join(timeout=5)
