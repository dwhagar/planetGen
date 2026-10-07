# tests/webpage_support.py

"""
`live_api(mysql_config)`: a pytest fixture that starts `planetgen.web.app.create_app()`
on a real background thread (`werkzeug.serving.make_server`, not Flask's
own `app.run()`, which blocks) bound to `127.0.0.1:0` (OS-assigned free
port), pointed at the fixture's throwaway database. Yields the API's base
URL; shuts the server down on teardown. For tests that need the API over
real HTTP; the pages themselves are tested with Flask's test client.
"""

import threading

import pytest
from werkzeug.serving import make_server


@pytest.fixture
def live_api(mysql_config):
    """Starts the real Flask API on a background thread against
    mysql_config's throwaway database; yields its base URL."""
    from planetgen.web.app import create_app
    from planetgen.api.config import Config

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
