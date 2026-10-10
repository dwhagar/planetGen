# tests/test_api_version.py

"""
The API version (API.22, `planetgen/api/version.py`): every `/api` answer carries it, `/api/health`
reports it, and it cannot sit still while the route table moves. `tests/fixtures/api_routes.json`
holds the version and the routes it was last written for; a change to either needs the other:
bump `API_VERSION` for a breaking change (and add its line to the change list in `docs/api.md`),
then rewrite the snapshot with `python tests/test_api_version.py --write`.
"""

import json
import pathlib
import sys

from planetgen.api.config import Config
from planetgen.api.version import API_VERSION, API_VERSION_HEADER
from planetgen.web.app import create_app

SNAPSHOT = pathlib.Path(__file__).parent / "fixtures" / "api_routes.json"


def _live_routes(app):
    return sorted(f"{rule.rule} {','.join(sorted(rule.methods - {'HEAD', 'OPTIONS'}))}"
                  for rule in app.url_map.iter_rules() if rule.rule.startswith("/api/"))


class _Config(Config):
    SECRET_KEY = "api-version-test"
    SESSION_COOKIE_SECURE = False


def _app():
    app = create_app(_Config)
    app.testing = True
    return app


def test_every_api_answer_carries_the_version():
    client = _app().test_client()
    for path in ("/api/health", "/api/nope"):
        assert client.get(path).headers[API_VERSION_HEADER] == str(API_VERSION)


def test_the_route_table_and_the_version_move_together():
    snapshot = json.loads(SNAPSHOT.read_text(encoding="utf-8"))
    live = _live_routes(_app())
    assert live == snapshot["routes"] and API_VERSION == snapshot["api_version"], (
        "The API's routes changed. If that is a breaking change, bump API_VERSION in planetgen/api/version.py and "
        "add a line to docs/api.md; either way rewrite the snapshot: python tests/test_api_version.py --write\n"
        f"added: {sorted(set(live) - set(snapshot['routes']))}\nremoved: {sorted(set(snapshot['routes']) - set(live))}"
        f"\nsnapshot version {snapshot['api_version']}, code version {API_VERSION}")


if __name__ == "__main__":
    if "--write" in sys.argv:
        SNAPSHOT.write_text(json.dumps({"api_version": API_VERSION, "routes": _live_routes(_app())}, indent=1) + "\n",
                            encoding="utf-8")
        print(f"wrote {SNAPSHOT}")
