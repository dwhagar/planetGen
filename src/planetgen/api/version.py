# planetgen/api/version.py

"""
The API version (API.22): one sequential integer, like the database schema number.

A pull request that makes a breaking change to an endpoint (a removed or renamed
route or field, a changed meaning or type, a new required parameter) adds one
here, adds a line to the change list in `docs/reference/api.md`, adds its row to
`docs/reference/api-versions.json` (`"pending": true`, `"first_release": null`; the compatibility page
`docs/reference/api-compatibility.md` is built from it, API.4), puts `api-version: changed` in its release note (the stamp names
the release) and updates the route snapshot `tests/fixtures/api_routes.json` (`tests/test_api_version.py` fails
until the snapshot and this number agree, `tests/test_api_compatibility.py` until the history and the page do). Additive changes (a new route, a new
optional field or parameter) leave it alone. The integer is returned by
`GET /api/health`, sent in the `X-PlanetGen-API-Version` header of every `/api`
response, and shown on the admin Stats page. The remote-run handshake (API.17)
and the compatibility data (API.4, `docs/reference/api-compatibility.md`) compare it.
"""

API_VERSION = 12
"""int: The version of the API this code serves."""

API_VERSION_HEADER = "X-PlanetGen-API-Version"
"""str: The response header carrying `API_VERSION`."""
