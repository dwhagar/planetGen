# planetgen/api/version.py

"""
The API version (API.22): one sequential integer, like the database schema number.

A pull request that makes a breaking change to an endpoint (a removed or renamed
route or field, a changed meaning or type, a new required parameter) adds one
here, adds a line to the change list in `docs/api.md` and updates the route
snapshot `tests/fixtures/api_routes.json` (`tests/test_api_version.py` fails
until the snapshot and this number agree). Additive changes (a new route, a new
optional field or parameter) leave it alone. The integer is returned by
`GET /api/health`, sent in the `X-PlanetGen-API-Version` header of every `/api`
response, and shown on the admin Stats page. The remote-run handshake (API.17)
and the compatibility data (API.4) compare it.
"""

API_VERSION = 3
"""int: The version of the API this code serves."""

API_VERSION_HEADER = "X-PlanetGen-API-Version"
"""str: The response header carrying `API_VERSION`."""
