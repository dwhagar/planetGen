#!/usr/bin/env python3
"""
Builds the API compatibility page (API.4) from the history of the API version.

    python scripts/build_api_compatibility.py            # write docs/reference/api-compatibility.md
    python scripts/build_api_compatibility.py --check    # build in memory and fail if the committed page differs

The history is `docs/reference/api-versions.json`: one row for every value `planetgen.api.version.API_VERSION` has had, with the release
that first carried it, the TODO items behind it, a one-line summary and the routes and payload formats it changed. The page is copied
from it and never typed by hand, so it cannot drift; `src/tests/test_api_compatibility.py` fails when it does, when the history does
not end at `API_VERSION`, or when a bump has no row. The stamp (`scripts/bump_version.py`) names the release of a row a PR added as
`pending` when it stamps the note that carries `api-version: changed`.
"""

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
HISTORY = os.path.join(HERE, "..", "docs", "reference", "api-versions.json")
OUT = os.path.join(HERE, "..", "docs", "reference", "api-compatibility.md")

CLIENTS = (
    ("The web pages", "`static/apiversion.js` (adds the version header to each call a page makes), `web/lib/apiclient.py`",
     "the API version of the release they ship in",
     "the server they were served by: pages and server are one release, so they always agree"),
    ("Remote generate (`generate.py --remote`, API.3)", "not built yet",
     "the API version of the release it is run from",
     "a server of the same API version (API.5 widens this to the current major and, for one release cycle after a major bump, the "
     "previous one, for upload payloads only)"),
    ("The command line (`planetgen ...`)", "`src/planetgen/cli/`", "none: it opens the database directly",
     "no server involved"),
)
"""tuple: `(client, where, speaks, talks to)` for every client of the API, as of this page."""


def load(path=HISTORY):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)["versions"]


def _release(row):
    return "next release" if row.get("pending") or not row.get("first_release") else row["first_release"]


def _releases_served(rows, index):
    """The releases whose pages speak `rows[index]`: from its first release to the one before the next version's."""
    row = rows[index]
    start = _release(row)
    if index + 1 == len(rows):
        return f"{start} and later"
    later = rows[index + 1]
    if later.get("pending"):
        return f"{start} up to the current release"
    if later["first_release"] == row["first_release"]:
        return f"none alone (replaced by version {later['version']} in {start})"
    return f"{start} up to the release before {later['first_release']}"


def render(rows):
    lines = [
        "# API compatibility",
        "",
        "<!-- Built from docs/reference/api-versions.json by scripts/build_api_compatibility.py. Do not edit by hand. -->",
        "",
        "The API has a version of its own, one integer that is not the release number (`src/planetgen/api/version.py`, API.22). Every "
        "`/api` answer carries it in the `X-PlanetGen-API-Version` header, `GET /api/health` reports it as `api_version`, and the "
        "`User-Agent` of everything the software sends names it (`PlanetGen-API/<n>`). It goes up by one, by hand, in the pull request "
        "that makes a breaking change to an endpoint: a removed or renamed route or field, a changed meaning or type, a new required "
        "parameter. Additive changes (a new route, a new optional field or parameter) leave it alone. That pull request also adds the "
        "row below to `docs/reference/api-versions.json`, a line to the change list in [api.md](api.md#api-version-and-object-ids), the "
        "route snapshot `src/tests/fixtures/api_routes.json` and `api-version: changed` to its release note; the stamp fills in the "
        "release. `src/tests/test_api_compatibility.py` fails when any of them is missing.",
        "",
        "## Versions",
        "",
        "| API version | First release | Items | Change | Routes and payload formats changed |",
        "|---|---|---|---|---|",
    ]
    for row in reversed(rows):
        changed = "<br>".join(row["changes"]) if row["changes"] else "none (the first numbered version)"
        lines.append(f"| **{row['version']}** | {_release(row)} | {', '.join(row['items'])} | {row['summary']} | {changed} |")
    lines += [
        "",
        "## Clients and the servers they talk to",
        "",
        "| Client | Where | Speaks | Talks to |",
        "|---|---|---|---|",
    ]
    for client, where, speaks, talks in CLIENTS:
        lines.append(f"| {client} | {where} | {speaks} | {talks} |")
    lines += [
        "",
        "## Which releases speak which version",
        "",
        "Every client that exists today ships with the server, so the release a page or a command comes from decides the version it speaks:",
        "",
        "| API version | Releases that speak it |",
        "|---|---|",
    ]
    for index in range(len(rows) - 1, -1, -1):
        lines.append(f"| **{rows[index]['version']}** | {_releases_served(rows, index)} |")
    lines += [
        "",
        "A client of one API version used against a server of another is not supported until API.5 adds the check and the converters "
        "for older upload payloads (`GET /api/version`, the `PlanetGen-API-Version` header both ways, a refusal `400 api_version_unsupported`).",
        "",
    ]
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("--check", action="store_true", help="fail if the committed page differs from the history")
    args = parser.parse_args(argv)
    text = render(load())
    if args.check:
        current = open(OUT, encoding="utf-8").read() if os.path.exists(OUT) else ""
        if current != text:
            print("docs/reference/api-compatibility.md is out of date; run python scripts/build_api_compatibility.py", file=sys.stderr)
            return 1
        return 0
    with open(OUT, "w", encoding="utf-8") as handle:
        handle.write(text)
    print(f"wrote {os.path.relpath(OUT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
