# tests/test_api_compatibility.py

"""
The API compatibility data (API.4): `docs/reference/api-versions.json` holds one row for every value of `planetgen.api.version.API_VERSION`,
`docs/reference/api-compatibility.md` is built from it by `scripts/build_api_compatibility.py`, and the change list of
`docs/reference/api.md` has a line for every version. A PR that raises `API_VERSION` fails here until all of them move with it.
"""

import importlib.util
import json
import os
import re

from planetgen.api.version import API_VERSION

ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
_spec = importlib.util.spec_from_file_location("build_api_compatibility", os.path.join(ROOT, "scripts", "build_api_compatibility.py"))
build = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(build)

RELEASE_RE = re.compile(r"^\d+\.\d+\.\d+$")
ROUTE_CHANGE_RE = re.compile(r"^(?P<method>GET|POST|PUT|PATCH|DELETE) (?P<path>/\S+) removed$")


def _release_key(release):
    return tuple(int(part) for part in release.split("."))


def test_the_history_runs_from_one_to_the_api_version_without_a_gap():
    rows = build.load()
    assert [row["version"] for row in rows] == list(range(1, API_VERSION + 1)), (
        "docs/reference/api-versions.json must have one row for each API version, 1 to API_VERSION "
        "(a PR that raises the version adds its row, with first_release null and pending true, and 'api-version: changed' in its note)")


def test_every_row_says_what_it_is_and_only_the_newest_may_wait_for_its_release():
    rows = build.load()
    for index, row in enumerate(rows):
        assert row["items"] and row["summary"] and isinstance(row["changes"], list), row
        if row.get("pending"):
            assert index == len(rows) - 1 and row["first_release"] is None, f"only the newest row may be pending: {row}"
        else:
            assert RELEASE_RE.match(row["first_release"]), row
    released = [_release_key(row["first_release"]) for row in rows if not row.get("pending")]
    assert released == sorted(released), "the first releases must not go backwards"


def test_the_compatibility_page_is_built_from_the_history():
    with open(build.OUT, encoding="utf-8") as handle:
        assert handle.read() == build.render(build.load()), "run python scripts/build_api_compatibility.py"


def test_the_change_list_of_the_api_reference_has_a_line_for_every_version():
    with open(os.path.join(ROOT, "docs", "reference", "api.md"), encoding="utf-8") as handle:
        text = handle.read()
    listed = {int(number) for number in re.findall(r"^- \*\*(\d+)\*\*", text, re.M)}
    assert listed == set(range(1, API_VERSION + 1)), "docs/reference/api.md needs a '- **N**' line for each API version"


def test_a_route_the_history_says_was_removed_is_not_in_the_route_snapshot():
    with open(os.path.join(os.path.dirname(__file__), "fixtures", "api_routes.json"), encoding="utf-8") as handle:
        routes = set(json.load(handle)["routes"])
    for row in build.load():
        for change in row["changes"]:
            match = ROUTE_CHANGE_RE.match(change)
            if match:
                assert f"{match['path']} {match['method']}" not in routes, f"version {row['version']} removed {change}, but the route is served"
