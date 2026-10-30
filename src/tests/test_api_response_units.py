# tests/test_api_response_units.py

"""API.33: every number an API answer returns states what it measures, over the real answers of the routes that serve a
galaxy (the request side is `test_api_field_units.py`). A number's key either ends in a unit suffix (`radius_km`, `age_gy`,
`luminosity_sol`) or is one of the plain kinds below, each with the noun it counts or the thing it names. A key that is
neither fails here, so a new field with a bare number is named or listed with a reason before it ships."""

import pytest

from planetgen.db import column_units

from tests.publicids import pid
from tests.test_api import _save_wide_binary_with_moons, client, seeded_sector  # noqa: F401 -- fixtures

NAME_UNITS = tuple(suffix for suffix, _unit in column_units.SUFFIXES) + ("_days", "_sol", "_earth", "_bytes", "_pc", "_pct", "_solar")

PLAIN = {
    # counts of the things the answer holds or paged past
    "count": "how many of the thing named beside it", "total": "how many rows match", "limit": "rows asked for in a page",
    "offset": "rows skipped before a page", "star_count": "stars", "system_count": "star systems",
    "interstellar_debris_count": "interstellar debris objects",
    # ids, grid addresses and versions: names, not quantities
    "id": "a row or object id", "shell_index": "grid shell", "band_index": "grid band", "wedge_index": "grid wedge",
    "api_version": "the API version", "schema_version": "the database schema version",
    "generator_version": "the generator epoch", "is_binary": "a 0 or 1 flag", "is_active": "a 0 or 1 flag",
    # shapes with no physical unit: a ratio with its range, a direction as the three parts of a unit vector
    "eccentricity": "orbital eccentricity: 0 is a circle, between 0 and 1 an ellipse",
    "binary_eccentricity": "orbital eccentricity of the pair: 0 is a circle, between 0 and 1 an ellipse",
    "spin_breakup_ratio": "spin rate as a share of the rate at which the body breaks up (1 is break-up)",
    "aerosol_draw": "the stored random number in [0, 1) that picked the aerosol optics within their brackets",
    "share": "a part's share of the whole it belongs to, 0 to 1",
    "spin_axis_x": "x part of the unit vector of the spin axis", "spin_axis_y": "y part of the unit vector of the spin axis",
    "spin_axis_z": "z part of the unit vector of the spin axis",
    # a unit the name states in the astronomer's own notation
    "metallicity_fe_h": "[Fe/H] in dex (the log of the iron to hydrogen ratio against the Sun)",
    "child_m": "how many sectors wide a cell of the galaxy map's drill is (1, 3, 9, 27, 81, 243)",
}
"""dict: The numbers an answer may carry with no unit suffix, each with what it is."""

PLAIN_ARRAYS = {"cells": "grid addresses [shell, band, wedge]"}
"""dict: Keys whose list holds lists of numbers that are addresses."""


REGISTERED = {}
for (_table, _column), _unit in column_units.COLUMNS.items():
    REGISTERED.setdefault(_column, set()).add(_unit)
UNSAID = ("ratio",)
"""tuple: Registry units that say nothing on their own: a column in one of them must be named in `PLAIN` with its reason."""


def _numbers(node, path, found, unit_missing):
    """Collects `(path of keys)` of every number in `node`, and the described quantities with a value and no unit."""
    if isinstance(node, dict):
        if isinstance(node.get("value"), (int, float)) and not isinstance(node["value"], bool) and "unit" not in node and "text" not in node:
            unit_missing.append("/".join(path))
        for key, value in node.items():
            _numbers(value, path + [key], found, unit_missing)
    elif isinstance(node, list):
        for value in node:
            _numbers(value, path, found, unit_missing)
    elif isinstance(node, (int, float)) and not isinstance(node, bool):
        found.add(tuple(path))


def _says_what(key):
    """Whether a key says what its number measures: a unit suffix, a column the registry gives a unit that says something
    (a fraction, a score of 9, an optical depth, a pH ...), or a listed kind."""
    if key.endswith(NAME_UNITS) or key in PLAIN:
        return True
    return bool(REGISTERED.get(key)) and not REGISTERED[key] <= set(UNSAID)


def _unstated(path):
    """A number is stated by its own key or by the key of a map above it (`gas_partials_pa` over its gases), and the two
    halves of a described quantity are checked beside their unit."""
    return not any(_says_what(key) for key in path) and path[-1] not in ("sol", "value") and not (
        len(path) > 1 and path[-2] in PLAIN_ARRAYS)


def test_every_number_an_answer_returns_names_its_unit_or_is_a_listed_kind(client, seeded_sector):
    config, sector_id, system_ids = seeded_sector
    sector = pid("sector", sector_id)
    paths = ["/api/sectors", f"/api/sectors/{sector}", f"/api/systems?sector_id={sector}", "/api/search?spectral=G",
             "/api/galaxy/shape", "/api/galaxy/stage", "/api/galaxy/sectors", "/api/health", "/api/galaxy/superbubbles",
             "/api/galaxy/cell?x=9000&y=-120&z=4000", "/api/galaxy/bright-stars?shell=3&band=0&wedge=1&all=1",
             "/api/galaxy/uncharted?shell=3&band=0&wedge=1"]
    for system_id in system_ids:
        paths += [f"/api/systems/{pid('system', system_id)}", f"/api/systems/{pid('system', system_id)}/description"]
    with_bodies = pid("system", _save_wide_binary_with_moons(config))
    paths += [f"/api/systems/{with_bodies}", f"/api/systems/{with_bodies}/description"]
    system = client.get(f"/api/systems/{with_bodies}").get_json()
    assert system["planets"] and system["stars"], "the system read has planets and stars, whose columns are the bulk of the numbers"
    unstated, valueless = [], []
    for path in paths:
        response = client.get(path)
        assert response.status_code == 200, path
        found = set()
        _numbers(response.get_json(), [], found, valueless)
        unstated += [f"{path}: {'.'.join(keys)}" for keys in sorted(found) if _unstated(keys)]
    assert unstated == [], "a number in an answer that names no unit and is not a listed kind: " + "; ".join(unstated)
    assert valueless == [], "a described quantity with no unit beside its value: " + "; ".join(valueless)


def test_every_listed_kind_is_a_noun_or_a_notation_not_an_empty_reason():
    assert all(len(reason) > 3 for reason in {**PLAIN, **PLAIN_ARRAYS}.values())
    assert not [key for key in PLAIN if key.endswith(NAME_UNITS)], "a key that already names a unit does not need listing"


def test_a_backup_archive_answers_its_sizes_in_bytes():
    from types import SimpleNamespace

    from planetgen.api import backups

    part = SimpleNamespace(name="settings", size=1234, files=[1, 2], known=True)
    archive = SimpleNamespace(name="a", taken_at="t", taken_by="u", reason="r", version="8.0", size=5678, complete=True,
                              readable=True, parts=[part])
    view = backups._archive_view(archive)
    assert view["size_bytes"] == 5678 and "size" not in view
    assert view["parts"][0]["size_bytes"] == 1234 and "size" not in view["parts"][0]
