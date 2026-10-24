# tests/test_property_sheet.py

"""
UX.168: the property sheet registry (`describe/sheet.py`) lists every stored column of a star, planet, moon, belt and
comet, so a column a later change adds cannot be forgotten, and `GET /api/objects/<ref>/sheet` serves it.
"""

import pytest

from planetgen.db import column_units, models
from planetgen.describe import sheet

NUMERIC = ("INTEGER", "FLOAT", "DOUBLE", "TINYINT", "SMALLINT", "BIGINT", "NUMERIC", "DECIMAL")


def _columns(kind):
    return [column.name for column in models.metadata.tables[sheet.TABLES[kind]].columns]


@pytest.mark.parametrize("kind", sorted(sheet.SPECS))
def test_every_column_is_on_the_sheet_or_says_why_it_is_not(kind):
    shown = {spec.stored_column() for spec in sheet.SPECS[kind]}
    hidden = sheet.NOT_SHOWN[kind]
    columns = set(_columns(kind))
    missing = columns - shown - set(hidden)
    assert not missing, (
        f"{sheet.TABLES[kind]} has columns the {kind} property sheet neither shows nor excuses: {sorted(missing)}. "
        f"Add a row to SPECS[{kind!r}] in planetgen/describe/sheet.py, or an entry with the reason to NOT_SHOWN[{kind!r}].")
    assert not shown & set(hidden), f"on the sheet and in NOT_SHOWN: {sorted(shown & set(hidden))}"
    assert not (shown | set(hidden)) - columns, f"not columns of {sheet.TABLES[kind]}: {sorted((shown | set(hidden)) - columns)}"
    assert all(reason.strip() for reason in hidden.values())


@pytest.mark.parametrize("kind", sorted(sheet.SPECS))
def test_a_row_has_a_label_a_known_group_and_one_place(kind):
    keys = [spec.key for spec in sheet.SPECS[kind]]
    assert len(keys) == len(set(keys))
    for spec in sheet.SPECS[kind]:
        assert spec.label and spec.group in sheet.GROUPS
        if spec.quantity == sheet.AUTO or spec.quantity not in (sheet.TEXT, sheet.FLAG, sheet.TIER):
            assert (sheet.TABLES[kind], spec.stored_column()) in column_units.COLUMNS, f"{kind}.{spec.key} has no unit"


def test_the_kinds_are_the_bodies_and_phenomena_a_page_exists_for():
    from planetgen.galaxy import objectref

    assert sorted(sheet.SPECS) == sorted(objectref.BODY_KINDS + objectref.PHENOMENON_KINDS) == sorted(sheet.TABLES)


def test_a_unit_decides_how_a_number_is_shown():
    spec = sheet.PropertySpec("mass_kg", "Mass", "physical")
    assert sheet._value(spec, "planets", 5.97e24)[0]["qty"]["kind"] == "mass"
    assert sheet._value(spec, "planets", None) == [sheet.NOT_STORED]
    dose = sheet.PropertySpec("surface_dose_msv_yr", "Dose", "radiation")
    assert sheet._value(dose, "planets", 2.4) == ["2.4 mSv/yr"]
    flag = sheet.PropertySpec("ozone_loss_flag", "Ozone", "atmosphere", sheet.FLAG)
    assert sheet._value(flag, "planets", 1) == ["Yes"] and sheet._value(flag, "planets", 0) == ["No"]
    tier = sheet.PropertySpec("tier_air", "Air", "habitability", sheet.TIER)
    assert sheet._value(tier, "planets", 1) == ["Green"]
    star = sheet.PropertySpec("mass_kg", "Mass", "physical", solar=True)
    qty = sheet._value(star, "stars", 1.98841e30)[0]["qty"]
    assert qty["sol"] == pytest.approx(1.0, rel=1e-4)


def test_an_unknown_kind_has_no_sheet():
    with pytest.raises(ValueError):
        sheet.specs_for("system")


# The API, over a real database ---------------------------------------------------------------------------------

from planetgen.db import store  # noqa: E402
from tests.publicids import pid  # noqa: E402
from tests.test_api import _save_wide_binary_with_moons, client  # noqa: E402,F401


def _ids(mysql_config, sql, params=()):
    conn = store.get_connection(mysql_config)
    try:
        return [row["id"] for row in conn.execute(sql, params).fetchall()]
    finally:
        conn.close()


def _rows(body):
    return {row["key"]: row for group in body["groups"] for row in group["rows"]}


def test_the_sheet_of_a_star_a_planet_and_a_moon_returns_every_property_with_units(client, mysql_config):
    system_id = _save_wide_binary_with_moons(mysql_config)
    for kind, table in (("star", "stars"), ("planet", "planets"), ("moon", "moons")):
        row_id = _ids(mysql_config, f"SELECT id FROM {table} WHERE star_system_id = ? ORDER BY id", (system_id,))[0]
        response = client.get(f"/api/objects/{kind}:{pid(kind, row_id)}/sheet")
        assert response.status_code == 200, kind
        body = response.get_json()
        rows = _rows(body)
        assert body["kind"] == kind
        ordered = [spec.key for key in sheet.GROUPS for spec in sheet.SPECS[kind] if spec.group == key]
        assert ordered == list(rows), "every spec is a row, grouped in the order of GROUPS"
        assert [group["key"] for group in body["groups"]] == [key for key in sheet.GROUPS
                                                              if any(s.group == key for s in sheet.SPECS[kind])]
        mass = rows["mass_kg"]
        assert mass["stored"] and mass["value"][0]["qty"]["kind"] == "mass" and mass["value"][0]["qty"]["value"] > 0
        assert mass["text"] and mass["text"] != sheet.NOT_STORED
        assert all(row["text"] for row in rows.values())


def test_a_property_with_no_stored_value_says_so(client, mysql_config):
    system_id = _save_wide_binary_with_moons(mysql_config)
    planet_id = _ids(mysql_config, "SELECT id FROM planets WHERE star_system_id = ? ORDER BY id", (system_id,))[0]
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("UPDATE planets SET ocean_depth_km = NULL WHERE id = ?", (planet_id,))
        conn.commit()
    finally:
        conn.close()
    row = _rows(client.get(f"/api/objects/planet:{pid('planet', planet_id)}/sheet").get_json())["ocean_depth_km"]
    assert row == {"key": "ocean_depth_km", "label": "Ocean depth", "value": [sheet.NOT_STORED], "text": sheet.NOT_STORED,
                   "stored": False}


def test_a_kind_with_no_sheet_and_a_missing_object_are_404s(client, mysql_config):
    system_id = _save_wide_binary_with_moons(mysql_config)
    for ref in (f"system:{pid('system', system_id)}", "moon:FFFFFFFFFF-0000005-099", "ship:1", "planet"):
        response = client.get(f"/api/objects/{ref}/sheet")
        assert response.status_code == 404, ref
        assert isinstance(response.get_json()["error"], str)


def test_every_phenomenon_has_a_sheet_over_its_stored_row(client, mysql_config):
    """UX.175: the eight standalone phenomenon kinds serve the same sheet as a body, row for row from `SPECS`."""
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.phenomena.asteroid_field import AsteroidField
    from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
    from planetgen.generation.phenomena.nebula import Nebula
    from planetgen.generation.phenomena.quasar import Quasar
    from planetgen.generation.phenomena.rogue import InterstellarComet, RoguePlanet
    from planetgen.generation.phenomena.supernova_remnant import SupernovaRemnant

    config = SystemConfig()
    made = {"nebula": (Nebula, "nebula"), "asteroid_field": (AsteroidField, "asteroid-field"),
            "black_hole": (BlackHole, "black-hole"), "neutron_star": (NeutronStar, "neutron-star"),
            "supernova_remnant": (SupernovaRemnant, "supernova-remnant"), "rogue_planet": (RoguePlanet, "rogue-planet"),
            "interstellar_comet": (InterstellarComet, "comet"), "quasar": (Quasar, "quasar")}
    for kind, (cls, stored_as) in made.items():
        row_id = store.save_phenomenon(cls(config), config, stored_as, config=mysql_config)
        response = client.get(f"/api/objects/{kind}:{pid(kind, row_id)}/sheet")
        assert response.status_code == 200, kind
        body = response.get_json()
        assert body["kind"] == kind
        rows = _rows(body)
        assert list(rows) == [spec.key for key in sheet.GROUPS for spec in sheet.SPECS[kind] if spec.group == key]
        assert all(row["text"] for row in rows.values()), kind
        assert any(row["stored"] for row in rows.values()), kind
