# tests/test_sensor.py

"""
UX.275: the Sensor registry (`web/lib/sensor.py`). Every kind of the property sheet has an entry (a painter, or the reason it has none), the
painter ids agree with `static/sensorview.js`, each builder makes its record from a stored row alone, and `GET /api/objects/<ref>/sensor`
serves it.
"""

import json
import re

import pytest

from planetgen.describe import sheet
from planetgen.web.lib import fmt, sensor
from planetgen.web.maps import bodylook, phenomenonrender


def test_every_kind_of_the_property_sheet_has_an_entry_and_no_other_kind_does():
    missing = sorted(set(sheet.SPECS) - set(sensor.REGISTRY))
    assert not missing, (
        f"the property sheet has kinds the Sensor registry does not: {missing}. Add an Entry for each to REGISTRY in planetgen/web/lib/sensor.py, "
        "with its painter and look builder, or with the reason it has no picture yet.")
    assert not set(sensor.REGISTRY) - set(sheet.SPECS), sorted(set(sensor.REGISTRY) - set(sheet.SPECS))
    assert all(entry.kind == key for key, entry in sensor.REGISTRY.items())


def test_an_entry_has_a_painter_and_a_builder_or_a_reason():
    for entry in sensor.REGISTRY.values():
        assert (entry.painter is None) == (entry.build is None) == bool(entry.reason), entry.kind
    with pytest.raises(ValueError):
        sensor.Entry("star", painter="sphere")
    with pytest.raises(ValueError):
        sensor.Entry("star")
    with pytest.raises(ValueError):
        sensor.Entry("star", painter="hologram", build=lambda conn, row: {})
    with pytest.raises(ValueError):
        sensor.entry_for("system")


def test_the_facility_has_no_picture_before_9_3():
    entry = sensor.entry_for("facility")
    assert entry.painter is None and "9.3" in entry.reason


def test_a_kind_with_no_painter_needs_no_row():
    record = sensor.sensor_for(None, "facility", 1)
    assert record["painter"] is None and record["look"] is None and record["caption"] is None and record["reason"] == sensor.FACILITY_REASON
    assert record["kind"] == "facility" and record["id"] == 1 and record["v"] == sensor.SENSOR_VERSION
    with pytest.raises(ValueError):
        sensor.sensor_for(None, "system", 1)


def _script():
    with open(f"{fmt._STATIC_DIR}/sensorview.js", encoding="utf-8") as handle:
        return handle.read()


def test_the_painter_ids_are_the_scripts_and_each_names_a_module_that_exists():
    source = _script()
    block = source[source.index("export const PAINTERS"):]
    ids = re.findall(r"^  (\w+): \{ module: \"([\w.]+\.js)\"", block, re.M)
    assert [name for name, _ in ids] == list(sensor.PAINTERS)
    for _, module in ids:
        with open(f"{fmt._STATIC_DIR}/{module}", encoding="utf-8"):
            pass
    assert f"export const SENSOR_VERSION = {sensor.SENSOR_VERSION};" in source
    assert {entry.painter for entry in sensor.REGISTRY.values() if entry.painter} <= set(sensor.PAINTERS)


# The builders, from a stored row alone ----------------------------------------------------------------------

def test_a_phenomenon_record_is_the_parameters_its_painter_reads():
    rows = {
        "neutron_star": {"id": 3, "spin_period_ms": 33.0, "pulsar_type": "radio", "magnetic_field_gauss": 1e12, "surface_temperature_k": 1e6},
        "black_hole": {"id": 4, "has_accretion_disk": 1, "mass_class": "stellar", "spin": 0.5, "temperature_k": 1e7},
        "quasar": {"id": 5, "is_radio_loud": 1},
        "interstellar_comet": {"id": 6, "is_active": 0},
        "rogue_planet": {"id": 7, "planet_type": "t", "surface_temperature_k": 90.0, "radius_km": 6000.0, "has_internal_heat": 0},
    }
    for kind, row in rows.items():
        record = sensor.REGISTRY[kind].build(None, row)
        assert record["look"]["kind"] == kind and record["caption"] and record["estimate"] == []
        params, caption = phenomenonrender.view_params(kind, row)
        assert record["look"] == {"kind": kind, **params} and record["caption"] == caption
        json.dumps(record)
        assert record == sensor.REGISTRY[kind].build(None, dict(row)), "the same row gives the same record"
    with pytest.raises(ValueError):
        phenomenonrender.view_params("asteroid_field", {})


def test_a_nebula_record_names_the_volume_it_is_drawn_from():
    record = sensor.REGISTRY["nebula"].build(None, {"id": 31, "nebula_type": "emission"})
    look = record["look"]
    assert look["nebulaId"] == "31" and look["volumePath"] == "/galaxy/nebula/31/volume" and look["surroundingsPath"].endswith("/31/surroundings")
    assert look["color"].startswith("#") and 0 <= look["edgeOpacity"] <= look["coreOpacity"] <= 1


def test_a_star_record_is_the_look_of_its_row():
    row = {"id": 12, "star_type": "G", "temperature_k": 5772.0, "luminosity_w": 3.828e26, "radius_km": 695700.0, "mass_kg": 1.989e30,
           "yerkes_class": "V", "rotation_period_hours": 600.0}
    record = sensor.REGISTRY["star"].build(None, row)
    assert record["look"] == bodylook.star_look(row)
    assert record["look"]["seed"] == bodylook.seed_for("star", 12)


# The API, over a real database -----------------------------------------------------------------------------

from planetgen.db import store  # noqa: E402
from tests.publicids import pid  # noqa: E402
from tests.test_api import _save_wide_binary_with_moons, client  # noqa: E402,F401


def _first(mysql_config, table, system_id):
    conn = store.get_connection(mysql_config)
    try:
        return conn.execute(f"SELECT id FROM {table} WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchone()["id"]
    finally:
        conn.close()


@pytest.mark.db
def test_the_sensor_of_a_star_a_planet_and_a_moon_is_the_look_of_its_stored_row(client, mysql_config):
    system_id = _save_wide_binary_with_moons(mysql_config)
    for kind, table in (("star", "stars"), ("planet", "planets"), ("moon", "moons")):
        row_id = _first(mysql_config, table, system_id)
        response = client.get(f"/api/objects/{kind}:{pid(kind, row_id)}/sensor")
        assert response.status_code == 200, kind
        body = response.get_json()
        assert body["v"] == sensor.SENSOR_VERSION and body["kind"] == kind and body["painter"] == "sphere" and body["reason"] is None
        assert body["look"]["kind"] == kind and body["look"]["seed"] == bodylook.seed_for(kind, row_id) and body["caption"]
        assert body["inputs"] and isinstance(body["estimate"], list)
        assert response.get_json() == client.get(f"/api/objects/{kind}:{pid(kind, row_id)}/sensor").get_json(), "from the row alone: the same twice"


@pytest.mark.db
def test_a_kind_with_no_painter_says_why_and_a_missing_object_is_a_404(client, mysql_config):
    system_id = _save_wide_binary_with_moons(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        belt = conn.execute("SELECT id FROM asteroid_belts WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchone()
    finally:
        conn.close()
    if belt is not None:
        body = client.get(f"/api/objects/belt:{pid('belt', belt['id'])}/sensor").get_json()
        assert body["painter"] is None and body["look"] is None and body["reason"]
    assert client.get("/api/objects/planet:FE81000A2B-0000005-999/sensor").status_code == 404
    assert client.get("/api/objects/nonsense/sensor").status_code == 404
