# tests/test_sheet_sources.py

"""
UX.187: a physics row of the property sheet carries a one-line source (`describe/sheet_sources.py`): how the value is worked
out and the item that does it. The table is held to the sheet's rows, each line is one sentence, a row with no line says
nothing, and a page and the API show the line beside the value once.
"""

import re

import pytest
from markupsafe import escape

from planetgen.db import store

from planetgen.describe import sheet
from planetgen.describe.sheet_sources import SOURCES, source_of
from tests.publicids import pid
from tests.test_api import _save_wide_binary_with_moons, client  # noqa: F401
from tests.test_property_sheet import _ids, _rows
from tests.test_uncharted_phenomenon import db_app  # noqa: F401

AUDIT_ROWS = {
    "planet": ("eccentricity", "orbital_inclination_deg", "flux_variation", "tidal_flux_w_m2", "internal_heat_flux_w_m2",
               "urey_ratio", "lid_regime", "c_th0_ppm", "c_u0_ppm", "rotation_period_hours", "axial_tilt_deg",
               "magnetic_moment_a_m2", "ring_outer_km", "ring_optical_depth", "life_stage", "life_age_abiogenesis_gyr",
               "life_age_photosynthesis_gyr", "life_age_complex_cells_gyr", "life_age_multicellularity_gyr",
               "life_age_technological_civilization_gyr"),
    "moon": ("eccentricity", "flux_variation", "tidal_flux_w_m2", "internal_heat_flux_w_m2", "urey_ratio", "lid_regime",
             "rotation_period_hours", "axial_tilt_deg", "magnetic_moment_a_m2", "life_stage"),
    "star": ("rotation_period_hours", "axial_tilt_deg"),
}


@pytest.mark.parametrize("kind", sorted(SOURCES))
def test_every_source_belongs_to_a_row_of_the_sheet_and_reads_as_one_sentence(kind):
    keys = {spec.key for spec in sheet.SPECS[kind]} | (set(sheet.CHEMISTRY_KEYS) if kind == "star" else set())    # UX.254: rows worked out from the stored values
    for key, line in SOURCES[kind].items():
        assert key in keys, f"{kind} has no row {key}"
        assert "\n" not in line and line.endswith(".") and 20 < len(line) < 240, (kind, key)
        assert "—" not in line and not re.search(r"\s{2}", line), (kind, key)


@pytest.mark.parametrize("kind", sorted(AUDIT_ROWS))
def test_the_audited_physics_rows_each_have_a_source(kind):
    assert [key for key in AUDIT_ROWS[kind] if source_of(kind, key) is None] == []


def test_a_row_with_no_line_written_has_none_and_a_moons_inclination_is_one():
    assert source_of("planet", "mass_kg") is None and source_of("system", "x") is None
    assert source_of("planet", "orbital_inclination_deg") and source_of("moon", "orbital_inclination_deg") is None


def test_the_sheet_api_carries_the_source_of_a_row_and_none_for_the_rest(client, mysql_config):  # noqa: F811
    system_id = _save_wide_binary_with_moons(mysql_config)
    for kind, table in (("star", "stars"), ("planet", "planets"), ("moon", "moons")):
        row_id = _ids(mysql_config, f"SELECT id FROM {table} WHERE star_system_id = ? ORDER BY id", (system_id,))[0]
        rows = _rows(client.get(f"/api/objects/{kind}:{pid(kind, row_id)}/sheet").get_json())
        assert all("source" in row for row in rows.values())
        assert {key: row["source"] for key, row in rows.items() if row["source"]} == {
            key: SOURCES[kind][key] for key in SOURCES[kind] if key in rows}
        assert rows["mass_kg"]["source"] is None


def _set_stored(mysql_config, table, row_id, **columns):
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            assignments = ", ".join(f"{name} = ?" for name in columns)
            conn.execute(f"UPDATE {table} SET {assignments} WHERE id = ?", (*columns.values(), row_id))
    finally:
        conn.close()


def test_each_physics_row_of_a_page_stands_once_with_its_source_beside_it(db_app, mysql_config):  # noqa: F811
    """UX.187's walk, on the fixture system: with the audited columns stored, the page of a star, a planet and a moon lists each
    row with its one-line source once, in the value cell and nowhere else, and the cell still has its unit."""
    system_id = _save_wide_binary_with_moons(mysql_config)
    physics = {
        "planet": dict(eccentricity=0.05, orbital_inclination_deg=1.2, urey_ratio=0.4, internal_heat_flux_w_m2=0.08,
                       rotation_period_hours=24.0, axial_tilt_deg=23.4, magnetic_moment_a_m2=8e22),
        "moon": dict(eccentricity=0.01, urey_ratio=0.3, rotation_period_hours=90.0, axial_tilt_deg=0.5),
        "star": dict(rotation_period_hours=600.0, axial_tilt_deg=3.0),
    }
    client = db_app.test_client()
    for kind, table in (("star", "stars"), ("planet", "planets"), ("moon", "moons")):
        row_id = _ids(mysql_config, f"SELECT id FROM {table} WHERE star_system_id = ? ORDER BY id", (system_id,))[0]
        _set_stored(mysql_config, table, row_id, **physics[kind])
        page = client.get(f"/object/{kind}:{pid(kind, row_id)}").get_data(as_text=True)
        cells = re.findall(r'<th scope="row">([^<]*)</th><td[^>]*>(.*?)</td>', page, flags=re.S)
        for key in physics[kind]:
            spec = next(spec for spec in sheet.SPECS[kind] if spec.key == key)
            line = source_of(kind, key)
            found = [cell for name, cell in cells if name == spec.label]
            assert found, f"{kind} {spec.label} is not on the page"
            if line:
                marker = f'<small class="sheet-source">{escape(line)}</small>'
                holders = [cell for cell in found if marker in cell]    # the At a glance panel may repeat a row; its copy has no source
                assert len(holders) == 1 and page.count(marker) == 1, f"{kind} {key}"
                found = holders
            else:
                found = found[:1]
            value = found[0].split("<small")[0]
            assert re.search(r"[A-Za-z\u00b0%\u00b1]", re.sub(r"<[^>]+>", "", value)), f"{kind} {key}: {value!r} has no unit"
