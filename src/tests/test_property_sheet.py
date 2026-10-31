# tests/test_property_sheet.py

"""
UX.168: the property sheet registry (`describe/sheet.py`) lists every stored column of a star, planet, moon, belt and
comet, so a column a later change adds cannot be forgotten, and `GET /api/objects/<ref>/sheet` serves it.
"""

import pytest

from planetgen.db import column_units, models
from planetgen.describe import model, sheet
from planetgen.describe.quantity import inline_text
from planetgen.physics import unit_registry

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
        if spec.quantity == sheet.AUTO or spec.quantity not in (sheet.TEXT, sheet.FLAG, sheet.TIER, sheet.LIMIT):
            assert (sheet.TABLES[kind], spec.stored_column()) in column_units.COLUMNS, f"{kind}.{spec.key} has no unit"


def test_the_kinds_are_the_bodies_phenomena_and_facilities_a_page_exists_for():
    from planetgen.galaxy import objectref

    assert sorted(sheet.SPECS) == sorted(objectref.BODY_KINDS + objectref.PHENOMENON_KINDS + objectref.FACILITY_KINDS
                                      + objectref.CONTAINER_KINDS) \
        == sorted(sheet.TABLES)


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


def _length_rows():
    """Every sheet row that shows a length on the ladder: `(kind, spec, table, unit)`."""
    rows = []
    for kind, specs in sorted(sheet.SPECS.items()):
        table = sheet.TABLES[kind]
        for spec in specs:
            if spec.quantity in (sheet.TEXT, sheet.FLAG, sheet.RATING, sheet.TIER, sheet.LIMIT, "radius"):
                continue
            unit = column_units.COLUMNS[(table, spec.stored_column())]
            if column_units.unit_info(unit).kind == "length" and unit in unit_registry.UNITS:
                rows.append((kind, spec, table, unit))
    return rows


def test_a_thin_ice_shell_reads_in_meters_not_a_fraction_of_a_kilometer():
    """GEN.769: Boss saw an ice shell of 0.00109 km; the sheet writes 1.09 m, and 3,500 km stays in km."""
    spec = next(spec for spec in sheet.SPECS["planet"] if spec.key == "ice_shell_km")
    assert inline_text(sheet._value(spec, "planets", 0.00109)) == "1.09 m"
    assert inline_text(sheet._value(spec, "planets", 2.4e-7)) == "0.24 mm"
    assert inline_text(sheet._value(spec, "planets", 3500.0)) == "3,500 km"


@pytest.mark.parametrize("kind, spec, table, unit", _length_rows(), ids=lambda value: getattr(value, "key", None))
def test_no_length_row_of_the_sheet_shows_a_number_under_one_in_its_unit(kind, spec, table, unit):
    """GEN.769: whatever size a length column holds, the unit is chosen so the number is at least 1 (the smallest unit,
    the millimeter, excepted: below it the number is scientific)."""
    factor = unit_registry.unit(unit).factor
    for exponent in range(-9, 13):
        if 10.0 ** exponent * factor >= 1e9:    # from a million km the ladder moves on to the astronomical units
            break
        text = inline_text(sheet._value(spec, table, 10.0 ** exponent))
        number = float(text.split(" ")[0].replace(",", "").split("×")[0])
        assert number >= 1 or text.endswith(" mm") or text.endswith(" in"), f"{kind}.{spec.key} {10.0 ** exponent} {unit}: {text}"


def test_the_day_and_night_swing_is_a_temperature_difference_not_a_temperature():
    """GEN.773: Boss saw a swing of 9.07 K converted to -264 °C, the 273.15 offset applied to a difference. A difference
    in kelvin is the same size in °C and 1.8 times it in °F."""
    for kind, table in (("planet", "planets"), ("moon", "moons")):
        spec = next(spec for spec in sheet.SPECS[kind] if spec.key == "temp_mean_swing_k")
        assert column_units.COLUMNS[(table, "temp_mean_swing_k")] == "dK"
        item = sheet._value(spec, table, 9.07)[0]["qty"]
        assert item["kind"] == "temperature difference"
        assert inline_text(sheet._value(spec, table, 9.07)) == "9.07 K (9.07 °C, 16.3 °F)"
        assert "-264" not in inline_text(sheet._value(spec, table, 9.07))
    assert "(hottest to coldest" in next(spec for spec in sheet.SPECS["planet"] if spec.key == "temp_mean_swing_k").label


def test_every_temperature_column_is_a_place_on_the_scale_unless_it_is_declared_a_difference():
    """GEN.773: a column whose name says swing, range, delta, spread or variation and whose unit is kelvin is a
    temperature difference (`dK`), never an absolute temperature."""
    wrong = [key for key, unit in column_units.COLUMNS.items()
             if unit == "K" and any(word in key[1] for word in ("swing", "range", "delta", "spread", "variation", "amplitude"))]
    assert not wrong, f"declare these in dK, not K: {wrong}"


@pytest.mark.parametrize("key, stored, text", [
    ("water_mass_fraction", 0.0004, "400 ppm"), ("water_mass_fraction", 0.0123, "1.23 %"),
    ("ocean_fraction", 0.71, "71 %"), ("land_fraction", 0.29, "29 %"),
    ("water_activity", 0.71, "0.71 of 1 (life needs about 0.6 or more)"), ("ocean_ph", 8.1, "8.1 pH (7 is neutral)"),
])
def test_the_water_rows_state_their_unit(key, stored, text):
    """GEN.770: Boss, 9:02 PM EDT 2026-10-24: water activity, surface under ocean, surface above sea level and the water
    share of the mass had no unit. A share is a percent from 1 percent and parts per million below, the activity says
    its scale, the pH says where neutral is."""
    for kind, table in (("planet", "planets"), ("moon", "moons")):
        spec = next(spec for spec in sheet.SPECS[kind] if spec.key == key)
        assert inline_text(sheet._value(spec, table, stored)) == text





def _numeric_rows():
    """Every sheet row that shows a number: `(kind, spec, table, unit)`."""
    rows = []
    for kind, specs in sorted(sheet.SPECS.items()):
        table = sheet.TABLES[kind]
        for spec in specs:
            if spec.quantity in (sheet.TEXT, sheet.FLAG, sheet.RATING, sheet.TIER, sheet.LIMIT):
                continue
            rows.append((kind, spec, table, column_units.COLUMNS[(table, spec.stored_column())]))
    return rows


def test_no_sheet_row_is_declared_in_a_unit_that_says_nothing():
    """GEN.771: the audit of 2026-10-24 listed 43 ratio and 7 index rows that printed a bare number. None is left, and
    `_plain_number` refuses one, so a column a later change adds cannot print bare."""
    bare = [(kind, spec.key, unit) for kind, spec, _table, unit in _numeric_rows() if unit in sheet.BARE_UNITS]
    assert not bare, f"declare these in a unit that says what they are, in db/column_units.py: {bare}"
    with pytest.raises(ValueError, match="bare number"):
        sheet._plain_number(0.5, "ratio")


@pytest.mark.parametrize("kind, spec, table, unit", _numeric_rows(), ids=lambda value: getattr(value, "key", None))
def test_every_number_of_the_sheet_has_a_unit_a_scale_or_a_noun_beside_it(kind, spec, table, unit):
    """GEN.771, Boss 9:03 PM EDT 2026-10-24: never a number without units or context."""
    for stored in (0.0004, 0.5, 3):
        text = inline_text(sheet._value(spec, table, stored))
        assert not model.is_bare_number(text), f"{kind}.{spec.key} prints {text!r}, a number alone"


@pytest.mark.parametrize("key, stored, text", [
    ("phi4", 0.82, "0.82 of 1"), ("l_rad", 0.9, "0.9 of 1"), ("life_flora", 3, "3 of 9"),
    ("tau_haze", 0.35, "0.35 optical depth"), ("cloud_cover", 0.6, "60 %"), ("bond_albedo", 0.3, "30 %"),
    ("rock_mg_si_ratio", 1.1, "1.1 Mg atoms per Si atom"), ("al26_ratio", 5e-5, "5.00 × 10⁻⁵ ²⁶Al atoms per ²⁷Al atom"),
    ("dose_helio_mult", 1.8, "1.8× the open-space cosmic-ray dose"), ("greenhouse_factor", 0.55, "0.55 (0 is no greenhouse)"),
    ("urey_ratio", 0.8, "80 %"), ("metal_mass_fraction", 0.33, "33 %"), ("envelope_z_fraction", 0.003, "3,000 ppm"),
])
def test_the_ratio_and_index_rows_read_as_decided(key, stored, text):
    spec = next(spec for spec in sheet.SPECS["planet"] if spec.key == key)
    assert inline_text(sheet._value(spec, "planets", stored)) == text


def test_the_swing_reads_as_a_range_and_the_mean_temperature_says_whether_it_has_one():
    """GEN.795, GEN.796, GEN.848, GEN.1108: the mean swing is plus or minus half with the whole in all; the mean temperature shows its coldest and hottest hour with the
    stored value in the middle; on a locked world, which has no day, it says that. The places are the bands, in the table of the bands."""
    row = {"temp_mean_swing_k": 6.0, "surface_temperature_k": 288.0}
    by_key = {spec.key: spec for spec in sheet.SPECS["planet"]}
    text = {key: inline_text(sheet._range_value(by_key[key], "planets", row[key], row)) for key in row}
    assert text["temp_mean_swing_k"] == "\u00b13 K, 6 K in all from the coldest hour to the hottest"
    assert text["surface_temperature_k"].endswith(", from 285 K to 291 K over a day")
    locked = {"temp_mean_swing_k": 0.0, "surface_temperature_k": 288.0}
    assert inline_text(sheet._range_value(by_key["surface_temperature_k"], "planets", 288.0, locked)).endswith(f"({sheet.NO_DAY})")
    assert inline_text(sheet._range_value(by_key["surface_temperature_k"], "planets", 288.0, {})) == inline_text(sheet._value(by_key["surface_temperature_k"], "planets", 288.0))
    assert inline_text(sheet._range_value(by_key["surface_temperature_k"], "planets", None, row)) == sheet.NOT_STORED


def test_the_table_of_the_bands_lists_each_band_with_its_figures():
    """GEN.1108: seven bands side by side, a figure a band does not have (a locked world has no hours) an empty cell, a body with no bands none."""
    from planetgen.physics import latitude_bands, regional_temperature
    earth = regional_temperature.compute(dict(surface_temperature_k=288.0, atmospheric_pressure_pa=101325.0, gravity_g=1.0, rotation_period_hours=24.0, axial_tilt_deg=23.4,
                                              period_years=1.0))
    table = sheet.bands_table(earth)
    assert not table["locked"] and [band["key"] for band in table["bands"]] == list(latitude_bands.NAMES)
    assert [row["key"] for row in table["rows"]] == [field for field, _label in sheet.BAND_ROW_LABELS]
    assert all(len(row["cells"]) == 7 and all(row["cells"]) for row in table["rows"])
    moon = regional_temperature.compute(dict(surface_temperature_k=250.0, atmospheric_pressure_pa=0.0, gravity_g=0.165, rotation_period_hours=655.7, axial_tilt_deg=1.5,
                                             period_years=0.0748))
    rings = sheet.bands_table(moon)
    assert rings["locked"] and [band["key"] for band in rings["bands"]] == list(latitude_bands.RING_NAMES)
    assert [row["key"] for row in rings["rows"]] == ["mean_k", "swing_k", "state"]
    assert sheet.bands_table({}) is None and sheet.bands_table(None) is None


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
        shown = [key for key in rows if not (kind != "star" and key.startswith("p_")) and key != sheet.ECLIPSE_KEY and key not in sheet.CHEMISTRY_KEYS]   # UX.205, UX.254: derived, not a spec
        assert ordered == shown, "every spec is a row, grouped in the order of GROUPS"
        assert all(rows[key]["value"][0]["qty"]["kind"] == "pressure" for key in rows if kind != "star" and key.startswith("p_"))
        assert [group["key"] for group in body["groups"]] == [key for key in sheet.GROUPS
                                                              if any(s.group == key for s in sheet.SPECS[kind]) or (kind, key) == ("star", "composition")]    # UX.254: the star's derived chemistry
        mass = rows["mass_kg"]
        assert mass["stored"] and mass["value"][0]["qty"]["kind"] == "mass" and mass["value"][0]["qty"]["value"] > 0
        assert mass["text"] and mass["text"] != sheet.NOT_STORED
        assert all(row["text"] for row in rows.values())


def test_the_air_of_a_planet_is_one_row_per_gas_the_most_abundant_first(client, mysql_config):
    system_id = _save_wide_binary_with_moons(mysql_config)
    planet_id = _ids(mysql_config, "SELECT id FROM planets WHERE star_system_id = ? ORDER BY id", (system_id,))[0]
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("DELETE FROM planet_gases WHERE planet_id = ?", (planet_id,))
        conn.executemany("INSERT INTO planet_gases (planet_id, gas, partial_pa) VALUES (?, ?, ?)",
                         [(planet_id, "O2", 21000.0), (planet_id, "N2", 78000.0)])
        conn.commit()
    finally:
        conn.close()
    rows = _rows(client.get(f"/api/objects/planet:{pid('planet', planet_id)}/sheet").get_json())
    air = [key for key in rows if key.startswith("p_")]
    assert air == ["p_N2", "p_O2"]
    assert rows["p_N2"]["label"] == "Nitrogen partial pressure" and rows["p_N2"]["value"][0]["qty"]["kind"] == "pressure"


def test_a_moon_in_its_planets_shadow_has_a_row_for_the_longest_eclipse(client, mysql_config):
    """UX.205: the longest eclipse is worked out from the moon's orbit and its planet's radius when the sheet is made, in the rotation
    group; a moon never in the shadow, and every planet, has no such row."""
    from planetgen.physics import daylight

    system_id = _save_wide_binary_with_moons(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        moon = conn.execute("SELECT m.id, m.parent_id FROM moons m JOIN planets p ON p.id = m.parent_id AND m.parent_type = 'planet' "
                            "WHERE m.star_system_id = ? ORDER BY m.id", (system_id,)).fetchone()
        conn.execute("UPDATE planets SET radius_km = 71492.0 WHERE id = ?", (moon["parent_id"],))
        conn.execute("UPDATE moons SET distance_km = 421700.0, period_years = ?, parent_shadow_share = 0.05 WHERE id = ?",
                     (42.46 / unit_registry.converter("yr", "h")(1.0), moon["id"]))
        conn.commit()
    finally:
        conn.close()
    ref = f"moon:{pid('moon', moon['id'])}"
    rows = _rows(client.get(f"/api/objects/{ref}/sheet").get_json())
    eclipse = rows[sheet.ECLIPSE_KEY]
    assert eclipse["stored"] and eclipse["raw"] == pytest.approx(daylight.longest_eclipse_hours(421700.0, 71492.0, 42.46), rel=1e-6)
    assert eclipse["value"][0]["qty"]["kind"] == "time" and eclipse["label"] == "Longest eclipse in the planet's shadow"
    body = client.get(f"/api/objects/{ref}/sheet").get_json()
    assert sheet.ECLIPSE_KEY in [row["key"] for row in next(g for g in body["groups"] if g["key"] == "rotation")["rows"]]
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("UPDATE moons SET parent_shadow_share = 0 WHERE id = ?", (moon["id"],))
        conn.commit()
    finally:
        conn.close()
    assert sheet.ECLIPSE_KEY not in _rows(client.get(f"/api/objects/{ref}/sheet").get_json())
    assert sheet.ECLIPSE_KEY not in _rows(client.get(f"/api/objects/planet:{pid('planet', moon['parent_id'])}/sheet").get_json())


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
                   "raw": None, "stored": False, "absent": True, "source": None, "show_absent": False}


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


def test_the_glance_and_show_when_absent_lists_name_rows_the_kind_has():
    for kind, keys in sheet.GLANCE.items():
        have = {spec.key for spec in sheet.SPECS[kind]}
        assert set(keys) <= have, f"{kind}: {sorted(set(keys) - have)}"
    for kind, keys in sheet.SHOW_WHEN_ABSENT.items():
        have = {spec.key for spec in sheet.SPECS[kind]}
        assert set(keys) <= have, f"{kind}: {sorted(set(keys) - have)}"


def test_a_star_that_varies_by_its_phase_has_a_variability_row_and_one_that_does_not_has_none(client, mysql_config):
    """GEN.944: the kind of variable is worked out from the stored state, class, temperature and mass when the sheet is made."""
    system_id = _save_wide_binary_with_moons(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        star = conn.execute("SELECT id FROM stars WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchone()
        conn.execute("UPDATE stars SET state = 'cepheid' WHERE id = ?", (star["id"],))
        conn.commit()
    finally:
        conn.close()
    ref = f"star:{pid('star', star['id'])}"
    rows = _rows(client.get(f"/api/objects/{ref}/sheet").get_json())
    assert rows[sheet.VARIABILITY_KEY]["text"] == "classical Cepheid" and rows[sheet.VARIABILITY_KEY]["label"] == "Variability"
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("UPDATE stars SET state = 'main-sequence', yerkes_class = 'V' WHERE id = ?", (star["id"],))
        conn.commit()
    finally:
        conn.close()
    assert sheet.VARIABILITY_KEY not in _rows(client.get(f"/api/objects/{ref}/sheet").get_json())


def test_an_asymptotic_giant_that_has_dredged_carbon_up_has_a_carbon_star_row(client, mysql_config):
    """GEN.1026: the carbon or S star row is worked out from the stored state, initial mass, mass, [Fe/H] and age when the sheet is made."""
    from planetgen.physics import constants, post_main_sequence
    window = post_main_sequence.timeline(2.5, 0.0).window(post_main_sequence.ASYMPTOTIC_GIANT)
    age_gy = window.start_gy + 0.9 * (window.end_gy - window.start_gy)
    system_id = _save_wide_binary_with_moons(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        star = conn.execute("SELECT id FROM stars WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchone()
        conn.execute("UPDATE stars SET state = 'asymptotic-giant', yerkes_class = 'III', lifespan_gy = ?, mass_kg = ?, metallicity_fe_h = 0.0, "
                     "age_gy = ? WHERE id = ?", (window.end_gy, 1.2 * constants.SOLAR_MASS_TO_KG, age_gy, star["id"]))
        conn.commit()
    finally:
        conn.close()
    ref = f"star:{pid('star', star['id'])}"
    rows = _rows(client.get(f"/api/objects/{ref}/sheet").get_json())
    assert rows[sheet.CARBON_KEY]["text"] == "carbon star" and rows[sheet.CARBON_KEY]["label"] == "Carbon or S star"
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("UPDATE stars SET state = 'red-giant' WHERE id = ?", (star["id"],))
        conn.commit()
    finally:
        conn.close()
    assert sheet.CARBON_KEY not in _rows(client.get(f"/api/objects/{ref}/sheet").get_json())


def _set_star(mysql_config, star_id, **columns):
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("UPDATE stars SET " + ", ".join(f"{name} = ?" for name in columns) + " WHERE id = ?", (*columns.values(), star_id))
        conn.commit()
    finally:
        conn.close()


def test_a_star_shows_its_chemistry_from_the_stored_values_with_units_and_sources(client, mysql_config):
    """UX.254: the disc population, [alpha/Fe] and the metal mass fraction come from `stellar_composition.star_makeup` on the stored class, mass, [Fe/H]
    and age (no column), and read the same on every request; a red giant adds the helium its first dredge-up added."""
    from planetgen.physics import constants, stellar_abundances, stellar_composition
    system_id = _save_wide_binary_with_moons(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        star = conn.execute("SELECT id, star_system_id FROM stars WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchone()
        distance_ly = sheet._star_distance_ly(conn, star["star_system_id"])
    finally:
        conn.close()
    ref = f"star:{pid('star', star['id'])}"
    seen = set()
    for fe_h, age_gy, mass_sol in ((0.0, 4.6, 1.0), (-0.8, 9.0, 0.8), (-1.4, 11.0, 0.7), (0.3, 1.0, 1.4)):
        mass_kg = mass_sol * constants.SOLAR_MASS_TO_KG
        _set_star(mysql_config, star["id"], state="main-sequence", yerkes_class="V", mass_kg=mass_kg, metallicity_fe_h=fe_h, age_gy=age_gy)
        rows = _rows(client.get(f"/api/objects/{ref}/sheet").get_json())
        key = stellar_composition.abundance_key(mass=mass_kg, distance_ly=distance_ly)
        makeup = stellar_composition.star_makeup(stellar_composition.MAIN_SEQUENCE, fe_h, age_gy, mass_sol, key)
        assert rows["chemistry_branch"]["text"] == sheet._BRANCH_WORDS[makeup.branch]
        assert rows["chemistry_alpha_fe"]["raw"] == pytest.approx(stellar_abundances.alpha_enhancement(makeup.birth))
        assert rows["chemistry_alpha_fe"]["text"].endswith(" dex")
        assert rows["chemistry_metal_fraction"]["raw"] == pytest.approx(makeup.birth.z)
        assert rows["chemistry_metal_fraction"]["value"][0]["qty"]["unit"] == "frac"
        assert all(rows[name]["source"] for name in ("chemistry_branch", "chemistry_alpha_fe", "chemistry_metal_fraction"))
        assert "chemistry_dredge_up_helium" not in rows
        seen.add(makeup.branch)
    assert seen == {stellar_abundances.BRANCH_THIN, stellar_abundances.BRANCH_HIGH_ALPHA}    # the sample reaches both sequences
    _set_star(mysql_config, star["id"], state="red-giant", yerkes_class="III", mass_kg=1.5 * constants.SOLAR_MASS_TO_KG, metallicity_fe_h=0.0, age_gy=2.5)
    rows = _rows(client.get(f"/api/objects/{ref}/sheet").get_json())
    assert rows["chemistry_dredge_up_helium"]["raw"] == pytest.approx(stellar_abundances.first_dredge_up(1.5, rows["chemistry_metal_fraction"]["raw"]).delta_y)
    assert rows["chemistry_dredge_up_helium"]["source"]


def test_a_star_with_no_metallicity_shows_no_chemistry(client, mysql_config):
    system_id = _save_wide_binary_with_moons(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        star = conn.execute("SELECT id FROM stars WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchone()
    finally:
        conn.close()
    _set_star(mysql_config, star["id"], metallicity_fe_h=None)
    rows = _rows(client.get(f"/api/objects/star:{pid('star', star['id'])}/sheet").get_json())
    assert not set(sheet.CHEMISTRY_KEYS) & set(rows)


def _peculiar_star_values(kind):
    from planetgen.physics import chemically_peculiar
    radius_km = 1.8 * chemically_peculiar.SOLAR_RADIUS_KM
    for step in range(2000):
        mass_kg, age_gy = 3.0e30 * (1.0 + step * 1.0e-3), 0.4
        found = chemically_peculiar.of_star("V", mass_kg, radius_km, 8500.0, age_gy)
        if (found.kind if found else None) == kind:
            return mass_kg, radius_km, age_gy, found
    raise AssertionError(f"no {kind} star in 2000 steps")


def _make_a_dwarf(mysql_config, mass_kg, radius_km, age_gy):
    system_id = _save_wide_binary_with_moons(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        star = conn.execute("SELECT id FROM stars WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchone()
        conn.execute("UPDATE stars SET state = 'main-sequence', yerkes_class = 'V', mass_kg = ?, radius_km = ?, temperature_k = 8500.0, age_gy = ? WHERE id = ?",
                     (mass_kg, radius_km, age_gy, star["id"]))
        conn.commit()
    finally:
        conn.close()
    return f"star:{pid('star', star['id'])}"


def test_an_ap_dwarf_has_a_kind_row_and_a_polar_field_row(client, mysql_config):
    """GEN.1027: the Ap or Am kind and the Ap star's polar field are worked out from the stored class, mass, radius, temperature and age when the sheet is made."""
    mass_kg, radius_km, age_gy, found = _peculiar_star_values("Ap")
    ref = _make_a_dwarf(mysql_config, mass_kg, radius_km, age_gy)
    rows = _rows(client.get(f"/api/objects/{ref}/sheet").get_json())
    assert rows[sheet.PECULIAR_KEY]["text"] == "Ap" and rows[sheet.PECULIAR_KEY]["label"] == "Chemically peculiar star"
    assert rows[sheet.PECULIAR_FIELD_KEY]["raw"] == pytest.approx(found.field_gauss) and 300.0 <= found.field_gauss <= 3.0e4


def test_an_am_dwarf_has_a_kind_row_and_no_field_row(client, mysql_config):
    mass_kg, radius_km, age_gy, found = _peculiar_star_values("Am")
    ref = _make_a_dwarf(mysql_config, mass_kg, radius_km, age_gy)
    rows = _rows(client.get(f"/api/objects/{ref}/sheet").get_json())
    assert rows[sheet.PECULIAR_KEY]["text"] == "Am" and sheet.PECULIAR_FIELD_KEY not in rows


def test_a_star_that_is_not_a_peculiar_dwarf_has_no_peculiar_row(client, mysql_config):
    from planetgen.physics import chemically_peculiar
    radius_km = 1.8 * chemically_peculiar.SOLAR_RADIUS_KM
    mass_kg = next(3.0e30 * (1.0 + k * 1.0e-3) for k in range(2000)
                   if chemically_peculiar.of_star("V", 3.0e30 * (1.0 + k * 1.0e-3), radius_km, 8500.0, 0.4) is None)
    ref = _make_a_dwarf(mysql_config, mass_kg, radius_km, 0.4)
    rows = _rows(client.get(f"/api/objects/{ref}/sheet").get_json())
    assert sheet.PECULIAR_KEY not in rows and sheet.PECULIAR_FIELD_KEY not in rows
