"""
TEST.215: the Solar System calibration targets live in test files only, and no generator table reads a measured Solar System shape unless its row is
a flagged stand-in with a replacement item (Boss, contradiction sweep fork 8, "Flagged stand-ins only"; docs/design/composition/body-composition-research.md 14.4, R5).

1. A `calibration_only` row cannot be in a registered table: `Registry.register` refuses it. The calibration targets (Earth, Mars, the Moon, the
   Galilean moons, Sirius B) are rows here, in the test file.
2. Every registered row whose source or limit text names a Solar System body is classified in `tests/data/solar_system_rows.json`: a flagged stand-in
   (checked: the state is `stand-in` and a replacement item is named), a row that is `pending_decision` (it rests on a Solar System body and is
   not flagged; reported to the TODO thread, which decides whether to flag or re-source it), or `reviewed_not_solar_system` with the reason. A new row that names a body
   and is in none of them fails, and so does an entry that no longer matches its row.
"""

import json
import pathlib
import re

import pytest

from planetgen.generation import possibility_table as pt
from tests.audit import sweep_tables

DATA = json.loads((pathlib.Path(__file__).resolve().parent / "data" / "solar_system_rows.json").read_text(encoding="utf-8"))

BODY = re.compile(r"solar system|solar-system|\bearth\b|\bmars\b|martian|\bmoon\b|lunar|\bvenus\b|\btitan\b|\bpluto\b|\btriton\b|\bjupiter\b|\bsaturn\b|europa|ganymede|"
                  r"callisto|\bio\b|halley|\b67p\b|meteorite|chondrite|\bsirius\b|\bceres\b|\bvesta\b|bennu|ryugu|\bneptune\b|\buranus\b|enceladus|\bsun\b|\bsolar\b", re.I)
UNITS = re.compile(r"(earth|jupiter|solar)[ -](mass|masses|radius|radii|luminosity|metallicity|abundance|composition|value|units|oxides)|earth.mass|\bMe\b", re.I)
"""The scan: a body's name in a row's source or limit, once the units that use the name (Earth masses, Jupiter radii) are taken out."""

CALIBRATION_TARGETS = (
    ("earth", "core mass fraction 0.325 (PREM)"),
    ("mars", "core mass fraction about 0.24 (InSight)"),
    ("moon", "core mass fraction about 0.02 (Apollo seismics, GRAIL)"),
    ("io", "mean density 3528 kg/m3"),
    ("europa", "mean density 3013 kg/m3"),
    ("ganymede", "mean density 1942 kg/m3"),
    ("callisto", "mean density 1834 kg/m3"),
    ("sirius-b", "1.018 solar masses at 5,840 km [R0454]"),
)
"""tuple: The measured bodies' target rows, as a test would hold them (a test point, never a table source)."""


def _rows():
    registry = sweep_tables.load_tables()
    return {(table.name, row.name): row for table in registry.tables() for row in table.rows}


def _listed(section):
    return {(table, row): why for table, rows in DATA[section].items() for row, why in rows.items()}


def _hits(rows):
    return {key for key, row in rows.items() if BODY.search(UNITS.sub("", f"{row.source} | {row.limit}"))}


def _calibration_table(name="calibration"):
    return pt.Table(name, 1, [pt.Row(body, (pt.ANY_CLASS,), 1.0, "a measured target for a test", notes=target, calibration_only=True)
                              for body, target in CALIBRATION_TARGETS])


def test_a_table_of_calibration_rows_is_a_test_object_and_the_registry_refuses_it():
    table = _calibration_table()
    assert [row.name for row in table.rows] == [body for body, _ in CALIBRATION_TARGETS]
    registry = pt.Registry()
    with pytest.raises(pt.TableError, match="calibration-only"):
        registry.register(table)
    assert "calibration" not in registry


def test_one_calibration_row_among_the_rows_of_a_table_is_refused_too():
    rows = [pt.Row("a", (pt.ANY_CLASS,), 1.0, "src"), pt.Row("earth", (pt.ANY_CLASS,), 1.0, "src", calibration_only=True)]
    with pytest.raises(pt.TableError, match="earth"):
        pt.Registry().register(pt.Table("mixed", 1, rows))


def test_no_registered_row_is_calibration_only():
    assert [key for key, row in _rows().items() if row.calibration_only] == []


def test_every_row_that_names_a_solar_system_body_is_classified():
    rows = _rows()
    classified = {*_listed("stand_ins"), *_listed("pending_decision"), *_listed("reviewed_not_solar_system")}
    unclassified = sorted(_hits(rows) - classified)
    assert not unclassified, f"rows naming a Solar System body that tests/data/solar_system_rows.json does not classify (stand-in, pending or reviewed): {unclassified}"


def test_a_flagged_stand_in_is_a_stand_in_with_a_replacement_item():
    rows = _rows()
    for key in _listed("stand_ins"):
        assert key in rows, f"{key} is listed as a stand-in and is not a registered row"
        assert rows[key].state == pt.STAND_IN and rows[key].replacement.strip(), f"{key} rests on a Solar System body and is not a flagged stand-in with a replacement item"


def test_the_pending_rows_exist_and_are_not_yet_flagged():
    """A pending row that has been flagged (or removed) leaves the list, so the list is only what still needs a decision."""
    rows = _rows()
    for key in _listed("pending_decision"):
        assert key in rows, f"{key} is pending and is not a registered row: remove it from tests/data/solar_system_rows.json"
        assert rows[key].state != pt.STAND_IN, f"{key} is flagged as a stand-in now: move it to stand_ins"


def test_the_reviewed_rows_still_name_a_body():
    rows = _rows()
    hits = _hits(rows)
    stale = sorted(key for key in _listed("reviewed_not_solar_system") if key not in hits)
    assert not stale, f"reviewed rows that no longer name a Solar System body: {stale}"


def test_no_row_is_classified_twice():
    sections = [set(_listed(name)) for name in ("stand_ins", "pending_decision", "reviewed_not_solar_system")]
    assert not (sections[0] & sections[1] or sections[0] & sections[2] or sections[1] & sections[2])
