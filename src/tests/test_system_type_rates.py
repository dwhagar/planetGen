# tests/test_system_type_rates.py

"""
GEN.268: `SYSTEM_TYPE_RATES` (`planetgen/generation/phenomena/system_types.py`)
has one entry per row of the catalogue of system types in
docs/design/multi-star-systems.md (sections 13.1 to 13.6). These tests fail
when a row of the note has no entry, an entry has no generator, owner or
audit, or a number is out of order.
"""

import importlib
import importlib.util
import pathlib
import re

import pytest

from planetgen.generation.phenomena import system_types as types

_NOTE = pathlib.Path(types.__file__).resolve().parents[4] / "docs" / "design" / "multi-star-systems.md"
_TAGS = {"S", "C", "R", "I"}


def _catalogue_rows():
    """`{section: [first cell of each table row]}` for 13.1 to 13.6 of the note."""
    rows, section = {}, None
    for line in _NOTE.read_text(encoding="utf-8").splitlines():
        heading = re.match(r"### (13\.[1-6]) ", line)
        if heading:
            section = heading.group(1)
            rows[section] = []
        elif line.startswith("## "):
            section = None
        elif section and line.startswith("| ") and not line.startswith(("| Type", "|---")):
            rows[section].append(line.split("|")[1].strip())
    return rows


def test_every_row_of_the_catalogue_has_an_entry_and_no_entry_is_extra():
    note = _catalogue_rows()
    assert sum(len(cells) for cells in note.values()) == 57
    for section, cells in note.items():
        entries = [entry.row for entry in types.SYSTEM_TYPE_RATES.values() if entry.section == section]
        assert sorted(entries) == sorted(cells), section
    assert len({entry.row for entry in types.SYSTEM_TYPE_RATES.values()}) == len(types.SYSTEM_TYPE_RATES)


@pytest.mark.parametrize("key", list(types.SYSTEM_TYPE_RATES))
def test_an_entry_is_complete_and_in_order(key):
    entry = types.SYSTEM_TYPE_RATES[key]
    assert entry.unit in types.UNITS and entry.mode in types.MODES
    assert entry.reference and entry.note is not None
    assert entry.populations and set(entry.populations) <= set(types.POPULATIONS)
    assert set(entry.tag.split(",")) <= _TAGS
    assert entry.audit and set(entry.audit) <= set(types.AUDITS) and 10 in entry.audit
    assert entry.owner is None or (isinstance(entry.owner, int) and 1 <= entry.owner < 1000)
    if entry.default is None:
        assert entry.derived, f"{key}: no default and no statement of where the number comes from"
    else:
        assert entry.default >= 0
        if entry.unit == "fraction":
            assert entry.default <= 1
    if entry.default is not None and entry.low is not None:
        assert entry.low <= entry.default
    if entry.default is not None and entry.high is not None:
        assert entry.default <= entry.high
    for value in (entry.by_class or {}).values():
        assert isinstance(value, (int, float)) and value >= 0


def test_every_generator_is_named_and_exists_once_its_module_does():
    module_exists = importlib.util.find_spec(types.GENERATOR_MODULE.rpartition(".")[0]) is not None \
        and importlib.util.find_spec(types.GENERATOR_MODULE) is not None
    for key, entry in types.SYSTEM_TYPE_RATES.items():
        assert entry.generator, key
        if ":" in entry.generator:
            assert entry.owner is None, key
            module, _, attr = entry.generator.partition(":")
            assert hasattr(importlib.import_module(module), attr), f"{key}: {entry.generator}"
        else:
            assert entry.generator.isidentifier() and entry.owner is not None, key
            if module_exists:
                assert hasattr(importlib.import_module(types.GENERATOR_MODULE), entry.generator), key


def test_the_audits_the_note_names_are_assigned():
    rates = types.SYSTEM_TYPE_RATES
    assert {k for k, v in rates.items() if 12 in v.audit} == {"quintuple-to-septuple", "eight-to-twelve-stars"}
    assert {k for k, v in rates.items() if 14 in v.audit} == {"young-group", "decayed-trapezium-remnant"}
    for key in ("contact-binary", "algol", "cataclysmic-variable", "am-cvn", "high-mass-xray-binary"):
        assert 13 in rates[key].audit, key


def test_the_defaults_are_the_notes_figures():
    rates = types.SYSTEM_TYPE_RATES
    assert sum(rates[key].default for key in ("close-binary", "intermediate-binary", "wide-binary",
                                              "ultra-wide-pair")) == pytest.approx(1.0)
    assert rates["contact-binary"].default == 0.002
    assert rates["close-double-white-dwarf"].default == 0.095
    assert rates["brown-dwarf-wide"].default == 0.032
    assert rates["thorne-zytkow"].theoretical and rates["thorne-zytkow"].default == 100
    assert types.fraction_of_systems("thorne-zytkow") == pytest.approx(3.4e-10, rel=0.01)
    assert rates["runaway-star"].low == 0.10 and rates["runaway-star"].high == 0.30
    assert {k for k, v in rates.items() if v.theoretical} >= {
        "thorne-zytkow", "magnetar-binary", "planetary-mass-companion", "compact-triple"}


def test_the_n_tail_follows_the_notes_ratios():
    tail = types.SYSTEM_TYPE_RATES["eight-to-twelve-stars"].by_class
    counts = [tail[f"N{n}"] * types.GALAXY_SYSTEMS for n in range(8, 13)]
    assert counts == pytest.approx([1.1e6, 1.1e5, 1.1e4, 1.1e3, 1.1e2])
    quint = types.SYSTEM_TYPE_RATES["quintuple-to-septuple"].by_class
    assert quint == {"N5": 2.4e-3, "N6": 2.9e-4, "N7": 3.0e-5}


def test_the_tuning_knobs_are_the_notes_defaults():
    knobs = types.TUNING_KNOBS
    assert knobs["tight_tail_weight"] == 1.0 and knobs["post_ce_fraction"] == 0.10
    assert knobs["merger_product_fraction"] == {"massive": 0.08, "field": 0.005, "white-dwarf": 0.10}
    assert min(knobs["accreting_fraction"].values()) == 0.01 and max(knobs["accreting_fraction"].values()) == 0.40


def test_a_type_without_a_number_of_its_own_says_so():
    with pytest.raises(ValueError, match="no default of its own"):
        types.fraction_of_systems("hot-subdwarf-stripped")
