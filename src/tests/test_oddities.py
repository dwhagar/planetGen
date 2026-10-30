"""
GEN.1036: `planetgen.physics.oddities`, the oddity table and its loader (docs/design/systems/weird-worlds.md 4.1 and 4.4). Four groups: the
shipped file loads and its vocabulary is sound; the schema refuses every row that lacks a source, a gate, limits, a tier or a tweak
level (one test per refusal); the gate, the limits and the rate (with the weirdness dial) read as the design says; and the draws are
isolated: a row's roll uses a stream of its own, so adding, removing or turning off a row moves no other draw.
"""

import copy
import json
import pathlib
import random
import sys

import pytest

from planetgen import tuning
from planetgen.physics import oddities
from planetgen.util import draw

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import build_stand_in_values as builder  # noqa: E402

PLAN = (ROOT / "docs" / "design" / "research" / "gap-filling-plan.md").read_text(encoding="utf-8")


def _document():
    with open(oddities.DATA_FILE, encoding="utf-8") as handle:
        return json.load(handle)


def _row(**changes):
    row = {
        "id": "W10", "name": "Iron-rich core", "kind": "planet",
        "requires": [{"quantity": "mass_earth", "lt": 5.0}, {"quantity": "collision_count", "min": 1}],
        "rate": {"value": 0.03, "source": "Plotnykov and Valencia 2020, recalled", "label": "RECALLED", "tier": "T2", "tweak": "L0",
                 "registry": "OD1"},
        "apply": {"lever": "L3", "input": "core_mass_fraction", "effect": "strip_mantle", "params": {"max_fraction": 0.85}},
        "limits": [{"quantity": "core_mass_fraction", "lt": 0.85}],
        "source": "Marcus et al. 2010; Benz et al. 2007",
        "text": "Its mantle was knocked off by a giant impact, so it is mostly iron.",
    }
    row.update(changes)
    return row


def _with(*rows):
    document = _document()
    document["rows"] = list(rows)
    return document


# --- the shipped file -----------------------------------------------------------------------------------------------------------------------

def test_the_shipped_table_loads_and_holds_the_rows_the_modules_have_added():
    rows, levers, quantities = oddities.load()
    assert rows == oddities.ROWS and {row.id for row in rows} == set(oddities.BY_ID)
    assert set(levers) == {"L1", "L2", "L3", "L4", "L5", "L6", "L7"} and quantities


def test_every_lever_has_a_label_and_inputs_and_no_input_belongs_to_two_levers():
    inputs = [name for lever in oddities.LEVERS.values() for name in lever["inputs"]]
    assert all(lever["label"] and lever["inputs"] for lever in oddities.LEVERS.values())
    assert len(inputs) == len(set(inputs))


def test_every_quantity_a_gate_may_read_has_a_label_and_a_note():
    assert all(entry.get("label") and entry.get("note") for entry in oddities.QUANTITIES.values())


def test_every_rate_of_the_shipped_table_has_a_registry_row():
    ids = {row[0] for _number, (_header, rows) in builder.registry_tables(PLAN).items() for row in rows}
    assert all(row.rate.registry in ids for row in oddities.ROWS)


def test_the_dial_and_its_scale_are_registered_defaults():
    assert tuning.WEIRDNESS_DIAL == 1.0 and tuning.WEIRDNESS_RATE_SCALE == {}
    tables = builder.registry_tables(PLAN)
    assert {"OD1", "OD2"} <= {row[0] for _header, rows in tables.values() for row in rows}


# --- the schema refuses ---------------------------------------------------------------------------------------------------------------------

def test_a_good_row_loads_with_everything_it_holds():
    rows, _levers, _quantities = oddities.parse(_with(_row()))
    (row,) = rows
    assert (row.id, row.kind, row.lever, row.input, row.effect) == ("W10", "planet", "L3", "core_mass_fraction", "strip_mantle")
    assert row.params == {"max_fraction": 0.85} and row.rate.value == 0.03 and row.rate.registry == "OD1"
    assert oddities.words(row).startswith("Odd: Iron-rich core. Its mantle")


def _drop(row, *path):
    row = copy.deepcopy(row)
    target = row
    for key in path[:-1]:
        target = target[key]
    del target[path[-1]]
    return row


@pytest.mark.parametrize("bad,match", [
    (_drop(_row(), "source"), "'source'"),
    (_drop(_row(), "text"), "'text'"),
    (_drop(_row(), "name"), "'name'"),
    (_drop(_row(), "requires"), "'requires' needs at least one condition"),
    (_row(requires=[]), "'requires' needs at least one condition"),
    (_drop(_row(), "limits"), "'limits' needs at least one condition"),
    (_row(limits=[]), "'limits' needs at least one condition"),
    (_drop(_row(), "rate"), "'rate' is missing"),
    (_drop(_row(), "rate", "source"), "'source'"),
    (_drop(_row(), "rate", "registry"), "'registry'"),
    (_drop(_row(), "rate", "tier"), "rate.tier"),
    (_drop(_row(), "rate", "tweak"), "rate.tweak"),
    (_row(rate=dict(_row()["rate"], tier="T0")), "rate.tier"),
    (_row(rate=dict(_row()["rate"], tweak="L5")), "rate.tweak"),
    (_row(rate=dict(_row()["rate"], label="HEARD")), "rate.label"),
    (_row(rate=dict(_row()["rate"], value=1.5)), "rate.value"),
    (_row(rate=dict(_row()["rate"], value=-0.1)), "rate.value"),
    (_row(rate=dict(_row()["rate"], value="few")), "rate.value"),
    (_row(rate=dict(_row()["rate"], value=True)), "rate.value"),
    (_row(kind="galaxy"), "kind"),
    (_row(id="weird"), "the id is not"),
    (_row(apply={"lever": "L9", "input": "eccentricity", "effect": "x"}), "apply.lever"),
    (_row(apply={"lever": "L3", "input": "eccentricity", "effect": "strip_mantle"}), "not an input of lever L3"),
    (_row(apply={"lever": "L3", "input": "core_mass_fraction", "effect": "Strip Mantle"}), "apply.effect"),
    (_row(apply={"lever": "L3", "input": "core_mass_fraction", "effect": "strip_mantle", "params": [1]}), "apply.params"),
    (_row(requires=[{"quantity": "luck", "min": 1}]), "does not declare"),
    (_row(limits=[{"quantity": "mass_earth", "above": 1}]), "outside"),
    (_row(requires=[{"quantity": "mass_earth"}]), "is not a quantity with"),
    (_row(requires=["mass_earth"]), "is not a quantity with"),
])
def test_the_loader_refuses_a_row_without_what_the_design_requires(bad, match):
    with pytest.raises(oddities.OddityError, match=match):
        oddities.parse(_with(bad))


def test_the_loader_refuses_two_rows_with_one_id():
    with pytest.raises(oddities.OddityError, match="used by two rows"):
        oddities.parse(_with(_row(), _row(name="Another")))


def test_the_loader_refuses_a_file_without_its_vocabulary():
    for key in ("kinds", "levers", "quantities"):
        document = _document()
        del document[key]
        with pytest.raises(oddities.OddityError):
            oddities.parse(document)
    document = _document()
    document["levers"]["L1"]["inputs"] = []
    with pytest.raises(oddities.OddityError, match="lever L1"):
        oddities.parse(document)


def test_a_refused_file_is_a_data_error():
    assert issubclass(oddities.OddityError, oddities.DataError)


# --- the gate, the limits and the rate ------------------------------------------------------------------------------------------------------

def _table():
    return oddities.parse(_with(_row(), _row(id="W14", name="Disintegrating planet", requires=[{"quantity": "mass_earth", "lt": 0.1}])))[0]


def test_a_body_that_passes_every_condition_of_a_gate_is_a_candidate_and_one_that_misses_any_is_not():
    rows = _table()
    assert [row.id for row in oddities.candidates("planet", {"mass_earth": 2.0, "collision_count": 3}, rows)] == ["W10"]
    assert oddities.candidates("planet", {"mass_earth": 2.0, "collision_count": 0}, rows) == []
    assert oddities.candidates("planet", {"mass_earth": 8.0, "collision_count": 3}, rows) == []
    assert [row.id for row in oddities.candidates("planet", {"mass_earth": 0.05}, rows)] == ["W14"]
    assert oddities.candidates("moon", {"mass_earth": 2.0, "collision_count": 3}, rows) == []


def test_a_quantity_the_body_does_not_store_fails_the_gate_instead_of_raising():
    rows = _table()
    assert oddities.candidates("planet", {"mass_earth": 2.0}, rows) == []


def test_the_limits_name_what_does_not_hold_after_the_effect():
    (row, _other) = _table()
    assert oddities.limit_failures(row, {"core_mass_fraction": 0.7}) == []
    assert oddities.limit_failures(row, {"core_mass_fraction": 0.9}) == [{"quantity": "core_mass_fraction", "lt": 0.85}]
    assert oddities.limit_failures(row, {}) == [{"quantity": "core_mass_fraction", "lt": 0.85}]


def test_the_probability_is_the_rate_times_the_dial_and_the_rows_own_scale_held_to_zero_to_one():
    (row, _other) = _table()
    assert oddities.probability(row) == pytest.approx(0.03)
    assert oddities.probability(row, dial=2.0) == pytest.approx(0.06)
    assert oddities.probability(row, dial=0.0) == 0.0
    assert oddities.probability(row, dial=1.0, scale={"W10": 0.5}) == pytest.approx(0.015)
    assert oddities.probability(row, dial=1.0, scale={"W14": 0.0}) == pytest.approx(0.03)
    assert oddities.probability(row, dial=1000.0) == 1.0


def test_the_dial_defaults_come_from_the_tuning_constants(monkeypatch):
    (row, _other) = _table()
    monkeypatch.setattr(tuning, "WEIRDNESS_DIAL", 3.0)
    monkeypatch.setattr(tuning, "WEIRDNESS_RATE_SCALE", {"W10": 0.5})
    assert oddities.probability(row) == pytest.approx(0.045)


def test_a_roll_follows_the_probability_over_many_bodies():
    (row, _other) = _table()
    hits = sum(oddities.roll(row, body_key=key, dial=10.0) for key in range(20000))
    assert hits / 20000 == pytest.approx(0.3, abs=0.02)
    assert sum(oddities.roll(row, body_key=key, dial=0.0) for key in range(2000)) == 0
    assert all(oddities.roll(row, body_key=key, dial=1000.0) for key in range(200))


# --- the draws are isolated -----------------------------------------------------------------------------------------------------------------

def test_a_roll_is_a_function_of_the_row_and_the_body_alone():
    (row, other) = _table()
    first = [oddities.roll(row, body_key=key, dial=10.0) for key in range(300)]
    again = [oddities.roll(row, body_key=key, dial=10.0) for key in range(300)]
    assert first == again
    assert oddities.stream("W10", 7).random() == oddities.stream("W10", 7).random()
    assert oddities.stream("W10", 7).random() != oddities.stream("W14", 7).random()
    assert oddities.stream("W10", 7).random() != oddities.stream("W10", 8).random()
    assert first != [oddities.roll(other, body_key=key, dial=10.0) for key in range(300)]


def test_a_roll_does_not_touch_the_run_stream():
    (row, _other) = _table()
    with draw.bound(1036):
        before = draw.getstate()
        for key in range(50):
            oddities.roll(row, body_key=key)
        assert draw.getstate() == before


def test_the_other_draws_of_a_run_are_the_same_with_the_oddity_rolls_or_without_them_and_at_any_dial():
    (row, other) = _table()

    def run(rolls, dial):
        out = []
        with draw.bound(2026):
            for key in range(40):
                if rolls:
                    for odd in (row, other):
                        oddities.roll(odd, body_key=key, dial=dial)
                out.append(draw.random())
        return out

    assert run(False, 1.0) == run(True, 1.0) == run(True, 0.0) == run(True, 5.0)


def test_adding_a_row_does_not_change_the_rolls_of_the_rows_already_there():
    (row, other) = _table()
    alone = [oddities.roll(row, body_key=key, dial=10.0) for key in range(500)]
    random.Random(1).shuffle(rows := [other, row])
    together = [[oddities.roll(odd, body_key=key, dial=10.0) for odd in rows] for key in range(500)]
    assert [pair[rows.index(row)] for pair in together] == alone


def test_the_stream_does_not_depend_on_the_hash_seed():
    # a string seed is hashed with SHA-512 by random.Random, not by hash(), so it is the same in every process
    import subprocess
    code = "from planetgen.physics import oddities; print(oddities.stream('W10', 'planet-7').random())"
    outs = {subprocess.run([sys.executable, "-c", code], env={"PYTHONHASHSEED": seed, "PYTHONPATH": str(ROOT / "src"), "PATH": "/usr/bin:/bin"},
                           capture_output=True, text=True, check=True).stdout for seed in ("1", "2")}
    assert len(outs) == 1
