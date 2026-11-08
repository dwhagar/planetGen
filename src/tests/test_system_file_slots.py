# tests/test_system_file_slots.py

"""
GEN.688: a system file's planet slot names conditions (a letter is only sugar), and a slot that cannot be read fails in plain words before anything is generated:
`run_system.slots_problem`, the command line's usage error and the Generate page's form error.
"""

import json
import sys

import pytest

from planetgen.cli import generate as generate_cli
from planetgen.generation import run_system
from planetgen.generation.config import SystemConfig
from planetgen.web import generate_page, system_page

WARM = {"kind": "t", "air": True, "temperature_k": [273, 320]}


@pytest.mark.parametrize("slots", [
    None, [], [None], [None, None],
    [{"type": "planet"}],
    [{"type": "planet", "planet_class": "M", "moons": 0}],
    [{"type": "planet", "conditions": WARM, "moons": 2}],
    [{"conditions": {"kind": "g"}}],
    [{"type": "asteroid_belt"}, None, {"type": "planet", "planet_class": "J"}],
])
def test_a_slot_list_that_can_be_read_has_no_problem(slots):
    assert run_system.slots_problem(slots) is None


@pytest.mark.parametrize("slots, message", [
    ({"type": "planet"}, "must be a list"),
    ("M", "must be a list"),
    (["M"], "slot 1 must be null or an object"),
    ([None, {"type": "moon"}], "slot 2 has the type 'moon'"),
    ([{"type": "planet", "class": "M"}], "slot 1 has 'class'"),
    ([{"type": "planet", "planet_class": "M", "conditions": WARM}], "not both"),
    ([{"type": "planet", "planet_class": "Q"}], "'Q' is not a planet class"),
    ([{"type": "planet", "conditions": {"flavor": [1, 2]}}], "is not a quantity"),
    ([{"type": "planet", "conditions": {"temperature_k": 300}}], "must be a range"),
    ([{"type": "planet", "conditions": {"temperature_k": [320, 273]}}], "low end under its high end"),
    ([{"type": "planet", "conditions": {"kind": "x"}}], "kind must be"),
    ([{"type": "planet", "conditions": [1]}], "must be an object"),
    ([{"type": "planet", "moons": -1}], "moons must be a whole number from 0"),
    ([{"type": "planet", "moons": True}], "moons must be a whole number from 0"),
    ([{"type": "asteroid_belt", "planet_class": "M"}], "asteroid belt, which takes no planet_class"),
    ([{"type": "asteroid_belt", "conditions": WARM}], "takes no conditions"),
])
def test_a_slot_list_that_cannot_be_read_says_what_is_wrong(slots, message):
    assert message in run_system.slots_problem(slots)


def test_the_message_names_the_slot_by_its_place():
    assert run_system.slots_problem([None, None, {"type": "planet", "moons": "two"}]).startswith("slot 3:")


def test_the_slot_keys_are_the_ones_the_config_documents():
    import inspect
    assert set(run_system.SLOT_KEYS) == {"type", "planet_class", "conditions", "moons"}
    assert all(f'"{key}"' in inspect.getsource(SystemConfig) for key in run_system.SLOT_KEYS)


def test_the_pages_example_slots_can_be_read():
    example = json.loads(system_page.SYSTEM_FILE_EXAMPLE)
    assert run_system.slots_problem(example["slots"]) is None
    assert any("conditions" in slot for slot in example["slots"] if slot)


def _form(slots):
    return {"system_file": json.dumps({"slots": slots})}


def test_the_page_accepts_a_slot_that_names_conditions():
    spec = system_page.system_request(_form([{"type": "planet", "conditions": WARM}]))
    assert spec["system_file"]["slots"][0]["conditions"] == WARM


def test_the_page_says_what_is_wrong_with_a_slot():
    with pytest.raises(generate_page.FormError, match="The system file's slot 1: .* not both"):
        system_page.system_request(_form([{"type": "planet", "planet_class": "M", "conditions": WARM}]))
    with pytest.raises(generate_page.FormError, match="is not a quantity"):
        system_page.system_request(_form([{"type": "planet", "conditions": {"flavor": [1, 2]}}]))


def _arguments(monkeypatch, tmp_path, slots):
    path = tmp_path / "system.json"
    path.write_text(json.dumps({"slots": slots}))
    monkeypatch.setattr(sys, "argv", ["planetgen", "system", "--system-file", str(path)])
    return generate_cli.process_args()


def test_the_command_line_accepts_conditions(monkeypatch, tmp_path):
    arguments = _arguments(monkeypatch, tmp_path, [{"type": "planet", "conditions": WARM}])
    assert arguments.system_file


def test_the_command_line_refuses_a_bad_slot_in_plain_words(monkeypatch, tmp_path, capsys):
    with pytest.raises(SystemExit) as exit_:
        _arguments(monkeypatch, tmp_path, [{"type": "planet", "planet_class": "M", "conditions": WARM}])
    assert exit_.value.code == 2
    assert "slot 1" in capsys.readouterr().err
