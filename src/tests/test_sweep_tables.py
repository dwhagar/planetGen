# tests/test_sweep_tables.py

"""
TEST.214: the possibility-table sweep. Walks every table of `possibility_table.REGISTRY` (a table added later is covered without a change
here) with seeded subjects of every class the table names, and checks what the composition model promises: every row has a source and a
state; every named class reaches a row; the gated weights of every subject add up to 100 percent within 1e-12; no subject fails but with an
empty gate (apart from the tables whose rows are independent flags, listed with a note); the fractions of every row add up to 1; `Table.draw`
picks the rows as often as their probabilities say; every row is reached for every class it names; the open rows are the ones listed with
a note; and a table that is not registered fails.

The first half tests each check on synthetic tables, where the answer is known.
"""

import json

import pytest

from planetgen.generation import possibility_table as pt
from planetgen.generation.possibility_table import Condition, Registry, Row, Table
from planetgen.util import draw
from tests.audit import sweep_tables as sweep


def _row(name, weight=1.0, classes=("rocky",), **kw):
    kw.setdefault("source", "Test et al. 2020")
    return Row(name, classes, weight, **kw)


def _open(name, classes=("rocky",)):
    return Row(name, classes, pt.OPEN, "open", state=pt.OPEN_STATE)


def _registry(*tables):
    registry = Registry()
    for table in tables:
        registry.register(table)
    return registry


def _outcomes(table, registry=None, **kw):
    registry = registry or _registry(table)
    return {c: sweep.run_samples(table, c, registry, count=20, **kw) for c in sweep.table_classes(table)}


# --- Each check on a table made for it ---------------------------------------------------------------------------------------


def test_a_sound_table_raises_no_flag():
    table = Table("sound", 1, [_row("a", 0.3), _row("b", 0.7)])
    outcomes = _outcomes(table)
    assert sweep.check_row_states(table) == []
    assert sweep.check_samples(table, outcomes) == ([], [], [])
    assert sweep.check_draw_follows_weights(table, outcomes) == []
    assert sweep.unreached(table, outcomes, _registry(table)) == []


def test_a_row_without_a_source_or_a_limit_or_a_replacement_is_reported():
    bad = Row("x", ("rocky",), 1.0, "Test et al. 2020")
    object.__setattr__(bad, "source", " ")
    partly = Row("y", ("rocky",), 1.0, "Test", state=pt.PARTLY_VERIFIED, limit="a limit")
    object.__setattr__(partly, "limit", "")
    stand_in = Row("z", ("rocky",), 1.0, "Test", state=pt.STAND_IN, replacement="GEN.1")
    object.__setattr__(stand_in, "replacement", "")
    found = sweep.check_row_states(Table("bad", 1, [bad, partly, stand_in]))
    assert len(found) == 3 and "no source" in found[0] and "no limit" in found[1] and "no replacement" in found[2]


def test_an_open_state_and_an_open_weight_must_go_together():
    row = _row("odd")
    object.__setattr__(row, "state", pt.OPEN_STATE)
    assert any("open" in finding for finding in sweep.check_row_states(Table("odd", 1, [row])))


def test_a_class_that_no_subject_reaches_is_reported_as_an_empty_gate():
    table = Table("gated", 1, [_row("a", classes=("rocky",), conditions=(Condition("x", minimum=5.0, maximum=1.0),))])
    gate, _sums, _failed = sweep.check_samples(table, _outcomes(table))
    assert gate and "rocky" in gate[0]


def test_weights_that_go_over_one_beside_an_open_row_are_a_failed_subject():
    table = Table("over", 1, [_row("a", 0.8), _row("b", 0.7), _open("c")])
    _gate, _sums, failed = sweep.check_samples(table, _outcomes(table))
    assert failed and "over 1" in failed[0]


def test_an_all_zero_answer_fails_unless_the_table_and_class_are_allowed_with_a_note():
    table = Table("zero", 1, [_row("a", 0.0)])
    outcomes = _outcomes(table)
    assert sweep.check_samples(table, outcomes)[2]
    allowed = {"all_zero_allowed": {"zero": {"classes": ["rocky"], "note": "a flag"}}}
    assert sweep.check_samples(table, outcomes, allowed)[2] == []
    assert sweep.check_samples(table, outcomes, {"all_zero_allowed": {"zero": {"classes": ["rocky"], "note": " "}}})[2]
    assert sweep.check_samples(table, outcomes, {"all_zero_allowed": {"zero": {"classes": ["icy"], "note": "x"}}})[2]
    assert sweep.check_samples(table, outcomes, {"all_zero_allowed": {"zero": {"classes": ["*"], "note": "x"}}})[2] == []


def test_a_weight_that_raises_is_a_failed_subject():
    def explode(subject):
        raise ArithmeticError("boom")

    table = Table("raises", 1, [_row("a", explode)])
    _gate, _sums, failed = sweep.check_samples(table, _outcomes(table))
    assert failed and "boom" in failed[0]


def test_probabilities_that_do_not_add_up_are_reported():
    outcome = sweep.Outcome("rocky", {}, "ok", probabilities=(("a", 0.5), ("b", 0.4)))
    table = Table("t", 1, [_row("a")])
    assert sweep.check_samples(table, {"rocky": [outcome]})[1]
    outcome = sweep.Outcome("rocky", {}, "ok", probabilities=(("a", 1.2), ("b", -0.2)))
    assert sweep.check_samples(table, {"rocky": [outcome]})[1]


def test_a_dead_row_is_reported_unreached_and_one_with_a_narrow_gate_is_reached():
    narrow = _row("narrow", classes=("rocky",), conditions=(Condition("x", minimum=300.0, maximum=301.0),))
    dead = _row("dead", 0.0)
    table = Table("reach", 1, [_row("live"), narrow, dead])
    registry = _registry(table)
    found = sweep.unreached(table, _outcomes(table, registry), registry)
    assert found == ["reach/dead/rocky"]


def test_the_fractions_of_a_row_must_add_up_and_name_known_components():
    materials = Table("materials", 1, [_row("iron", classes=("*",)), _row("rock", classes=("*",))])
    good = Table("good", 1, [_row("a", fractions={"iron": 1.0, "rock": 3.0}, concentration=5.0)], components="materials")
    assert sweep.check_fractions(good, _registry(materials, good)) == []
    unknown = Table("unknown", 1, [_row("a", fractions={"iron": 1.0, "unobtainium": 1.0})], components="materials")
    assert any("unobtainium" in finding for finding in sweep.check_fractions(unknown, _registry(materials, unknown)))


def test_a_row_whose_drawn_fractions_do_not_add_up_is_reported(monkeypatch):
    materials = Table("materials", 1, [_row("iron", classes=("*",))])
    table = Table("bent", 1, [_row("a", fractions={"iron": 1.0})], components="materials")
    monkeypatch.setattr(Row, "draw_fractions", lambda self: {"iron": 0.9})
    assert sweep.check_fractions(table, _registry(materials, table))


def test_a_draw_that_ignores_the_weights_is_reported(monkeypatch):
    table = Table("skew", 1, [_row("a", 0.9), _row("b", 0.1)])
    outcomes = _outcomes(table)
    monkeypatch.setattr(draw, "choices", lambda population, weights=None, k=1: [population[-1]] * k)
    assert sweep.check_draw_follows_weights(table, outcomes)


def test_an_open_row_without_a_note_is_new_and_a_note_for_a_row_that_is_gone_is_stale():
    table = Table("opened", 1, [_row("a", 0.5), _open("b")])
    registry = _registry(table)
    new, stale = sweep.check_open_rows(registry, {"open_rows": {}})
    assert new == ["opened/b"] and stale == []
    new, stale = sweep.check_open_rows(registry, {"open_rows": {"opened": {"b": "waits on a study"}, "old": {"z": "gone"}}})
    assert new == [] and stale == ["old/z"]
    new, _stale = sweep.check_open_rows(registry, {"open_rows": {"opened": {"b": " "}}})
    assert new == ["opened/b (no note)"]


def test_an_unreached_row_needs_a_note_and_a_note_for_a_reached_row_is_stale():
    assert sweep.check_unreached(["t/r/c"], {}) == (["t/r/c"], [])
    assert sweep.check_unreached(["t/r/c"], {"unreached_allowed": {"t/r/c": "why"}}) == ([], [])
    assert sweep.check_unreached([], {"unreached_allowed": {"t/r/c": "why"}}) == ([], ["t/r/c"])


def test_a_table_that_is_not_the_registered_one_of_its_name_fails():
    registered = Table("one", 1, [_row("a")])
    registry = _registry(registered)
    stray = Table("two", 1, [_row("a")])
    twin = Table("one", 1, [_row("a")])
    assert sweep.unregistered([registered], registry) == []
    assert sweep.unregistered([registered, stray, twin], registry) == ["one", "two"]


def test_the_samples_repeat():
    table = Table("repeat", 1, [_row("a", lambda subject: subject["x"] + 1.0, conditions=(Condition("x", minimum=0.0),))])
    first = [(o.kind, o.quantities) for o in sweep.run_samples(table, "rocky", _registry(table), count=10)]
    assert first == [(o.kind, o.quantities) for o in sweep.run_samples(table, "rocky", _registry(table), count=10)]


# --- The tables of the generator ---------------------------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def registry():
    return sweep.load_tables()


@pytest.fixture(scope="module")
def notes():
    return sweep.load_notes()


@pytest.fixture(scope="module")
def outcomes(registry):
    return {table.name: {c: sweep.run_samples(table, c, registry) for c in sweep.table_classes(table)} for table in registry.tables()}


def test_the_walk_finds_the_tables_of_the_generator(registry):
    names = {table.name for table in registry.tables()}
    assert len(registry) >= 33 and {f"atmosphere-{letter}" for letter in "abcdefghijklmnop"} <= names
    assert {"materials", "moon-origin", "star-count", "bulk-rock-kind", "comet-carbon-chain"} <= names


def test_every_row_of_every_table_has_a_source_and_a_state(registry):
    found = [finding for table in registry.tables() for finding in sweep.check_row_states(table)]
    assert not found, "\n".join(found)


def test_every_class_of_every_table_reaches_a_row_and_every_gated_sum_is_one(registry, outcomes, notes):
    gate, sums, failed = [], [], []
    for table in registry.tables():
        found = sweep.check_samples(table, outcomes[table.name], notes)
        gate += found[0]
        sums += found[1]
        failed += found[2]
    assert not (gate or sums or failed), "\n".join(["empty gates:", *gate[:5], "sums:", *sums[:5], "failed subjects:", *failed[:5]])


def test_every_all_zero_allowance_carries_a_note_and_names_a_table_of_the_registry(registry, notes):
    for name, entry in notes["all_zero_allowed"].items():
        assert name in registry and entry["classes"] and entry["note"].strip(), name


def test_the_fractions_of_every_row_add_up_to_one(registry):
    found = [finding for table in registry.tables() for finding in sweep.check_fractions(table, registry)]
    assert not found, "\n".join(found)


def test_the_registry_names_only_components_of_the_shared_table(registry):
    registry.validate()


def test_the_draw_picks_each_row_as_often_as_its_probability(registry, outcomes):
    found = [finding for table in registry.tables() for finding in sweep.check_draw_follows_weights(table, outcomes[table.name])]
    assert not found, "\n".join(found)


def test_every_row_is_reached_for_every_class_it_names_or_has_a_note(registry, outcomes, notes):
    reached = [item for table in registry.tables() for item in sweep.unreached(table, outcomes[table.name], registry)]
    new, stale = sweep.check_unreached(reached, notes)
    assert not new, "rows no subject reaches (add a note to tests/data/possibility_sweep_notes.json if it is meant):\n" + "\n".join(new)
    assert isinstance(stale, list)


def test_the_open_rows_are_the_ones_listed_with_a_note(registry, notes):
    new, stale = sweep.check_open_rows(registry, notes)
    assert not new, "open rows with no note in tests/data/possibility_sweep_notes.json:\n" + "\n".join(new)
    assert isinstance(stale, list)


def test_the_notes_file_is_json_with_the_three_sections():
    text = sweep.NOTES_FILE.read_text(encoding="utf-8")
    assert set(json.loads(text)) == {"open_rows", "all_zero_allowed", "unreached_allowed"}


@pytest.mark.slow
def test_every_table_a_fresh_process_holds_is_registered():
    assert sweep.unregistered_in_fresh_process() == []
