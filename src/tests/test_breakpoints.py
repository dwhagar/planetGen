"""OPS.161: what a job declares when it enters the queue: the suggested split, where a cut may fall, the split flag and what it depends on (`planetgen.queue.breakpoints`)."""

import dataclasses

import pytest

from planetgen.queue import breakpoints as bp, commands, units


def test_a_declaration_holds_the_suggested_split_as_cuts_and_gives_the_sizes():
    decl = bp.even(1000, 300)
    assert decl.suggested == (300, 600, 900)
    assert decl.sizes() == [300, 300, 300, 100]
    assert decl.spans() == [(0, 300), (300, 600), (600, 900), (900, 1000)]
    assert decl.split == bp.SUGGESTED and not decl.dynamic, "a job that declares no flag uses its suggested split"


def test_the_split_flag_is_one_of_two_words():
    assert bp.even(10, 5, split=bp.BYTE_IDENTICAL).dynamic
    with pytest.raises(bp.DeclarationError):
        bp.even(10, 5, split="maybe")


def test_a_cut_may_not_fall_inside_a_span_that_stays_whole_or_off_the_boundaries():
    decl = bp.Declaration(100, (), whole=((20, 40),))
    assert 30 not in decl.legal_cuts and 20 in decl.legal_cuts and 40 in decl.legal_cuts and 19 in decl.legal_cuts
    decl = bp.Declaration(100, (), boundaries=(10, 50, 90))
    assert decl.legal_cuts == (10, 50, 90)
    decl = bp.Declaration(100, (), boundaries=(10, 50, 90), whole=((40, 60),))
    assert decl.legal_cuts == (10, 90)


def test_a_suggested_split_must_itself_be_legal():
    with pytest.raises(bp.DeclarationError, match="where a cut may not fall"):
        bp.Declaration(100, (30,), whole=((20, 40),))


@pytest.mark.parametrize("make", [
    lambda: bp.Declaration(-1),
    lambda: bp.Declaration(10, (0,)),
    lambda: bp.Declaration(10, (10,)),
    lambda: bp.Declaration(10, (5, 3)),
    lambda: bp.Declaration(10, (), whole=((5, 5),)),
    lambda: bp.Declaration(10, (), whole=((5, 11),)),
    lambda: bp.Declaration(10, (), floor_items=0),
    lambda: bp.Declaration(10, (), floor_items=5, ceiling_items=2),
    lambda: bp.even(10, 0),
])
def test_a_declaration_that_contradicts_itself_is_refused(make):
    with pytest.raises(bp.DeclarationError):
        make()


def test_the_split_a_stage_makes_today_is_declared_from_its_chunks():
    chunks = units.static_chunks([(i, 1.0) for i in range(100)], units=4)
    decl = bp.from_chunks(chunks)
    assert decl.items == 100 and decl.sizes() == [len(chunk) for chunk in chunks]


def test_a_job_with_no_item_or_one_unit_by_nature_has_nothing_to_cut():
    assert bp.Declaration(0).sizes() == [] and bp.Declaration(0).spans() == []
    assert bp.whole_job().sizes() == [1] and bp.whole_job().legal_cuts == ()


def test_the_unit_is_kept_between_a_quarter_and_four_times_the_suggested_one_unless_the_job_says_otherwise():
    decl = bp.even(1000, 100)
    assert decl.clamp() == (25, 400)
    assert bp.even(1000, 100, floor_items=10, ceiling_items=150).clamp() == (10, 150)
    assert bp.even(1000, 100, floor_items=500).clamp() == (500, 500), "a floor above the default ceiling lifts it"


def test_two_jobs_overlap_only_when_both_declared_and_neither_changes_what_the_other_uses():
    a = bp.Declaration(1, reads=("sectors",), writes=("links",))
    b = bp.Declaration(1, reads=("sectors",), writes=("paths",))
    c = bp.Declaration(1, reads=("links",))
    d = bp.Declaration(1, writes=("links",))
    assert bp.may_overlap(a, b), "both read the sectors and write different things"
    assert not bp.may_overlap(a, c), "c reads what a writes"
    assert not bp.may_overlap(a, d), "both write the links"
    assert not bp.may_overlap(a, bp.Declaration(1)), "a job that declares nothing is unknown, so it runs in plan order"


def test_a_whole_command_declares_itself_and_one_that_declares_nothing_keeps_its_typed_split():
    assert commands.declaration("plan-skeleton") == bp.whole_job()
    assert commands.declaration("scatter-mass") is None


def test_a_command_that_names_no_dependency_takes_the_ones_the_queue_already_holds(monkeypatch):
    assert commands.declaration("run-territories").after == ("run-population",)
    assert commands.declaration("finish-phenomena").after == ("scatter-phenomena",)


def _own(params):
    return bp.Declaration(params["items"], after=("scatter-mass",))


def _plain(params):
    return bp.even(params["items"], 10)


def test_a_command_that_declares_gets_its_declaration_for_the_params_of_the_stage(monkeypatch):
    spec = commands.REGISTRY["run-territories"]
    monkeypatch.setitem(commands.REGISTRY, "run-territories", dataclasses.replace(spec, breaks=("tests.test_breakpoints", "_own")))
    assert commands.declaration("run-territories", {"items": 7}) == bp.Declaration(7, after=("scatter-mass",)), "its own dependencies win"
    monkeypatch.setitem(commands.REGISTRY, "run-territories", dataclasses.replace(spec, breaks=("tests.test_breakpoints", "_plain")))
    assert commands.declaration("run-territories", {"items": 25}).after == ("run-population",)
    assert commands.declaration("run-territories", {"items": 25}).sizes() == [10, 10, 5]


def test_the_split_flag_of_a_command_is_its_spec_and_wins_over_what_its_declaration_says(monkeypatch):
    spec = commands.REGISTRY["run-territories"]
    assert spec.split == commands.SPLIT_SUGGESTED and commands.SPLIT_BYTE_IDENTICAL == bp.BYTE_IDENTICAL and commands.SPLIT_SUGGESTED == bp.SUGGESTED
    monkeypatch.setitem(commands.REGISTRY, "run-territories", dataclasses.replace(spec, breaks=("tests.test_breakpoints", "_plain"), split=commands.SPLIT_BYTE_IDENTICAL))
    assert commands.declaration("run-territories", {"items": 25}).dynamic
    monkeypatch.setitem(commands.REGISTRY, "run-territories", dataclasses.replace(spec, breaks=("tests.test_breakpoints", "_flagged")))
    assert not commands.declaration("run-territories", {"items": 25}).dynamic, "the spec's flag, not the declaration's own"


def _flagged(params):
    return bp.even(params["items"], 10, split=bp.BYTE_IDENTICAL)
