# tests/test_sweep_life.py

"""
TEST.452: the class-blind life test (GEN.915, Boss 2:46 PM EDT 2026-10-29: life depends only on the conditions, never on the class
letter). Every terrestrial body of the shared sweep is given each of the 20 class letters with its conditions unchanged, and the life model
and the habitability index must answer the same on every letter; the life code is scanned for reads of the class letter and each one must be
listed with its reason. The harness is proved able to fail by breaking the life model so that one letter changes a score.

The first half tests the rule on synthetic bodies.
"""

import types

import pytest

from planetgen import tuning
from planetgen.physics import atmosphere, constants, life_classes
from tests.audit import prevalence as census
from tests.audit import sweep
from tests.audit import sweep_life as rules

REPORT_NAME = "class-blind-life-latest.md"


def _body(**changes):
    """Earth's stored conditions as a body object; its class letter is M."""
    fields = dict(
        body_type="t", mass=constants.EARTH_MASS_TO_KG, gravity=1.0, surface_temperature=288.0, atmospheric_pressure=101325.0,
        hydrosphere="surface ocean", ocean_fraction=0.7, water_mass_fraction=2.3e-4, flux_variation=0.0, ocean_ph=None,
        planet_class="M", **{atmosphere.GAS_FIELD: {"O2": 21200.0, "N2": 78800.0, "Ar": 900.0, "CO2": 40.0}})
    fields.update(changes)
    return types.SimpleNamespace(**fields)


def _entry(body=None):
    return sweep.Entry(body or _body(), None, None, "x")


@pytest.fixture(scope="module")
def notes():
    return rules.load_notes()


# --- The rule on bodies made for it ----------------------------------------------------------------------------------------


def test_there_are_twenty_class_letters_to_try():
    assert len(rules.letters()) == 20 == len(tuning.PLANET_CLASSES)


def test_an_earth_like_body_scores_the_same_on_every_letter():
    assert rules._check_scores_class_blind(_entry()) == []


def test_a_gas_giant_is_not_judged():
    assert rules._check_scores_class_blind(_entry(_body(body_type="g"))) == []


# --- The source scan -------------------------------------------------------------------------------------------------------


def test_every_read_of_the_class_letter_in_the_life_code_is_listed_with_its_reason(notes):
    found = rules.class_references()
    unlisted = rules.unlisted_references(found, notes)
    assert not unlisted, "the life code reads the class letter where the notes list no reason: " + ", ".join(
        f"{place} (lines {lines})" for place, lines in unlisted.items())
    assert all(reason for reason in notes["allowed_references"].values())


def test_the_scan_finds_a_read_that_is_planted(monkeypatch, tmp_path, notes):
    planted = tmp_path / "planted_life.py"
    planted.write_text("def draw_life(body):\n    return 0.0 if body.planet_class == 'C' else 1.0\n\ndef harmless(body):\n    return body.mass\n")
    import importlib
    import sys

    monkeypatch.syspath_prepend(str(tmp_path))
    monkeypatch.setattr(rules, "LIFE_MODULES", ("planted_life",))
    found = rules.class_references()
    sys.modules.pop("planted_life", None)
    importlib.invalidate_caches()
    assert list(found) == ["planted_life:draw_life"] and rules.unlisted_references(found, notes) == found


def test_a_listed_place_that_no_longer_reads_the_letter_is_reported_stale(notes):
    assert rules.stale_references({}, notes) == sorted(notes["allowed_references"])
    assert rules.stale_references({place: [1] for place in notes["allowed_references"]}, notes) == []


# --- The sweep -------------------------------------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def everything():
    return sweep.everything()


@pytest.mark.slow
def test_the_life_model_and_the_index_answer_the_same_on_every_letter_for_every_body_of_the_sweep(everything, notes):
    found = sweep.run_rules(rules.RULES, everything)["scores-class-blind"]
    capable = rules.capable_by_class(everything)
    total = sum(count for count, _ in capable.values())
    lines = [f"scores-class-blind: {len(found)} of {len(everything)} bodies differ with the class letter",
             f"terrestrial bodies {total}; the life model gives a class to {sum(ok for _, ok in capable.values())}, by the letter each was drawn with:"]
    lines += [f"  {letter}: {ok} of {count}" for letter, (count, ok) in sorted(capable.items())]
    stale = rules.stale_references(rules.class_references(), notes)
    if stale:
        lines.append("STALE NOTES: these listed places no longer read the class letter; remove them from tests/data/class_blind_life_notes.json: "
                     + ", ".join(stale))
    census.write_report("\n".join(lines) + "\n", name=REPORT_NAME)
    assert not found, f"{len(found)} bodies whose life depends on the class letter, first {found[:5]}"


@pytest.mark.slow
def test_the_sweep_holds_life_capable_conditions_so_the_letters_are_tried_on_them(everything):
    capable = rules.capable_by_class(everything)
    assert sum(ok for _, ok in capable.values()) >= 50
    assert len([letter for letter, (_, ok) in capable.items() if ok]) >= 10


@pytest.fixture(scope="module")
def capable_entry(everything):
    """The sweep's terrestrial body to which the life model gives the highest score: a body the proofs below can see move."""
    best = max((entry for entry in everything if entry.body.body_type == "t"),
               key=lambda entry: max([result.score for result in life_classes.score_classes(entry.body).values()] or [0]))
    assert rules._scores(best.body)[1] is not None
    return best


@pytest.mark.slow
def test_a_life_model_that_reads_the_letter_is_reported_with_the_letter_and_the_class(monkeypatch, capable_entry):
    real = life_classes.score_classes
    letter = next(letter for letter in rules.letters() if letter != capable_entry.body.planet_class)
    monkeypatch.setattr(life_classes, "score_classes", lambda body: {} if body.planet_class == letter else real(body))
    problems = rules._check_scores_class_blind(capable_entry)
    assert problems and all(f"as class {letter}" in problem for problem in problems), problems
    assert any("the best life class is None" in problem for problem in problems)
    assert any(rules._scores(capable_entry.body)[1] in problem for problem in problems)


@pytest.mark.slow
def test_a_single_score_that_moves_with_the_letter_names_that_life_class(monkeypatch, capable_entry):
    real = life_classes.score_classes
    class_id = rules._scores(capable_entry.body)[1]
    letter = next(letter for letter in rules.letters() if letter != capable_entry.body.planet_class)

    def nudged(body):
        scores = real(body)
        if body.planet_class == letter:
            scores[class_id] = types.SimpleNamespace(id=class_id, likelihood=scores[class_id].likelihood, score=scores[class_id].score - 1)
        return scores

    monkeypatch.setattr(life_classes, "score_classes", nudged)
    problems = rules._check_scores_class_blind(capable_entry)
    assert any(f"life class {class_id} differ" in problem and f"as class {letter}" in problem for problem in problems), problems


@pytest.mark.slow
def test_an_index_that_reads_the_letter_is_reported(monkeypatch, capable_entry):
    from planetgen.physics import habitability_world

    real = habitability_world.compute
    letter = next(letter for letter in rules.letters() if letter != capable_entry.body.planet_class)

    def reads_the_letter(values, compact_host=None):
        result = dict(real(values, compact_host))
        if values.get("planet_class") == letter:
            result["phi_bio"] = -1.0
        return result

    monkeypatch.setattr(habitability_world, "compute", reads_the_letter)
    problems = rules._check_scores_class_blind(capable_entry)
    assert any(f"as class {letter} the habitability index differs" in problem and "phi_bio" in problem for problem in problems), problems
