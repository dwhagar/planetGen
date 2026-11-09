# tests/test_train_affects.py

"""A train's epoch alters the stages its notes say (`generation-affects:`), not every stage by default."""

import importlib.util
import json
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("train", ROOT / "scripts" / "train.py")
train = importlib.util.module_from_spec(spec)
spec.loader.exec_module(train)
bump = train.bump_version

NOTE = "**Headline:** A thing.\n\n### Changed\n- Something.\n"


def test_the_stage_list_is_the_epoch_files():
    stages = json.loads((ROOT / "src" / "tests" / "golden" / "epochs.json").read_text(encoding="utf-8"))["stages"]
    assert list(bump.GENERATION_STAGES) == stages


@pytest.mark.parametrize("line, expected", [
    ("generation-affects: fill", ("fill",)),
    ("generation-affects: population, scatter", ("scatter", "population")),
    ("generation-affects:   Plan , FILL  ", ("plan", "fill")),
    ("generation-affects: none", ()),
])
def test_a_note_line_gives_its_stages_in_epoch_order(line, expected):
    assert bump.parse_generation_affects(NOTE + line + "\n") == expected


def test_a_note_without_the_line_says_nothing():
    assert bump.parse_generation_affects(NOTE) is None


@pytest.mark.parametrize("line", ["generation-affects:", "generation-affects: fill, bogus", "generation-affects: none, fill"])
def test_a_wrong_line_is_a_problem_in_the_note_check(line):
    with pytest.raises(ValueError):
        bump.parse_generation_affects(NOTE + line + "\n")
    assert any("generation-affects" in problem for problem in bump.validate_fragment_body("x.patch.md", NOTE + line + "\n"))


def test_a_good_line_passes_the_note_check():
    assert bump.validate_fragment_body("x.patch.md", NOTE + "generation-affects: fill\n") == []


def test_the_union_of_the_declared_stages_is_the_epochs():
    notes = {"a.patch.md": NOTE + "generation-affects: fill\n", "b.patch.md": NOTE + "generation-affects: scatter, fill\n",
             "c.patch.md": NOTE + "generation-affects: none\n"}
    assert train.stages_of_notes(notes) == (("scatter", "fill"), [], [])


def test_a_note_with_no_line_counts_as_every_stage_and_is_named():
    notes = {"a.patch.md": NOTE + "generation-affects: fill\n", "b.patch.md": NOTE}
    stages, undeclared, problems = train.stages_of_notes(notes)
    assert stages == ("plan", "scatter", "fill", "population") and undeclared == ["b.patch.md"] and problems == []


def test_a_wrong_line_is_reported_and_does_not_widen_the_epoch():
    stages, undeclared, problems = train.stages_of_notes({"a.patch.md": NOTE + "generation-affects: nope\n"})
    assert stages == () and undeclared == [] and len(problems) == 1 and "a.patch.md" in problems[0]


def test_the_affects_line_is_not_release_text():
    """The line is for the train runner; the changelog entry made from a note leaves it out."""
    body = "**Headline:** A thing.\n\n### Changed\n- It changed.\n\ngeneration-affects: fill\n"
    merged = bump.merge_bodies(body, "")
    assert "generation-affects" not in merged and "- It changed." in merged
