"""API.33: a count the commands and the Generate page show says what it counts."""

from rich.progress import Progress

from planetgen.generation import run_common, steps
from planetgen.queue import progress_file
from planetgen.util.format import counted
from planetgen.web import generate_page


def test_every_kind_with_a_unit_is_a_known_step_kind():
    assert set(steps.STEP_UNITS) <= set(steps.STEP_KINDS)


def test_a_top_level_bar_counts_something_and_a_detail_bar_names_its_own_items():
    assert steps.unit_of("sector") == "sector" and steps.unit_of("population") == "pass"
    assert steps.unit_of("a-kind-nobody-named") == steps.DEFAULT_UNIT
    assert steps.unit_of("stages:Scanning planets") is None and steps.unit_of(None) is None


def test_the_noun_goes_with_the_count():
    assert [counted("sector", n) for n in (0, 1, 2)] == ["sectors", "sector", "sectors"]


def test_the_command_line_bar_counts_in_its_unit():
    progress = Progress(run_common._CountColumn())
    task = progress.add_task("Neighbours", total=40, completed=12, unit="sector")
    assert run_common._CountColumn().render(progress.tasks[task]).plain == "12/40 sectors"
    one = progress.add_task("Population pass", total=1, completed=0, unit="pass")
    assert run_common._CountColumn().render(progress.tasks[one]).plain == "0/1 pass"
    bare = progress.add_task("  Scanning planets", total=500, completed=5)
    assert run_common._CountColumn().render(progress.tasks[bare]).plain.strip() == "5/500"


def test_a_step_hands_its_unit_to_the_bar():
    progress = run_common._generation_progress(disable=True)
    with steps.Step("Sector paths", "paths", 40, progress=progress, force=True) as bar:
        bar.update(completed=3)
        assert progress.tasks[0].fields["unit"] == "sector"


def test_the_progress_file_carries_the_unit_and_the_page_says_it(tmp_path, monkeypatch):
    path = tmp_path / "progress.json"
    monkeypatch.setenv(progress_file.ENV_VAR, str(path))
    progress_file.report(3, 40, "Neighbours", force=True, unit="sector")
    import json
    assert json.loads(path.read_text())["unit"] == "sector"
    text = generate_page._progress_text
    assert text({"description": "Neighbours", "completed": 3, "total": 40, "unit": "sector"}) == "Neighbours: 3 of 40 sectors"
    assert text({"description": "Population pass", "completed": 0, "total": 1, "unit": "pass"}) == "Population pass: 0 of 1 pass"
    assert text({"description": "Sectors", "completed": 12, "total": 40, "unit": "sector"}) == "Sectors: 12 of 40"
    assert text({"description": "Scanning", "completed": 5, "total": 9}) == "Scanning: 5 of 9"
