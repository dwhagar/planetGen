"""PERF.202: every generating verb is a stage of the shared run (`followed.run_followed`), not a one-step run of its own."""

import pathlib

import planetgen

ALLOWED = {
    "generation/followed.py": "the one place a verb's stages are started as a run",
    "generation/finish_run.py": "carries a stopped run on from its own record",
    "cli/datasets.py": "one-unit maintenance command (datasets check)",
    "cli/reset.py": "one-unit maintenance command (reset)",
    "cli/warm.py": "one-unit maintenance command (warm)",
    "cli/checkdb.py": "one-unit maintenance command (check-db)",
    "cli/generate.py": "one-unit maintenance command (math check)",
}


def test_no_verb_starts_a_run_of_its_own_outside_the_shared_run():
    root = pathlib.Path(planetgen.__file__).parent
    starters = {path.relative_to(root).as_posix() for path in root.rglob("*.py")
                if "wrapper.run(" in path.read_text(encoding="utf-8")}
    assert starters <= set(ALLOWED), f"a verb starts its own run: {sorted(starters - set(ALLOWED))}"
