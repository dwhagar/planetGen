"""
DOC.72: the README lists every external file the program needs, and nothing the catalogue does not hold. The list is written from the catalogue
(`scripts/readme_external_files.py`); this test fails when a dataset is added (an AD item, GEN.594) or a URL changes and the README is not rewritten,
when the README names a set or URL the catalogue lacks, and when the statements about a missing file no longer match the code.
"""

import importlib.util
import pathlib
import re

from planetgen.datasets import catalogue, fetch
from planetgen.generation import run_galaxy

ROOT = pathlib.Path(__file__).resolve().parents[2]
README = (ROOT / "README.md").read_text(encoding="utf-8")


def _script():
    spec = importlib.util.spec_from_file_location("readme_external_files", ROOT / "scripts" / "readme_external_files.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _listed():
    module = _script()
    return README.split(module.START, 1)[1].split(module.END, 1)[0]


def test_the_readme_section_is_the_one_written_from_the_catalogue():
    module = _script()
    assert module.START in README and module.END in README
    assert README.split(module.START, 1)[1].split(module.END, 1)[0].strip() == module.section().split(module.START, 1)[1].split(module.END, 1)[0].strip(), \
        "run: python3 scripts/readme_external_files.py --write"


def test_every_dataset_and_every_url_of_the_catalogue_is_in_the_readme_and_nothing_else():
    listed = _listed()
    keys = re.findall(r"^### `([^`]+)`", listed, re.M)
    assert keys == [dataset.key for dataset in catalogue.ALL]
    urls = re.findall(r"<(https://[^>]+)>", listed)
    assert sorted(urls) == sorted(piece.url for dataset in catalogue.ALL for piece in dataset.pieces)
    for dataset in catalogue.ALL:
        assert f"python3 -m planetgen.cli.datasets fetch {dataset.key}" in listed
        assert dataset.replaces in listed


def test_each_set_says_what_it_falls_back_to_or_that_it_does_not():
    listed = _listed()
    for dataset in catalogue.ALL:
        block = listed.split(f"### `{dataset.key}`", 1)[1].split("\n### `", 1)[0]
        if dataset.fallback:
            assert dataset.fallback[0].upper() + dataset.fallback[1:] in block
        else:
            assert "nothing falls back" in block


def test_the_readme_statements_about_a_missing_file_match_the_code(monkeypatch, tmp_path):
    assert "--disregard-missing-data" in README
    assert fetch.CACHE_DIR_ENV in README and fetch.DEFAULT_CACHE_DIR in README
    for command in ("fetch", "status"):
        assert f"planetgen.cli.datasets {command}" in README
    # The refusal names the set, every file with its URL and the setup command: the README says it does.
    monkeypatch.setenv(fetch.CACHE_DIR_ENV, str(tmp_path))
    dataset = catalogue.ALL[0]
    text = fetch.locate(dataset).describe()
    assert dataset.key in text and "python3 -m planetgen.cli.datasets fetch" in text
    assert all(piece.url in text for piece in dataset.pieces)
    # The flag and the Generate page's box exist under the names the README gives.
    web = "".join(path.read_text(encoding="utf-8") for path in (ROOT / "src" / "planetgen" / "web").rglob("*") if path.suffix in (".html", ".py", ".js"))
    assert "--disregard-missing-data" in web and "Disregard missing data files" in web and "disregard_missing_data" in web


def test_the_deployment_guide_points_at_the_list():
    guide = (ROOT / "docs" / "deployment" / "README.md").read_text(encoding="utf-8")
    assert "README.md#external-files" in guide
