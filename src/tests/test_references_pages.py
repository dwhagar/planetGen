# tests/test_references_pages.py

"""
DOC.88 (References R2): the References section, Help > References (`/references`), the one place the site speaks out of universe.

Three parts are held here.

* The build (`scripts/build_reference_pages.py`): `docs/references/overview.json` is the same bytes every time, equals what the registries
  give now (`--check`), and its counts are the registries' own, counted again here a different way.
* The page (`web/references_pages.py`, `web/lib/references_site.py`): it renders without a database, opens with the Out of universe banner,
  tints the header, sits last in the Help menu after a divider, and its sentences keep the plain-words guide (the out-of-universe words
  are allowed on it).
* The setup step (`update.sh`, `deploy-common.sh`): the build runs at setup and a failure only warns.
"""

import json
import os
import pathlib
import re
import subprocess
import sys

import pytest

from planetgen.util import plainwords
from planetgen.web.lib import reference_urls, references_site

from tests.test_class_pages import _html, app, client  # noqa: F401

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import build_reference_pages as builder  # noqa: E402

REGISTRY = ROOT / "docs" / "design" / "research" / "stand-in-values.md"
BANNER_ID = 'id="out-of-universe-banner"'


# --- The build -------------------------------------------------------------------------------------------------------------------

def test_the_committed_pages_are_what_the_registries_give_now():
    assert builder.stale(builder.build()) == []
    assert builder.main(["--check"]) == 0


def test_a_build_is_the_same_bytes_every_time_and_ends_in_one_newline():
    first, second = builder.build(), builder.build()
    assert first == second
    for text in first.values():
        assert text.endswith("}\n") and not text.endswith("\n\n")
        assert json.loads(text)


def test_a_missing_or_changed_or_extra_file_is_stale(tmp_path):
    built = builder.build()
    assert builder.stale(built, str(tmp_path)) == list(built)
    for name, text in built.items():
        (tmp_path / name).write_text(text, encoding="utf-8")
    assert builder.stale(built, str(tmp_path)) == []
    name = next(iter(built))
    (tmp_path / name).write_text(built[name].replace("0", "1", 1), encoding="utf-8")
    (tmp_path / "left-over.json").write_text("{}\n", encoding="utf-8")
    assert builder.stale(built, str(tmp_path)) == [name, "left-over.json (nothing builds it)"]


def test_the_check_fails_and_names_the_file_when_the_committed_one_is_missing(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(builder, "OUT", str(tmp_path))
    assert builder.main(["--check"]) == 1
    assert "docs/references is out of date: overview.json" in capsys.readouterr().err
    assert builder.main([]) == 0 and (tmp_path / "overview.json").is_file()
    assert builder.main(["--check"]) == 0


def test_the_counts_are_the_registries_own_counted_a_second_way():
    """The list the registry prints (`stand-in-values.md`) states its row count and gives a tier by tweak table; both are read here as text."""
    overview = builder.build_overview()
    text = REGISTRY.read_text(encoding="utf-8")
    assert overview["stand_ins"]["rows"] == int(re.search(r"There are (\d+) rows below", text).group(1)) > 900
    cross = {}
    for line in text.splitlines():
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if re.fullmatch(r"T\d", cells[0]) and len(cells) >= 6 and all(cell.isdigit() for cell in cells[1:6]):
            cross[cells[0]] = [int(cell) for cell in cells[1:6]]
    assert sorted(cross) == [f"T{n}" for n in range(1, 10)]
    for code, number in overview["stand_ins"]["by_tier"].items():
        assert number == sum(cross[code]), code
    for index, (code, number) in enumerate(overview["stand_ins"]["by_tweak"].items()):
        assert number == sum(row[index] for row in cross.values()), code
    assert sum(overview["stand_ins"]["by_tweak"].values()) == overview["stand_ins"]["rows"]
    data = json.loads((ROOT / "docs" / "reference" / "references-data.json").read_text(encoding="utf-8"))
    assert overview["papers"]["entries"] == len(data["entries"])
    assert overview["papers"]["with_doi"] + overview["papers"]["without_doi"] == overview["papers"]["entries"]
    assert overview["data"] == {"downloaded_at_setup": len(data["setup_datasets"]), "archives": len(data["archives"])}


def test_the_verdict_adds_up_and_the_shares_round_half_up():
    verdict = builder.build_overview()["verdict"]
    total = builder.build_overview()["stand_ins"]["rows"]
    kinds = ("read_from_a_source", "unverified", "educated_guess_from_physics", "chosen_for_the_output")
    assert sum(verdict[kind]["rows"] for kind in kinds) == total
    assert 98 <= sum(verdict[kind]["percent"] for kind in kinds) <= 102
    assert builder.share(1, 2) == 50 and builder.share(1, 3) == 33 and builder.share(2, 3) == 67 and builder.share(0, 0) == 0
    assert builder.share(1, 8) == 13 and builder.share(7, 8) == 88


def test_the_tiers_and_tweaks_are_the_registrys_legends():
    overview = builder.build_overview()
    assert [row[0] for row in overview["tiers"]] == [f"T{n}" for n in range(1, 10)]
    assert [row[0] for row in overview["tweaks"]] == [f"L{n}" for n in range(5)]


# --- The page --------------------------------------------------------------------------------------------------------------------

def test_the_page_renders_with_the_banner_the_tinted_header_and_the_legends(client):
    html = _html(client, "/references")
    assert BANNER_ID in html and "Out of universe." in html
    assert '<body class="app references-mode">' in html
    assert "<title>References" in html
    for heading in ("Where the numbers stand", "How sure the stand-ins are", "What the numbers are read from", "The source tiers", "The tweak levels"):
        assert f">{heading}</h2>" in html, heading
    for code in [f"T{n}" for n in range(1, 10)] + [f"L{n}" for n in range(5)]:
        assert f'<th scope="row">{code}</th>' in html, code
    assert html.count(BANNER_ID) == 1


def test_the_other_pages_have_no_banner_and_no_tint(client):
    html = _html(client, "/help/life")
    assert BANNER_ID not in html and "references-mode" not in html


def test_the_help_menu_ends_with_references_after_a_divider(client):
    html = _html(client, "/classes")
    menu = re.search(r'<sl-dropdown class="site-help-menu".*?</sl-dropdown>', html, re.S).group(0)
    assert re.findall(r'data-href="([^"]*)"[^>]*>([^<]*)<', menu)[-1] == ("/references", "References")
    assert menu.index("<sl-divider>") < menu.index('data-href="/references"')
    assert menu.count("<sl-divider>") == 1
    page = _html(client, "/references")
    assert 'data-href="/references" aria-current="page">References</sl-menu-item>' in page
    assert 'class="site-help-trigger" aria-haspopup="menu" aria-current="page"' in page


def test_the_route_and_the_address_scheme_agree(app):
    with app.test_request_context("/"):
        from planetgen.web.helpers import page_url
        assert page_url("references_page") == reference_urls.reference_path("overview") == "/references"


def test_the_sentences_keep_the_guide_and_may_use_the_out_of_universe_words():
    page = references_site.overview(references_site.load("overview.json"))
    texts = list(page["intro"]) + [p for section in page["sections"] for p in section["paragraphs"]]
    assert len(texts) >= 8
    broken = [(sentence, rule) for text in texts for sentence in plainwords.sentences(text)
              for rule, _detail in plainwords.check_sentence(sentence, references=True)]
    assert not broken, broken
    assert plainwords.exempt("/references") and any(plainwords.out_of_universe(text) for text in texts)
    labels = [head for section in page["sections"] if section["table"] for head in section["table"]["head"]]
    assert all(not plainwords.check_label(label) for label in labels), labels


def test_a_sentence_with_a_process_word_breaks_the_guide_off_references_and_not_on_it():
    sentence = "Each stand-in has a tweak level."
    assert {rule for rule, _detail in plainwords.check_sentence(sentence)} == {"banned"}
    assert plainwords.check_sentence(sentence, references=True) == []
    assert plainwords.check_sentence("It runs in order to start.", references=True)


def test_a_site_that_has_not_built_the_data_says_so_and_keeps_the_banner(client, monkeypatch):
    monkeypatch.setattr(references_site, "load", lambda name, folder=None: None)
    html = _html(client, "/references")
    assert BANNER_ID in html and "The pages are not built" in html and "build_reference_pages.py" in html
    assert 'id="standing"' not in html


def test_the_data_is_read_again_when_the_file_changes(tmp_path):
    (tmp_path / "overview.json").write_text('{"a": 1}', encoding="utf-8")
    assert references_site.load("overview.json", str(tmp_path)) == {"a": 1}
    path = tmp_path / "overview.json"
    path.write_text('{"a": 2}', encoding="utf-8")
    os.utime(path, ns=(path.stat().st_atime_ns, path.stat().st_mtime_ns + 5_000_000_000))
    assert references_site.load("overview.json", str(tmp_path)) == {"a": 2}
    assert references_site.load("nothing.json", str(tmp_path)) is None


# --- The setup step --------------------------------------------------------------------------------------------------------------

def _deploy_common_shell(python, script):
    code = f'set -euo pipefail\nSCRIPT_DIR="{ROOT}"\nPYTHON="{python}"\nsource "$SCRIPT_DIR/scripts/deploy-common.sh"\n{script}\n'
    return subprocess.run(["bash", "-c", code], capture_output=True, text=True, timeout=120)


def test_update_sh_builds_the_pages_and_a_current_checkout_is_left_as_it_is():
    before = (ROOT / "docs" / "references" / "overview.json").read_bytes()
    result = _deploy_common_shell(sys.executable, 'build_reference_pages; echo "failed=$REFERENCE_PAGES_FAILED"')
    assert result.returncode == 0, result.stderr
    assert "failed=0" in result.stdout and "Building the References pages" in result.stdout
    assert (ROOT / "docs" / "references" / "overview.json").read_bytes() == before


def test_update_sh_build_failure_only_warns(tmp_path):
    fake = tmp_path / "python"
    fake.write_text("#!/bin/sh\necho 'boom'\nexit 1\n")
    fake.chmod(0o755)
    result = _deploy_common_shell(str(fake), 'build_reference_pages; echo "failed=$REFERENCE_PAGES_FAILED"')
    assert result.returncode == 0, result.stderr
    assert "failed=1" in result.stdout and "warning: the References pages could not be built" in result.stderr


def test_update_script_builds_the_pages_in_step_3_and_repeats_a_failure_at_the_end():
    text = (ROOT / "update.sh").read_text(encoding="utf-8")
    step3 = text.split('UPDATE_STEP="step 3/9')[1].split('UPDATE_STEP="step 4/9')[0]
    assert "build_reference_pages" in step3
    assert 'REFERENCE_PAGES_FAILED' in text.split("DATASETS_MISSING")[-1]
