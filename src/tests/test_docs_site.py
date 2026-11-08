"""
The documentation website (`scripts/build_docs_site.py`, `docs/site/`): Markdown becomes pages for GitHub Pages with
working links and anchors, subscripted chemical formulas outside code, and math left for KaTeX; and the pages are current.
"""

import importlib.util
import os
import re

import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_spec = importlib.util.spec_from_file_location("build_docs_site", os.path.join(ROOT, "scripts", "build_docs_site.py"))
site = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(site)


@pytest.mark.parametrize("formula", ["H2O", "CO2", "CH4", "H2SO4", "Fe2O3", "Mg(ClO4)2", "O2", "NH4SH"])
def test_chemical_formulas_get_subscripts(formula):
    assert site.is_formula(formula)
    out = site.decorate(f"<p>the {formula} here</p>")
    assert "<sub>" in out and re.sub(r"</?sub>", "", out) == f"<p>the {formula} here</p>"


@pytest.mark.parametrize("label", ["GEN1377", "F5V", "BC13", "HP22", "C4", "R5203", "IC348"])
def test_labels_that_look_like_formulas_are_left_alone(label):
    assert not site.is_formula(label)


def test_the_water_formula_is_subscripted_in_prose_but_not_in_code_or_math():
    out = site.decorate('<p>H2O and <code>H2O</code> and <span class="math">\\(H2O\\)</span></p>')
    assert out.count("H<sub>2</sub>O") == 1


def test_math_is_protected_from_markdown_and_marked_for_katex():
    text, found = site.protect_math("Distance $x_c - x_p$ costs $5 and $10.\n\n$$a_1 + b_2$$\n\n`$y_1$`\n")
    assert sorted(tex for _display, tex in found) == ["a_1 + b_2", "x_c - x_p"]
    assert "`$y_1$`" in text and "costs $5 and $10" in text
    inline = next(i for i, (display, _tex) in enumerate(found) if not display)
    assert site.restore_math(f"\ue000{inline}\ue001", found) == '<span class="math">\\(x_c - x_p\\)</span>'


def test_headings_get_the_github_anchor():
    seen = {}
    assert site.slug("What (is) this? Part 2", seen) == "what-is-this-part-2"
    assert site.slug("What (is) this? Part 2", seen) == "what-is-this-part-2-1"


def test_the_site_has_a_page_for_every_markdown_file_and_every_link_resolves():
    files = site.build_files()
    pages = site.discover()
    assert len(pages) > 300 and all(page.out in files for page in pages)
    assert "index.html" in files and "nav.json" in files
    for path, text in files.items():
        if not path.endswith(".html"):
            continue
        body = re.sub(r"<script>.*?</script>", "", text, flags=re.S)
        for href in re.findall(r'<a [^>]*href="([^"]+)"', body):
            if href.startswith(("http", "mailto:", "#")):
                continue
            target = os.path.normpath(os.path.join(os.path.dirname(path), href.partition("#")[0]))
            if target.startswith(".."):
                real = os.path.normpath(os.path.join(site.SITE, target))
                assert os.path.exists(real) or href.endswith("webconfig.md"), f"{path}: {href}"
            elif "webconfig.md" not in href:  # CHANGELOG.md names a file that was removed
                assert target in files, f"{path}: {href}"


@pytest.mark.globaldocs
def test_the_site_is_up_to_date():
    assert site.main(["--check"]) == 0
