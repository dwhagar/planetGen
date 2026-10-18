"""
The docs index (DOC.52, `scripts/build_docs_index.py`, `docs/index.html`): every HTML page under `docs/` is listed
with a description and links back to the index, and the index is current.
"""

import importlib.util
import os
import re

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_spec = importlib.util.spec_from_file_location("build_docs_index", os.path.join(ROOT, "scripts", "build_docs_index.py"))
index = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(index)


def _text():
    with open(index.INDEX, encoding="utf-8") as handle:
        return handle.read()


def test_the_index_is_up_to_date():
    assert _text() == index.build(), "docs/index.html is out of date: run python scripts/build_docs_index.py"


def test_every_html_page_under_docs_is_listed_with_a_description():
    text = _text()
    pages = index.html_pages()
    assert len(pages) >= 7
    for path, title, description in pages:
        assert f'href="{path}"' in text, f"{path} is not in the docs index: run python scripts/build_docs_index.py"
        assert description, f"{path} has no <meta name=\"description\" content=\"...\"> (what it contains and who it is for)"
    assert "design/architecture-map.html" in text and "plan/version-plan.html" in text


def test_every_page_links_back_to_the_index():
    for path, _title, _description in index.html_pages():
        depth = path.count("/")
        with open(os.path.join(index.DOCS, path), encoding="utf-8") as handle:
            page = handle.read()
        assert f'href="{"../" * depth}index.html"' in page, f"{path} has no link back to the docs index"


def test_the_index_links_resolve_to_files_and_the_goal_issues():
    text = _text()
    for href in re.findall(r'<a href="([^":#]+)"', text):
        assert os.path.exists(os.path.join(index.DOCS, href)), href
    assert index.ISSUES_URL in text
    assert 'href="plan/issue-groups.md"' in text and 'href="design/README.md"' in text


def test_a_page_without_a_description_is_noticed(tmp_path, monkeypatch):
    docs = tmp_path / "docs"
    (docs / "plan").mkdir(parents=True)
    (docs / "plan" / "new.html").write_text("<!doctype html><title>New page</title>", encoding="utf-8")
    monkeypatch.setattr(index, "DOCS", str(docs))
    assert index.html_pages() == [("plan/new.html", "New page", "")]
