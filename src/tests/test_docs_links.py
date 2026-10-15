# tests/test_docs_links.py

"""
DOC.29: every relative link and anchor in the project's Markdown documents resolves.

Broken relative paths and 11 dead links in docs/design/architecture.md went unnoticed until the 2026-10-14 sweep.
The check covers README.md, INSTALL.md, CHANGELOG.md's links and every `.md` file under `docs/`; it reads the
files, so it needs no database or browser. Links to web addresses are not followed, and links inside code (a
fence, or `backticks`) are not links.
"""

import os
import re
from pathlib import Path
from urllib.parse import unquote

import pytest

ROOT = Path(__file__).resolve().parents[2]
FILES = sorted({ROOT / "README.md", ROOT / "INSTALL.md", *ROOT.glob("docs/**/*.md")})

LINK = re.compile(r"(?<!\!)\[(?:[^\]\[]|\[[^\]]*\])*\]\(\s*<?([^)\s>]+)>?(?:\s+\"[^\"]*\")?\s*\)")
HEADING = re.compile(r"^ {0,3}(#{1,6})\s+(.*?)\s*#*\s*$")
EXPLICIT_ANCHOR = re.compile(r"""<a\s+(?:name|id)=["']([^"']+)["']""")


def _prose(text):
    """`text` without fenced code blocks and inline code (what a reader sees as links and headings)."""
    lines, fence = [], None
    for line in text.splitlines():
        stripped = line.lstrip()
        marker = stripped[:3]
        if fence is None and marker in ("```", "~~~"):
            fence = marker
            lines.append("")
            continue
        if fence is not None:
            if stripped.startswith(fence):
                fence = None
            lines.append("")
            continue
        lines.append(re.sub(r"`[^`]*`", "", line))
    return "\n".join(lines)


def _slug(heading):
    """GitHub's anchor for a heading: lower case, formatting and punctuation dropped, spaces to hyphens."""
    heading = re.sub(r"!?\[([^\]]*)\]\([^)]*\)", r"\1", heading)
    heading = re.sub(r"[*~]", "", heading).strip().lower()
    heading = re.sub(r"[^\w\- ]", "", heading)
    return heading.replace(" ", "-")


_anchor_cache = {}


def _anchors(path):
    if path not in _anchor_cache:
        text = path.read_text(encoding="utf-8")
        seen, found = {}, set(EXPLICIT_ANCHOR.findall(text))
        raw = []
        fence = None
        for line in text.splitlines():
            stripped = line.lstrip()
            if fence is None and stripped[:3] in ("```", "~~~"):
                fence = stripped[:3]
                continue
            if fence is not None:
                if stripped.startswith(fence):
                    fence = None
                continue
            match = HEADING.match(line)
            if match:
                raw.append(match.group(2).replace("`", ""))
        for heading in raw:
            slug = _slug(heading)
            count = seen.get(slug, 0)
            seen[slug] = count + 1
            found.add(slug if count == 0 else f"{slug}-{count}")
        _anchor_cache[path] = found
    return _anchor_cache[path]


def broken_links(path):
    """`(line number, link, why)` for each relative link or anchor of `path` that does not resolve."""
    problems = []
    for number, line in enumerate(_prose(path.read_text(encoding="utf-8")).splitlines(), 1):
        for target in LINK.findall(line):
            if re.match(r"[a-zA-Z][a-zA-Z0-9+.-]*:", target) or target.startswith("//"):
                continue  # http, https, mailto
            name, _, anchor = target.partition("#")
            name = unquote(name)
            if name:
                resolved = (path.parent / name).resolve()
                if not resolved.exists():
                    problems.append((number, target, "no such file"))
                    continue
            else:
                resolved = path
            if anchor and resolved.suffix == ".md" and unquote(anchor).lower() not in _anchors(resolved):
                problems.append((number, target, "no such heading"))
    return problems


@pytest.mark.parametrize("path", FILES, ids=lambda p: os.path.relpath(p, ROOT))
def test_every_relative_link_and_anchor_resolves(path):
    problems = broken_links(path)
    assert not problems, "\n".join(f"{os.path.relpath(path, ROOT)}:{n}: {link} ({why})" for n, link, why in problems)


def test_the_check_finds_a_dead_path_a_dead_anchor_and_ignores_code(tmp_path):
    other = tmp_path / "other.md"
    other.write_text("# Real heading\n\n## Two words, with punctuation!\n", encoding="utf-8")
    page = tmp_path / "page.md"
    page.write_text(
        "[ok](other.md) [ok anchor](other.md#real-heading) [ok punctuation](other.md#two-words-with-punctuation) "
        "[same file](#top)\n# Top\n\n"
        "[dead path](missing.md) [dead anchor](other.md#nowhere) [web](https://example.com/x.md)\n\n"
        "`[in code](missing.md)`\n\n```\n[fenced](missing.md)\n```\n", encoding="utf-8")
    found = {(link, why) for _n, link, why in broken_links(page)}
    assert found == {("missing.md", "no such file"), ("other.md#nowhere", "no such heading")}
