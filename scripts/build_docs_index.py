#!/usr/bin/env python3
"""
Builds the docs index, `docs/index.html` (DOC.52).

Boss (2026-10-15, 7:01 PM EDT): "Add a TODO item in DOCS to build an index with links to pages with descriptions of
what each page contains under docs to link the HTML together." The HTML pages under `docs/` (the plans, the TODO
pages, the architecture pages and any added later) otherwise link to each other only by hand.

How it works
============
The index lists every `*.html` page under `docs/` with its `<title>` and its `<meta name="description">` (one or two
sentences: what the page contains and who it is for), so a new page cannot be forgotten: a page with no description,
or one missing from the index, fails `src/tests/test_docs_index.py`. It also links the reference notes (`docs/*.md`,
each with its first heading and first line), the design notes' own index (`docs/design/README.md`), the plan notes
and the GitHub issue list for the goals (`docs/plan/issue-groups.md` names the issues). Every link is relative, so
the index and the pages open offline. Each generated page links back to the index from its header.

Usage
=====
    python scripts/build_docs_index.py            # write docs/index.html
    python scripts/build_docs_index.py --check    # fail if docs/index.html is out of date
"""

import argparse
import html
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCS = os.path.join(ROOT, "docs")
INDEX = os.path.join(DOCS, "index.html")
ISSUES_URL = "https://github.com/dwhagar/planetGen/issues"

TITLE_RE = re.compile(r"<title>(.*?)</title>", re.S)
DESCRIPTION_RE = re.compile(r'<meta\s+name="description"\s+content="([^"]*)"', re.S)
BACK_LINK = '<a href="{}index.html">Docs index</a>'.format

SECTIONS = (
    ("generated", "Generated pages: plans, the TODO and the code map",
     "Rebuilt by scripts from docs/TODO.md, the plan notes and the code; never edited by hand."),
    ("design", "Architecture overview", "Written by hand: how the program fits together."),
)
"""tuple: `(folder under docs/, heading, lead)`: how the HTML pages are grouped; a folder not named here is listed
last under its own name."""


def html_pages():
    """`[(path relative to docs/, title, description)]` of every HTML page under `docs/` but the index itself."""
    found = []
    for folder, dirs, files in os.walk(DOCS):
        dirs[:] = sorted(d for d in dirs if d != "__pycache__")
        for name in sorted(files):
            if not name.endswith(".html"):
                continue
            path = os.path.relpath(os.path.join(folder, name), DOCS).replace(os.sep, "/")
            if path == "index.html":
                continue
            with open(os.path.join(folder, name), encoding="utf-8") as handle:
                text = handle.read()
            title = TITLE_RE.search(text)
            description = DESCRIPTION_RE.search(text)
            found.append((path, html.unescape(title.group(1).strip()) if title else path,
                          html.unescape(description.group(1).strip()) if description else ""))
    return found


def markdown_notes(folder):
    """`[(file, heading, first line)]` of the Markdown files directly in `docs/<folder>`."""
    base = os.path.join(DOCS, folder) if folder else DOCS
    found = []
    for name in sorted(os.listdir(base)):
        if not name.endswith(".md") or name == "TODO.md":
            continue
        with open(os.path.join(base, name), encoding="utf-8") as handle:
            lines = [line.strip() for line in handle.read().splitlines()]
        heading = next((line[2:].strip() for line in lines if line.startswith("# ")), name)
        paragraph, started = [], False
        for line in lines[lines.index(f"# {heading}") + 1 if f"# {heading}" in lines else 0:]:
            if not line:
                if started:
                    break
                continue
            if line.startswith(("#", "|", "```", "-", "*", ">", "<", "[", "!")) and not started:
                continue
            started = True
            paragraph.append(line)
        text = re.sub(r"`", "", re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", " ".join(paragraph)))
        sentence = re.match(r"(.+?[.!?])(\s|$)", text)
        first = sentence.group(1) if sentence else text
        found.append((f"{folder + '/' if folder else ''}{name}", heading, first if len(first) <= 240 else first[:237] + "..."))
    return found


def _entry(href, title, text):
    return (f'<li><a href="{html.escape(href)}">{html.escape(title, quote=False)}</a>'
            f'<span class="file"> {html.escape(href)}</span>'
            f'{"<p>" + html.escape(text, quote=False) + "</p>" if text else ""}</li>')


def build():
    """The index page's text."""
    pages = html_pages()
    grouped = {}
    for path, title, description in pages:
        grouped.setdefault(path.split("/")[0] if "/" in path else "", []).append((path, title, description))
    order = [key for key, _heading, _lead in SECTIONS] + sorted(k for k in grouped if k not in dict((s[0], 1) for s in SECTIONS))
    headings = {key: (heading, lead) for key, heading, lead in SECTIONS}
    body = []
    for key in order:
        if key not in grouped:
            continue
        heading, lead = headings.get(key, (key or "Pages", ""))
        body.append(f"<h2>{html.escape(heading)}</h2>" + (f'<p class="lead">{html.escape(lead)}</p>' if lead else ""))
        body.append("<ul>" + "".join(_entry(path, title, description) for path, title, description in grouped[key]) + "</ul>")
    reference = "<ul>" + "".join(_entry(f, h, t) for f, h, t in markdown_notes("reference")) + "</ul>"
    design = ("<ul>" + _entry("README.md", "Docs overview", "One line for every folder and key file under docs/.") +
              _entry("design/README.md", "Design notes", "The index of every design note, with its status.") +
              _entry("plan/issue-groups.md", "Goal issues",
                     "The GitHub issues that group the TODO items into goals, and the items each covers.") +
              _entry("TODO.md", "TODO", "Every open item, in full.") + "</ul>")
    return (TEMPLATE.replace("/*PAGES*/", "\n".join(body)).replace("/*REFERENCE*/", reference)
            .replace("/*NOTES*/", design).replace("/*ISSUES*/", ISSUES_URL))


TEMPLATE = """<!doctype html>
<html lang="en">
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>planetGen Docs Index</title>
<meta name="description" content="Every HTML page under docs/ with what it contains, the reference notes, the design notes and the goal issues; start here to find a page.">
<style>
/* Generated by scripts/build_docs_index.py (DOC.52): do not edit by hand. */
:root{--bg:#F3F5F8;--surface:#fff;--ink:#152030;--muted:#566275;--line:#D8DEE7;--link:#1F5FA8}
@media (prefers-color-scheme: dark){:root{--bg:#0C1118;--surface:#141B25;--ink:#E5EAF1;--muted:#97A3B5;--line:#273244;--link:#7DB4F0;color-scheme:dark}}
body{margin:0;background:var(--bg);color:var(--ink);font:16px/1.55 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif}
main{max-width:860px;margin:0 auto;padding:24px 16px 60px}
h1{font-size:26px;margin:0 0 4px} h2{font-size:19px;margin:30px 0 4px;border-bottom:1px solid var(--line);padding-bottom:4px}
a{color:var(--link);text-decoration:none} a:hover{text-decoration:underline}
.lead,.muted{color:var(--muted)} .lead{margin:2px 0 8px}
ul{list-style:none;padding:0;margin:8px 0} li{background:var(--surface);border:1px solid var(--line);border-radius:8px;padding:9px 13px;margin:7px 0}
li p{margin:3px 0 0;color:var(--muted);font-size:14.5px} .file{font:12px ui-monospace,Menlo,Consolas,monospace;color:var(--muted)}
</style>
<main>
<h1>planetGen docs index</h1>
<p class="lead">Where each page under <code>docs/</code> is and what it holds. Generated by <code>scripts/build_docs_index.py</code>; the goals are on <a href="/*ISSUES*/">GitHub</a>.</p>
/*PAGES*/
<h2>Reference notes</h2>
<p class="lead">How to install, run, configure and use the program.</p>
/*REFERENCE*/
<h2>Design notes, plans and goals</h2>
/*NOTES*/
</main>
"""


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--check", action="store_true", help="fail if docs/index.html is out of date")
    args = parser.parse_args(argv)
    text = build()
    if args.check:
        current = open(INDEX, encoding="utf-8").read() if os.path.exists(INDEX) else None
        if current != text:
            print("docs/index.html is out of date: run python scripts/build_docs_index.py")
            return 1
        print("docs index is up to date")
        return 0
    with open(INDEX, "w", encoding="utf-8") as handle:
        handle.write(text)
    print(f"wrote docs/index.html: {len(html_pages())} pages")
    return 0


if __name__ == "__main__":
    sys.exit(main())
