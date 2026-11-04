"""
DOC.102: the whole site is in American spelling (Boss, 10:33 AM EDT 2026-10-31). `util/plainwords.british_spellings` reads the words the
program shows: every string of the program's code that is a sentence or a piece of one (a docstring, a comment and a one-word name,
key, or stored value are not shown to a visitor), the text of every template, and the scripts and styles of the pages.
"""

import ast
import os
import re

from planetgen.util import plainwords

PACKAGE = os.path.dirname(os.path.dirname(os.path.abspath(plainwords.__file__)))

SKIPPED_FILES = ("util/plainwords.py", "db/models.py", "db/migrations/")
"""tuple: The files not read: the checker's own word list, the table definitions (their CHECK lists hold stored values) and the migrations."""

JINJA = re.compile(r"\{#.*?#\}|\{%.*?%\}|\{\{.*?\}\}", re.DOTALL)
"""re.Pattern: A Jinja comment, tag or expression: code, not the text a template shows."""


MARKUP_NAME = re.compile(r"data-[a-z-]+|[\w./-]*\.md")
"""re.Pattern: A data attribute name or a file name, which is code: it keeps the spelling the other code reads it by."""


def _skipped(path):
    relative = os.path.relpath(path, PACKAGE).replace(os.sep, "/")
    return any(relative.startswith(skipped) or relative == skipped for skipped in SKIPPED_FILES)


def _files(extensions):
    for folder, _dirs, names in os.walk(PACKAGE):
        if "__pycache__" in folder:
            continue
        for name in sorted(names):
            path = os.path.join(folder, name)
            if name.endswith(extensions) and not _skipped(path):
                yield path


def _shown_strings(path):
    """`[(line, string)]` of the strings of a Python file that read as text: not a docstring, and not a single word outside an f-string."""
    with open(path, encoding="utf-8") as handle:
        tree = ast.parse(handle.read())
    statements, pieces = set(), set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant) and isinstance(node.value.value, str):
            statements.add(id(node.value))
        if isinstance(node, ast.JoinedStr):
            pieces.update(id(value) for value in node.values)
    return [(node.lineno, node.value) for node in ast.walk(tree)
            if isinstance(node, ast.Constant) and isinstance(node.value, str) and id(node) not in statements
            and (" " in node.value.strip() or id(node) in pieces)]


def test_the_strings_of_the_code_are_in_american_spelling():
    found = []
    for path in _files((".py",)):
        for line, text in _shown_strings(path):
            for _rule, detail in plainwords.british_spellings(MARKUP_NAME.sub(" ", text)):
                found.append(f"{os.path.relpath(path, PACKAGE)}:{line}: {detail}")
    assert not found, "British spelling in text the program shows (DOC.102):\n  " + "\n  ".join(sorted(found))


def test_the_templates_and_pages_are_in_american_spelling():
    found = []
    for path in _files((".html", ".js", ".css", ".txt", ".md")):
        with open(path, encoding="utf-8", errors="ignore") as handle:
            text = JINJA.sub(lambda match: re.sub(r"[^\n]", " ", match.group(0)), handle.read())    # the lines stay where they were
        for number, line in enumerate(text.split("\n"), 1):
            for _rule, detail in plainwords.british_spellings(line):
                found.append(f"{os.path.relpath(path, PACKAGE)}:{number}: {detail}")
    assert not found, "British spelling in a template or page (DOC.102):\n  " + "\n  ".join(sorted(found))


def test_the_word_list_knows_the_spellings_the_sweep_found():
    for british in ("the cancelled run", "a disc of gas", "the catalogue", "moves towards", "a modelled value", "a labelling pass"):
        assert plainwords.british_spellings(british), british
    assert not plainwords.british_spellings("a canceled run in a disk of gas")
