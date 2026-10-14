# tests/test_no_todo_codes_shown.py

"""
No TODO code (`PERF.67`, `NAV.42`, ...) appears where a user can see it
(UX.93): in page templates, static files, the strings the code shows
(messages, settings descriptions, CLI output) or the CLI `--help` text.
Code comments, docstrings, the bare strings that document a constant,
`TODO.md`, `changes/` notes, the changelog and the design docs may keep
them. A failure names the file and line so the code can be cut from the
text.
"""

import argparse
import ast
import importlib
import re
import sys
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[1]
CATEGORIES = "ADM|API|DB|DOC|GEN|MAP|NAV|OPS|PERF|POP|SEC|TEST|USR|UX|VIEW"
TODO_CODE = re.compile(rf"\b(?:{CATEGORIES})\.\d+\b")
"""re.Pattern: A category from the TODO table, a dot and a number."""

ALLOWED = {}
"""dict: `{path relative to src: reason}` for a file where a code is wanted. None are."""


def _shown(path):
    return [p for p in path if "vendor" not in p.parts and "__pycache__" not in p.parts
            and str(p.relative_to(SRC)) not in ALLOWED]


def _hits(text, name):
    return [f"{name}:{text.count(chr(10), 0, m.start()) + 1}: {m.group(0)}" for m in TODO_CODE.finditer(text)]


def test_templates_show_no_todo_code():
    found = []
    for path in _shown(sorted((SRC / "planetgen" / "web" / "templates").rglob("*.html"))):
        text = path.read_text(encoding="utf-8")
        # Comments may keep codes; blank them (keeping line breaks) so line numbers hold.
        text = re.sub(r"\{#.*?#\}|<!--.*?-->", lambda m: re.sub(r"[^\n]", " ", m.group(0)), text, flags=re.S)
        found += _hits(text, str(path.relative_to(SRC)))
    assert not found, "TODO codes in template text:\n" + "\n".join(found)


def test_static_files_show_no_todo_code():
    found = []
    for path in _shown(sorted((SRC / "html").rglob("*"))):
        if path.suffix not in (".html", ".js", ".mjs", ".css", ".svg", ".txt", ".json", ".md"):
            continue
        text = path.read_text(encoding="utf-8", errors="ignore")
        if path.suffix in (".js", ".mjs", ".css"):
            text = re.sub(r"/\*.*?\*/", lambda m: re.sub(r"[^\n]", " ", m.group(0)), text, flags=re.S)
            text = re.sub(r"(?m)(^|[\s;{}),])//.*$", r"\1", text)
        if path.suffix in (".svg", ".html"):
            text = re.sub(r"<!--.*?-->", lambda m: re.sub(r"[^\n]", " ", m.group(0)), text, flags=re.S)
        found += _hits(text, str(path.relative_to(SRC)))
    assert not found, "TODO codes in static files:\n" + "\n".join(found)


def test_strings_the_code_shows_have_no_todo_code():
    """Every string literal in `planetgen` except docstrings and the bare
    string under a constant that documents it (comments are not strings)."""
    found = []
    for path in _shown(sorted((SRC / "planetgen").rglob("*.py"))):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        documentation = {id(n.value) for n in ast.walk(tree)
                         if isinstance(n, ast.Expr) and isinstance(n.value, ast.Constant)}
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and isinstance(node.value, str) and id(node) not in documentation \
                    and TODO_CODE.search(node.value):
                found.append(f"{path.relative_to(SRC)}:{node.lineno}: {TODO_CODE.search(node.value).group(0)}")
    assert not found, "TODO codes in strings the code shows:\n" + "\n".join(found)


COMMAND_LINE_PROGRAMS = ("generate", "query", "lockouts", "config", "worker", "orbits", "migrate", "reset", "dedupe",
                         "benchmark", "version_history")


@pytest.mark.parametrize("program", COMMAND_LINE_PROGRAMS)
def test_command_line_help_shows_no_todo_code(program, monkeypatch, capsys):
    """Each program's `--help` (and every subcommand's), as argparse prints it."""
    module = importlib.import_module(f"planetgen.cli.{program}")
    if not hasattr(module, "main"):
        pytest.skip(f"planetgen.cli.{program} has no main()")
    made = []
    real_init = argparse.ArgumentParser.__init__

    def init(self, *args, **kwargs):
        real_init(self, *args, **kwargs)
        made.append(self)

    monkeypatch.setattr(argparse.ArgumentParser, "__init__", init)
    monkeypatch.setattr(sys, "argv", [program, "--help"])
    try:
        module.main()
    except SystemExit:
        pass
    capsys.readouterr()
    found = []
    for parser in made:
        text = parser.format_help()
        found += [f"{parser.prog}: {m.group(0)}" for m in TODO_CODE.finditer(text)]
        found += [f"{parser.prog}: {name}" for name in re.findall(r"\b[A-Z][A-Z0-9]*(?:_[A-Z0-9]+){2,}\b", text)
                  if name.endswith("_THRESHOLD") or name.endswith("_LIMIT")]
    assert made, f"planetgen.cli.{program} built no parser"
    assert not found, "TODO codes or constant names in --help:\n" + "\n".join(found)
