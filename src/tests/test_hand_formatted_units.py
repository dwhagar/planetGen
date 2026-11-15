"""
UX.318 (Boss, 2026-11-14): every rate, duration, size and memory figure the generator logs or shows comes from the units code
(`planetgen.util.format`: `format_duration_seconds`, `format_bytes`, `format_megabytes`, `format_rate`, `format_byte_rate`),
so a log line never says "0.0 sectors/s" where "2.4 s a sector" would help, and memory is never bytes in one line, MB in
the next and GB in a third. This test fails on a hand-formatted unit: an f-string whose formatted value is followed by a
unit word (`f"{seconds:.1f} s"`, `f"{mb:,.0f} MB"`, `f"{rate:.2f} per s"`), and a template whose `{{ ... }}` is followed by one.

The units code itself and the generated browser mirrors are the only places allowed to spell a unit after a number. A few
exact counts stay as words on purpose and are listed in `ALLOWED`: a byte length that must read exactly (a seed is 32 bytes,
a row of 70,000 bytes against a packet budget) and a described fact ("Hours outside a year").
"""

import ast
import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1] / "planetgen"

UNITS = ("B", "KB", "kB", "MB", "GB", "TB", "bytes", "byte", "ms", "µs", "s", "sec", "secs", "second", "seconds", "min",
         "mins", "minute", "minutes", "h", "hr", "hrs", "hour", "hours", "day", "days", "/s", "per s", "per sec", "per second",
         "a second", "MB/s", "KB/s", "kB/s")
UNIT_START = re.compile(r"^ ?(" + "|".join(re.escape(unit) for unit in sorted(UNITS, key=len, reverse=True)) + r")(?![\w/])")

OWN_CODE = {"util/format.py", "util/ladders.py", "physics/unit_registry.py"}
"""The units code: the one place a number may be followed by its unit."""

ALLOWED = {
    ("galaxy/seed.py", "bytes"), ("galaxy/object_uid.py", "bytes"),                 # exact lengths of a key
    ("db/store.py", "bytes"), ("db/connection_loss.py", "bytes"), ("backup/parts.py", "bytes"), ("datasets/fetch.py", "bytes"),  # exact byte counts against a budget or a size on disk
    ("describe/build.py", "hours"),                                                  # a described fact, "Hours outside a year"
    ("web/templates/partials/day_counter.html", "hours"),                            # the solar day of a body, a described fact
    ("web/templates/slow_loading.html", "seconds"), ("web/templates/galaxy_generating.html", "seconds"),
    ("web/templates/admin_update.html", "seconds"), ("web/templates/admin_update.html", "minutes"),
    ("web/templates/maintenance.html", "minutes"),                                   # prose about a page's own refresh, a whole number of seconds
}
"""set[tuple[str, str]]: (file under planetgen/, unit word) pairs that are counts in words on purpose, not sizes or durations to read."""

SQL_KEYWORD = re.compile(r"^ ?s (WHERE|ORDER|GROUP|LEFT|JOIN|ON|INNER|AND)\b")
"""A table alias `s` in an SQL f-string is not a second."""


def _python_files():
    return sorted(path for path in ROOT.rglob("*.py") if "migrations" not in path.parts)


def _templates():
    return sorted((ROOT / "web" / "templates").rglob("*.html"))


def _unit_after_value(text):
    """The unit word `text` (the literal that follows a formatted value) starts with, or None."""
    if SQL_KEYWORD.match(text):
        return None
    found = UNIT_START.match(text)
    if found is None or (found.group(1) == "s" and not text.startswith(" ")):
        return None  # `{noun}s` is a plural
    return found.group(1)


def _hand_formatted_in_python(path):
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in ast.walk(tree):
        if not isinstance(node, ast.JoinedStr):
            continue
        values = node.values
        for index, part in enumerate(values[:-1]):
            following = values[index + 1]
            if isinstance(part, ast.FormattedValue) and isinstance(following, ast.Constant) and isinstance(following.value, str):
                unit = _unit_after_value(following.value)
                if unit:
                    yield node.lineno, unit, following.value[:40]


TEMPLATE_VALUE = re.compile(r"\}\}( ?)(" + "|".join(re.escape(unit) for unit in sorted(UNITS, key=len, reverse=True)) + r")(?![\w/])")


def _hand_formatted_in_template(path):
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        for found in TEMPLATE_VALUE.finditer(line):
            if found.group(2) == "s" and not found.group(1):
                continue  # `{{ noun }}s` is a plural
            yield number, found.group(2), line.strip()[:60]


@pytest.mark.parametrize("path", _python_files(), ids=lambda path: str(path.relative_to(ROOT)))
def test_no_python_module_formats_a_unit_by_hand(path):
    where = str(path.relative_to(ROOT)).replace("\\", "/")
    if where in OWN_CODE:
        pytest.skip("the units code spells units")
    found = [(line, unit, text) for line, unit, text in _hand_formatted_in_python(path) if (where, unit) not in ALLOWED]
    assert not found, (f"{where} formats a unit by hand at " + "; ".join(f"line {line}: {{...}}{text!r}" for line, _, text in found)
                       + ". Use format_duration_seconds, format_bytes, format_megabytes, format_rate or format_byte_rate from planetgen.util.format (UX.318).")


@pytest.mark.parametrize("path", _templates(), ids=lambda path: str(path.relative_to(ROOT)))
def test_no_template_formats_a_unit_by_hand(path):
    where = str(path.relative_to(ROOT)).replace("\\", "/")
    found = [(line, unit, text) for line, unit, text in _hand_formatted_in_template(path) if (where, unit) not in ALLOWED]
    assert not found, (f"{where} formats a unit by hand at " + "; ".join(f"line {line}: {text!r}" for line, _, text in found)
                       + ". Pass the seconds or bytes to duration_text() or bytes_text() instead (UX.318).")


def test_the_guard_sees_a_hand_formatted_unit(tmp_path):
    """The scanner itself: it catches the forms the sweep removed and lets plurals, SQL aliases and words through."""
    source = tmp_path / "x.py"
    source.write_text('a = f"{n:.1f} s"\nb = f"{mb:,.0f} MB"\nc = f"{rate:.2f} per s"\nd = f"{count} systems"\n'
                      'e = f"SELECT 1 FROM {table} s WHERE 1"\ng = f"{format_duration_seconds(x)} elapsed"\nh = f"{x}s"\ni = f"{x:.2f}ms"\n', encoding="utf-8")
    assert [(line, unit) for line, unit, _ in _hand_formatted_in_python(source)] == [(1, "s"), (2, "MB"), (3, "per s"), (8, "ms")]


def test_every_allowed_entry_is_still_needed():
    """An entry whose file no longer formats that unit by hand is stale and goes."""
    live = set()
    for path in _python_files():
        where = str(path.relative_to(ROOT)).replace("\\", "/")
        live.update((where, unit) for _, unit, _ in _hand_formatted_in_python(path))
    for path in _templates():
        where = str(path.relative_to(ROOT)).replace("\\", "/")
        live.update((where, unit) for _, unit, _ in _hand_formatted_in_template(path))
    assert ALLOWED <= live, f"stale ALLOWED entries: {sorted(ALLOWED - live)}"
