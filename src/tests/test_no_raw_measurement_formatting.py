"""
No browser file formats a measurement by hand (UX.22, UX.81).

A measurement on a page goes through the unit ladders (`unitladder.js`) or the
number rule (`numberformat.js`), so it reads like the same value on every
other page. `toPrecision`, `toExponential` and `toLocaleString` are how a file
used to round, group and switch to scientific notation by itself; this test
fails when one is used outside the two formatters. `toFixed` stays free: it
sets a map coordinate's decimals, a pixel size, a form field's value or a
shader constant, none of which is a measurement the reader sees in a unit.
"""
import os
import re

import pytest

STATIC = os.path.join(os.path.dirname(__file__), "..", "html", "static")

FORMATTERS = {"numberformat.js"}
"""set: The files whose job it is to round and group numbers."""

GENERATED = {"unitladders.js"}

ALLOWED = {
    # A ring, wedge or slot number is an address, never scientific (its own comment says so).
    ("galaxystages.js", "toLocaleString"),
}
"""set: `(file, method)` pairs that are not measurements."""

RAW = re.compile(r"\.(toPrecision|toExponential|toLocaleString)\(")


def _scripts():
    for name in sorted(os.listdir(STATIC)):
        if name.endswith(".js") and name not in FORMATTERS | GENERATED:
            yield name


@pytest.mark.parametrize("name", list(_scripts()))
def test_a_script_does_not_round_or_group_a_measurement_itself(name):
    with open(os.path.join(STATIC, name), encoding="utf-8") as handle:
        text = handle.read()
    found = [(match.group(1), text.count("\n", 0, match.start()) + 1) for match in RAW.finditer(text)]
    unexpected = [(method, line) for method, line in found if (name, method) not in ALLOWED]
    assert not unexpected, (
        f"{name} formats a number by hand ({unexpected}): use unitladder.js for a measurement or "
        f"numberformat.js for a count (docs/design/units-and-number-formatting.md)")
