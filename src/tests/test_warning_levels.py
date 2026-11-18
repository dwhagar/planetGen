"""OPS.167: a warning is logged as a warning. A line that says "Warning:" at normal level never reaches a warning reader (error.log, the summary's warning list)."""

import re
from pathlib import Path

import planetgen

SOURCE = Path(planetgen.__file__).parent
SEVERITY_WORD_AT_NORMAL = re.compile(r"""log\.(?:normal|info)\(\s*f?["'](?:Warning|WARNING|warning)\b""")


def test_no_line_says_warning_at_normal_level():
    offenders = []
    for path in SOURCE.rglob("*.py"):
        for number, line in enumerate(path.read_text().splitlines(), 1):
            if SEVERITY_WORD_AT_NORMAL.search(line):
                offenders.append(f"{path.relative_to(SOURCE)}:{number}")
    assert not offenders, f"use log.warning (the severity says it; leave the word out of the text): {offenders}"

