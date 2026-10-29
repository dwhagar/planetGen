"""
The shell's side of the logging standard (UX.123, `scripts/logrecord.sh`): update.sh's lines are the display lines
`logrecord.display_line` draws, so the Generate page, the command line and update.sh read alike.
"""

import os
import re
import subprocess

import pytest

from planetgen.util import logrecord

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
LIBRARY = os.path.join(REPO, "scripts", "logrecord.sh")
TIME = r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} \w+"


def _run(script, env=None):
    done = subprocess.run(["bash", "-c", f'source "{LIBRARY}"; {script}'], capture_output=True, text=True,
                          env={**os.environ, "TZ": "America/New_York", **(env or {})}, check=False)
    return done


def _without_time(line):
    match = re.match(rf"^({TIME}) (.*)$", line)
    assert match, line
    return match.group(2)


@pytest.mark.parametrize("severity, stream", [("info", "out"), ("notice", "out"), ("warning", "err"), ("error", "err")])
def test_a_shell_line_is_the_display_line_of_its_record(severity, stream):
    done = _run(f'log_{severity} "Placed the sectors"')
    line = (done.stdout if stream == "out" else done.stderr).rstrip("\n")
    other = done.stderr if stream == "out" else done.stdout
    assert other == "" and done.returncode == 0
    record = logrecord.Record.now(severity, "update", "Placed the sectors")
    assert _without_time(line) == _without_time(logrecord.display_line(record))
    assert re.match(rf"^{TIME} {severity.upper():<8} update: Placed the sectors$", line)


def test_a_stage_is_the_divider_of_the_stage_record():
    done = _run('log_stage 4 9 "Checking the generator\'s math"')
    record = logrecord.stage_record(logrecord.STAGE_START, "update", 4, 9, "Checking the generator's math")
    assert _without_time(done.stdout.rstrip("\n")) == _without_time(logrecord.display_line(record))
    assert "=== Stage 4 of 9: Checking the generator's math ===" in done.stdout


def test_a_message_of_several_lines_is_one_record_to_a_line_and_the_component_can_change():
    done = _run('LOG_COMPONENT=install log_warning "first\nsecond"')
    lines = done.stderr.splitlines()
    assert len(lines) == 2
    assert lines[0].endswith("WARNING  install: first") and lines[1].endswith("WARNING  install: second")


def test_update_sh_prints_its_own_lines_through_the_layer():
    text = open(os.path.join(REPO, "update.sh"), encoding="utf-8").read()
    assert 'source "$SCRIPT_DIR/scripts/logrecord.sh"' in text
    own = [line for line in text.splitlines() if re.match(r'\s*echo "', line)]
    assert own == [], f"update.sh's own lines go through scripts/logrecord.sh: {own}"
    assert all(f"log_stage {n} 9 " in text for n in range(1, 10))
