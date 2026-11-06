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


DEPLOY_SCRIPTS = ("update.sh", "install.sh", "scripts/deploy-common.sh", "scripts/install-python-deps.sh", "scripts/enable-swap.sh",
                  "examples/apache/setup-debug-log.sh", "examples/apache/set-permissions.sh", "examples/apache/create-cache-dir.sh")
"""tuple: The scripts whose lines are log records (OPS.98): the update, the install and the pieces they run."""


def test_every_deploy_script_prints_its_lines_through_the_layer():
    """OPS.98: a line a script prints carries a time, so no `echo "text"` is left but the ones that return a value to a caller."""
    for name in DEPLOY_SCRIPTS:
        text = open(os.path.join(REPO, name), encoding="utf-8").read()
        own = [line.strip() for line in text.splitlines() if re.match(r'\s*echo "', line) and "|" not in line and ";;" not in line]
        returns = [line for line in own if line.startswith(('echo "$candidate"', 'echo "${text# }', 'echo "$file"', 'echo "$2"', 'echo "$1"'))]
        assert [line for line in own if line not in returns] == [], f"{name}: lines that print without a time: {own}"
        assert "logrecord.sh" in text, name


def test_log_lines_makes_a_record_of_each_line_of_a_heredoc_and_skips_the_empty_ones():
    done = _run('log_lines notice <<EOF\nfirst\n\n  second\nEOF')
    lines = done.stdout.splitlines()
    assert len(lines) == 2 and lines[0].endswith("NOTICE   update: first") and lines[1].endswith("NOTICE   update:   second")


def test_a_duration_is_minutes_and_seconds_and_hours_from_an_hour_on():
    assert _run("log_duration 247").stdout == "4:07" and _run("log_duration 3847").stdout == "1:04:07"


def test_the_summary_ends_with_one_line_and_repeats_every_warning_but_not_the_lines_under_it():
    """OPS.100: the headline says how long it took and what it brought; every warning and error is repeated after it, headline only."""
    done = _run('log_summary_begin; log_info "fine"; log_warning "the math check failed\n  run it again"; log_error "no Redis" 2>/dev/null; '
                'log_summary Update "code 8.0.1"; echo "file:$LOG_WARNING_FILE"')
    out, err = done.stdout.splitlines(), done.stderr.splitlines()
    assert any(re.search(r"=== Update done after \d+:\d\d: code 8\.0\.1 ===$", line) for line in out)
    repeated = [line for line in err if "repeated here" in line]
    assert len(repeated) == 1 and "2 warnings during this update" in repeated[0]
    assert err[-2].endswith("WARNING  update:   the math check failed") and err[-1].endswith("WARNING  update:   no Redis")
    assert not any("run it again" in line for line in err[-2:])
    assert out[-1] == "file:", "the warnings file is forgotten"


def test_a_summary_without_warnings_says_so():
    done = _run('log_summary_begin; log_summary Install ""')
    assert done.stdout.splitlines()[-1].endswith("INFO     update: No warnings during this install.")


def test_the_shell_vocabulary_is_the_modules_and_the_committed_block_is_current():
    """OPS.104: the severities, the ones that go to stderr, the stage word and the divider are generated from `logrecord`."""
    text = open(LIBRARY, encoding="utf-8").read()
    assert logrecord.shell_vocabulary() in text
    done = subprocess.run(["python3", os.path.join(REPO, "scripts", "build_log_vocabulary.py"), "--check"], capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    shown = _run('echo "$LOG_SEVERITIES|$LOG_STDERR_SEVERITIES|$LOG_STAGE_WORD|$LOG_DIVIDER"').stdout.strip()
    assert shown == "|".join((" ".join(logrecord.SEVERITIES), " ".join(logrecord.STDERR_SEVERITIES), "Stage", logrecord.DIVIDER))


def test_a_severity_the_module_does_not_have_is_an_info_record():
    done = _run('log_record loud "Placed the sectors"')
    assert done.stdout.rstrip("\n").endswith("INFO     update: Placed the sectors") and done.stderr == ""


def test_a_script_log_keeps_the_records_from_before_its_folder_was_known_and_then_appends(tmp_path):
    """OPS.101: update.sh keeps every record from the start, writes them to update.log once the folder is known, and appends from then on."""
    target = tmp_path / "logs" / "update.log"
    done = _run(f'log_file_begin; log_info "before the folder"; log_file_to "{target}"; log_warning "after it"')
    lines = target.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 2 and lines[0].endswith("INFO     update: before the folder") and lines[1].endswith("WARNING  update: after it")
    assert done.stdout.count("before the folder") == 1 and "after it" in done.stderr, "the screen has them too"
    again = _run(f'log_file_begin; log_info "second run"; log_file_to "{target}"')
    assert again.returncode == 0 and len(target.read_text(encoding="utf-8").splitlines()) == 3, "a second run appends"


def test_a_script_log_that_cannot_be_written_leaves_the_screen_as_it_was(tmp_path):
    blocker = tmp_path / "file"
    blocker.write_text("")
    done = _run(f'log_file_begin; log_info "kept on the screen"; log_file_to "{blocker}/update.log"; log_info "still"')
    assert done.returncode == 0 and "kept on the screen" in done.stdout and "still" in done.stdout


def test_update_sh_starts_update_log_once_the_python_environment_is_known():
    text = open(os.path.join(REPO, "update.sh"), encoding="utf-8").read()
    assert "log_file_begin" in text and text.index("log_file_begin") < text.index("start_update_log") < text.index("ensure_update_guard required")
