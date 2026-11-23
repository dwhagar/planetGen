"""OPS.301: an update whose pull changed update.sh runs the new script from its start, under the same process id, and once only."""

import os
import subprocess

REPO_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _run(tmp_path, change_script):
    """A stand-in update.sh that sources the real shared functions, notes who it is, 'pulls' (rewrites itself when `change_script`) and asks to be run again."""
    script = tmp_path / "update.sh"
    body = f'''#!/bin/bash
SCRIPT_DIR="{tmp_path}"
source "{REPO_DIR}/scripts/deploy-common.sh"
echo "start pid=$$ reexec=${{PLANETGEN_UPDATE_REEXEC:-0}} held=$UPDATE_GUARD_HELD" >> "{tmp_path}/trace"
UPDATE_GUARD_HELD=1
before="$(script_sum "$SCRIPT_DIR/update.sh")"
if [[ "${{PLANETGEN_UPDATE_REEXEC:-0}}" != 1 && "{int(change_script)}" == 1 ]]; then
    cat "$SCRIPT_DIR/update.sh" > "$SCRIPT_DIR/update.sh.next"; echo "# new" >> "$SCRIPT_DIR/update.sh.next"
    mv "$SCRIPT_DIR/update.sh.next" "$SCRIPT_DIR/update.sh"
fi
reexec_when_script_changed "$SCRIPT_DIR/update.sh" "$before" "$@"
echo "end pid=$$" >> "{tmp_path}/trace"
'''
    script.write_text(body)
    done = subprocess.run(["bash", str(script), "--yes"], capture_output=True, text=True, timeout=60)
    assert done.returncode == 0, done.stderr
    return (tmp_path / "trace").read_text().splitlines()


def test_a_pull_that_changed_the_script_runs_the_new_one_with_the_same_pid_and_guard(tmp_path):
    lines = _run(tmp_path, True)
    first, second, end = lines
    assert first.startswith("start") and "reexec=0 held=0" in first
    assert "reexec=1 held=1" in second          # the guard the first run took is told to the second
    pids = {line.split("pid=")[1].split()[0] for line in lines}
    assert len(pids) == 1, lines                # exec keeps the process id, so the update's lock still names it



def test_a_pull_that_left_the_script_alone_does_not_run_it_again(tmp_path):
    lines = _run(tmp_path, False)
    assert len(lines) == 2 and lines[1].startswith("end")
