# tests/test_command_records.py

"""
UX.122: the job-step commands (`datasets`, `reset`, `orbits`) say what they do through the logging layer, so a line on
the console is the plain line of a record, the same record `--record-file` keeps: nothing the console showed is missing
from the job's `records.log`, and a warning or an error goes to stderr as these commands' messages always did.
"""

import os
import subprocess
import sys

import pytest

from planetgen.util import logrecord
from tests import datasets_support as support
from tests.bughunt_support import mysql_argv
from tests.test_update_guard import begin

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def run_command(module, arguments, tmp_path, env=None):
    """Runs `python -m module ... --record-file F`; returns (result, the records F holds)."""
    record_file = tmp_path / "records.log"
    result = subprocess.run([sys.executable, "-m", module, *arguments, "--record-file", str(record_file)],
                            capture_output=True, text=True, cwd=REPO, env=env or os.environ.copy(), timeout=180)
    lines = record_file.read_text(encoding="utf-8").splitlines() if record_file.exists() else []
    return result, [logrecord.parse_file_line(line) for line in lines]


def assert_console_is_the_records(result, records):
    """Stdout holds the plain lines of the records below `warning`, stderr those of `warning` and above, in order."""
    quiet = [logrecord.plain_line(record) for record in records
             if logrecord.SYSLOG_CODES[record.severity] > logrecord.SYSLOG_CODES["warning"]]
    loud = [logrecord.plain_line(record) for record in records
            if logrecord.SYSLOG_CODES[record.severity] <= logrecord.SYSLOG_CODES["warning"]]
    assert records, "the command wrote no records"
    assert result.stdout.splitlines() == quiet
    assert result.stderr.splitlines() == loud


def test_the_datasets_command_writes_what_it_prints_as_records(tmp_path, monkeypatch):
    support.reset_caches(tmp_path / "cache", monkeypatch)
    env = dict(os.environ, PLANETGEN_DATA_DIR=str(tmp_path / "cache"))
    result, records = run_command("planetgen.cli.datasets", ["status"], tmp_path, env)
    assert result.returncode == 1
    assert "isochrones: missing (SC1)" in result.stdout
    assert_console_is_the_records(result, records)


def test_the_orbital_update_refused_during_an_update_is_an_error_record_on_stderr(tmp_path):
    assert begin(sources=[])
    result, records = run_command("planetgen.cli.orbits", ["--mysql-host", "no-such-host.invalid"], tmp_path)
    assert result.returncode == 75
    assert [record.severity for record in records] == ["error"]
    assert result.stderr.startswith("error: The orbital update did not run and moved nothing")
    assert_console_is_the_records(result, records)


def test_a_reset_dry_run_writes_what_it_prints_as_records(mysql_config, tmp_path):
    result, records = run_command("planetgen.cli.reset", [*mysql_argv(mysql_config), "--dry-run"], tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    assert any(record.message.startswith("Would truncate") for record in records)
    assert_console_is_the_records(result, records)


def test_a_reset_that_cannot_reach_the_database_is_an_error_record_on_stderr(tmp_path):
    result, records = run_command("planetgen.cli.reset", ["--mysql-host", "no-such-host.invalid", "--yes"], tmp_path)
    assert result.returncode == 1
    assert [record.severity for record in records] == ["error"]
    assert_console_is_the_records(result, records)
