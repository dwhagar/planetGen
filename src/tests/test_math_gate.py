# tests/test_math_gate.py

"""
The math check in front of bulk generation (TEST.68): `planetgen
check-math`, the gate every bulk command and the Sector page's
neighbourhood generation pass first (refusing, naming the failed checks,
and writing nothing), the Generate page's first job step, and the warning
`update.sh`/`update.ps1` give.
"""

import os
import subprocess
import sys

import pytest

from planetgen.cli import generate as generate_cli
from planetgen.db import store
from planetgen.generation import run_galaxy
from planetgen.admin import activity_log
from planetgen.physics import mathcheck as mathCheck

pytestmark = pytest.mark.mathcheck

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _failed():
    check = mathCheck.Check("snow_line_1_lsun", "reference", "formation.snow_line_au(L_sun)",
                            lambda: 2.6, 2.7, 1e-6, "test")
    return mathCheck.run_check(check)


@pytest.fixture
def broken_math(monkeypatch):
    """The math check has failed in this process."""
    monkeypatch.setattr(mathCheck, "_STARTUP_RESULTS", [_failed()])


@pytest.fixture
def handlers(monkeypatch):
    """Records which subcommand handlers ran and every activity log line,
    instead of generating anything."""
    ran, events = [], []
    for name in list(generate_cli._COMMAND_HANDLERS):
        monkeypatch.setitem(generate_cli._COMMAND_HANDLERS, name, lambda args, name=name: ran.append(name))
    monkeypatch.setattr(activity_log, "event", lambda *a, **k: events.append((a, k)))
    return ran, events


def _main(monkeypatch, *argv):
    monkeypatch.setattr(sys, "argv", ["planetgen", *argv])
    generate_cli.main()


def test_check_math_command_passes_on_real_math(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["planetgen", "check-math"])
    generate_cli.main()


def test_check_math_command_exits_1_naming_the_failure(monkeypatch):
    bad = mathCheck.Check("bad_check", "invariant", "bad", lambda: 1.0, 0.0, 0.0, "test", mode="max")
    monkeypatch.setattr(mathCheck, "all_checks", lambda: [bad])
    monkeypatch.setattr(sys, "argv", ["planetgen", "check-math"])
    with pytest.raises(SystemExit) as exit_info:
        generate_cli.main()
    assert exit_info.value.code == 1


def test_check_math_runs_from_the_command_line():
    result = subprocess.run([sys.executable, "-m", "planetgen.cli.generate", "check-math", "-v"],
                            capture_output=True, text=True, cwd=REPO, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "Math check passed" in result.stdout + result.stderr
    assert "sun_luminosity" in result.stdout + result.stderr


@pytest.mark.parametrize("argv", [
    ["galaxy", "--ring", "5"],
    ["plan"],
    ["population"],
    ["sector", "--num-sectors", "3"],
])
def test_bulk_runs_refuse_and_write_nothing(broken_math, handlers, monkeypatch, capsys, argv):
    ran, events = handlers
    with pytest.raises(SystemExit) as exit_info:
        _main(monkeypatch, *argv)
    assert exit_info.value.code == 1
    assert ran == [] and events == []
    out = capsys.readouterr()
    assert "snow_line_1_lsun" in out.out + out.err


@pytest.mark.parametrize("argv, command", [
    (["sector"], "sector"),
    (["system", "--output", "-"], "system"),
    (["phenomenon", "--type", "nebula"], "phenomenon"),
])
def test_single_object_runs_are_not_gated(broken_math, handlers, monkeypatch, argv, command):
    ran, _ = handlers
    _main(monkeypatch, *argv)
    assert ran == [command]


def test_bulk_runs_go_ahead_when_the_math_is_right(handlers, monkeypatch):
    monkeypatch.setattr(mathCheck, "_STARTUP_RESULTS", [])
    ran, _ = handlers
    _main(monkeypatch, "galaxy", "--ring", "5")
    assert ran == ["galaxy"]


def test_neighborhood_generation_refuses_before_touching_the_database(broken_math, monkeypatch):
    def no_database(*args, **kwargs):
        raise AssertionError("the database was touched")

    monkeypatch.setattr(store, "get_connection", no_database)
    with pytest.raises(run_galaxy.MathCheckFailed, match="snow_line_1_lsun"):
        run_galaxy.generate_sector_neighborhood(1)


def test_neighborhood_estimate_is_not_gated(broken_math, monkeypatch):
    """An estimate writes nothing, so it isn't refused (it gets as far as
    the database)."""
    class Reached(Exception):
        pass

    def reached(*args, **kwargs):
        raise Reached

    monkeypatch.setattr(store, "get_connection", reached)
    with pytest.raises(Reached):
        run_galaxy.generate_sector_neighborhood(1, estimate_only=True)


def test_the_neighborhood_route_reports_the_failure(broken_math, monkeypatch):
    """The Sector page's button gets a 409 naming the failed check."""
    from planetgen.api import routes
    from planetgen.web.app import create_app
    from planetgen.api.config import Config

    class _Config(Config):
        SECRET_KEY = "test-secret"

    monkeypatch.setattr(routes, "_resolve_requested_write_db_config", lambda: None)
    app = create_app(_Config)
    view = app.view_functions["api.generate_sector_neighborhood_route"]
    while hasattr(view, "__wrapped__"):  # past the admin and rate-limit wrappers
        view = view.__wrapped__
    with app.test_request_context("/api/sectors/1/generate-neighborhood", method="POST", json={}):
        with pytest.raises(routes.ApiError) as error:
            view(1)
    assert error.value.status_code == 409
    assert "snow_line_1_lsun" in str(error.value)


def _deploy_common_shell(python, script):
    """Runs `script` in bash with scripts/deploy-common.sh sourced and
    PYTHON set to `python`."""
    code = f'set -euo pipefail\nSCRIPT_DIR="{REPO}"\nPYTHON="{python}"\n' \
           f'source "$SCRIPT_DIR/scripts/deploy-common.sh"\n{script}\n'
    return subprocess.run(["bash", "-c", code], capture_output=True, text=True, timeout=120)


@pytest.mark.skipif(os.name == "nt", reason="bash script")
def test_update_sh_check_math_passes(tmp_path):
    result = _deploy_common_shell(sys.executable, 'check_math; echo "failed=$MATH_CHECK_FAILED"')
    assert result.returncode == 0, result.stderr
    assert "failed=0" in result.stdout


@pytest.mark.skipif(os.name == "nt", reason="bash script")
def test_update_sh_check_math_failure_only_warns(tmp_path):
    fake = tmp_path / "python"
    fake.write_text("#!/bin/sh\necho 'FAIL bad_check'\nexit 1\n")
    fake.chmod(0o755)
    result = _deploy_common_shell(
        str(fake), 'check_math; echo "failed=$MATH_CHECK_FAILED"; POPULATION=1 offer_population_pass')
    assert result.returncode == 0, result.stderr
    assert "failed=1" in result.stdout
    assert "warning: the math check failed" in result.stderr
    assert "Skipping the population pass: the math check failed." in result.stdout


def test_update_scripts_run_the_math_check():
    with open(os.path.join(REPO, "update.sh"), encoding="utf-8") as handle:
        update_sh = handle.read()
    assert "check_math" in update_sh and "MATH_CHECK_FAILED" in update_sh
    with open(os.path.join(REPO, "update.ps1"), encoding="utf-8") as handle:
        update_ps1 = handle.read()
    assert "Test-MathCheck" in update_ps1 and "$mathOk" in update_ps1
