# tests/test_math_gate.py

"""
The math check's command and its place in an update (TEST.68, PERF.170): `planetgen
check-math`, the warning `update.sh` gives, and that generation does not run it (Boss,
2026-10-23: the check runs at every update and from the admin page, not before every
generation).
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


@pytest.mark.parametrize("argv, command", [
    (["galaxy", "--shell", "5"], "galaxy"),
    (["plan"], "plan"),
    (["population"], "population"),
    (["sector", "--num-sectors", "3"], "sector"),
    (["sector"], "sector"),
    (["system", "--output", "-"], "system"),
    (["phenomenon", "--type", "nebula"], "phenomenon"),
])
def test_no_run_checks_the_math_first(broken_math, handlers, monkeypatch, argv, command):
    """A failed math check does not stop any run: the check is not part of generating."""
    ran, _ = handlers
    monkeypatch.setattr(mathCheck, "run_all", lambda *a, **k: pytest.fail("the math check ran"))
    monkeypatch.setattr(run_galaxy, "require_datasets", lambda *a, **k: [])
    _main(monkeypatch, *argv)
    assert ran == [command]


def test_the_neighborhood_run_does_not_check_the_math(broken_math, monkeypatch):
    """It gets as far as the database."""
    class Reached(Exception):
        pass

    def reached(*args, **kwargs):
        raise Reached

    monkeypatch.setattr(run_galaxy, "require_datasets", lambda *a, **k: [])
    monkeypatch.setattr(store, "get_connection", reached)
    with pytest.raises(Reached):
        run_galaxy.generate_sector_neighborhood(1)


def test_the_neighborhood_route_does_not_check_the_math(broken_math, monkeypatch):
    """The Sector page's button is not refused for a failed math check."""
    from planetgen.api import routes
    from planetgen.web.app import create_app
    from planetgen.api.config import Config

    class _Config(Config):
        SECRET_KEY = "test-secret"

    class Reached(Exception):
        pass

    def reached(*args, **kwargs):
        raise Reached

    monkeypatch.setattr(routes, "_resolve_requested_write_db_config", lambda: None)
    monkeypatch.setattr(run_galaxy, "require_datasets", lambda *a, **k: [])
    monkeypatch.setattr(run_galaxy, "generate_sector_neighborhood", reached)
    app = create_app(_Config)
    view = app.view_functions["api.generate_sector_neighborhood_route"]
    while hasattr(view, "__wrapped__"):  # past the admin and rate-limit wrappers
        view = view.__wrapped__
    with app.test_request_context("/api/sectors/1/generate-neighborhood", method="POST", json={}):
        with pytest.raises(Reached):
            view(1)


def _deploy_common_shell(python, script):
    """Runs `script` in bash with scripts/deploy-common.sh sourced and
    PYTHON set to `python`."""
    code = f'set -euo pipefail\nSCRIPT_DIR="{REPO}"\nPYTHON="{python}"\n' \
           f'source "$SCRIPT_DIR/scripts/deploy-common.sh"\n{script}\n'
    return subprocess.run(["bash", "-c", code], capture_output=True, text=True, timeout=120)


def test_update_sh_check_math_passes(tmp_path):
    result = _deploy_common_shell(sys.executable, 'check_math; echo "failed=$MATH_CHECK_FAILED"')
    assert result.returncode == 0, result.stderr
    assert "failed=0" in result.stdout


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


def test_update_script_runs_the_math_check():
    with open(os.path.join(REPO, "update.sh"), encoding="utf-8") as handle:
        update_sh = handle.read()
    assert "check_math" in update_sh and "MATH_CHECK_FAILED" in update_sh

