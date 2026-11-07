# tests/test_update_button_units.py

"""
SEC.43, the root side of the admin pages' Update button: the web user only drops `<state folder>/update.request`; the root-owned
`planetgen-update.path` starts `planetgen-update.service`, which removes the request and runs `update.sh`.

These tests read the unit files and run `install_maintenance_units` (scripts/deploy-common.sh, OPS.110) against a throwaway unit
folder with a fake `systemctl`; no systemd is needed.
"""

import os
import re
import stat
import subprocess

import pytest

REPO_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MAINT = os.path.join(REPO_DIR, "examples", "maintenance")
REQUEST = "/var/lib/planetGen/state/update.request"


def _units(name):
    """The unit file's directives (comments and blank lines dropped), as a list of lines."""
    with open(os.path.join(MAINT, name), encoding="utf-8") as f:
        return [line.strip() for line in f if line.strip() and not line.lstrip().startswith("#")]


def test_the_path_unit_watches_only_the_fixed_request_file_and_starts_the_update_service():
    lines = _units("planetgen-update.path")
    assert f"PathExists={REQUEST}" in lines
    assert "Unit=planetgen-update.service" in lines
    assert "WantedBy=multi-user.target" in lines
    # a fixed name: no wildcard, no specifier, nothing read from the file
    assert not any(line.startswith(("PathExistsGlob", "PathChanged", "PathModified", "DirectoryNotEmpty")) for line in lines)


def test_the_service_removes_the_request_before_update_sh_and_takes_nothing_from_it():
    lines = _units("planetgen-update.service")
    exec_start = [line for line in lines if line.startswith("ExecStart=")]
    pre = [line for line in lines if line.startswith("ExecStartPre=")]
    assert exec_start == ["ExecStart=/var/lib/planetGen/update.sh"], "no argument from a request reaches root"
    assert pre == [f"ExecStartPre=/bin/rm -f -- {REQUEST}"], "a failed removal must stop the service, not be tolerated with '-'"
    for line in lines:
        if "=" in line:
            assert not re.search(r"[$%]", line.split("=", 1)[1]), f"no variable or specifier in {line}"
    assert lines.index(pre[0]) < lines.index(exec_start[0])
    assert "Type=oneshot" in lines


def test_the_service_has_no_input_the_web_user_can_write():
    text = " ".join(_units("planetgen-update.service"))
    assert "EnvironmentFile=-/etc/planetgen/maintenance.env" in text  # root's, mode 600
    assert "state" not in text.replace(REQUEST, "")


@pytest.fixture
def deployed(tmp_path):
    """A fake checkout with the unit files, a fake systemctl that logs its calls, and an empty unit folder."""
    checkout = tmp_path / "check out"
    checkout.mkdir()
    (checkout / "examples").symlink_to(os.path.join(REPO_DIR, "examples"))
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    log = tmp_path / "systemctl.log"
    fake = bin_dir / "systemctl"
    fake.write_text(f'#!/bin/bash\necho "$@" >> "{log}"\n[[ "$1" == is-enabled ]] && exit 1\nexit 0\n')
    fake.chmod(fake.stat().st_mode | stat.S_IXUSR)
    units = tmp_path / "units"
    units.mkdir()
    return checkout, bin_dir, log, units


def _run(deployed, state, checkout_name=None):
    checkout, bin_dir, log, units = deployed
    script_dir = checkout if checkout_name is None else checkout.parent / checkout_name
    script_dir.mkdir(exist_ok=True)
    if checkout_name is not None and not (script_dir / "examples").exists():
        (script_dir / "examples").symlink_to(os.path.join(REPO_DIR, "examples"))
    script = (f'SCRIPT_DIR="{script_dir}"; PYTHON=/nonexistent; source "{REPO_DIR}/scripts/deploy-common.sh"; '
              f'is_macos() {{ return 1; }}; UNIT_STATE_DIR="{state}" install_maintenance_units')
    env = dict(os.environ, PATH=f"{bin_dir}:{os.environ['PATH']}", SYSTEMD_UNIT_DIR=str(units))
    return subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=env, timeout=60)


def test_a_path_with_a_space_is_refused_not_written(deployed):
    checkout, _, log, units = deployed  # the checkout's folder name has a space
    result = _run(deployed, "/srv/state")
    assert result.returncode == 1                       # refused, and it says why
    assert "special character" in result.stdout + result.stderr
    assert not list(units.iterdir())
    assert not log.exists()


def test_the_units_are_written_with_the_real_paths_and_the_path_unit_is_enabled(deployed):
    checkout, _, log, units = deployed
    result = _run(deployed, "/srv/planetgen/state", checkout_name="deployed")
    assert result.returncode == 0, result.stderr
    path_unit = (units / "planetgen-update.path").read_text()
    service = (units / "planetgen-update.service").read_text()
    root = str(checkout.parent / "deployed")
    assert "PathExists=/srv/planetgen/state/update.request" in path_unit
    assert "ExecStartPre=/bin/rm -f -- /srv/planetgen/state/update.request" in service
    assert f"ExecStart={root}/update.sh" in service
    assert f"WorkingDirectory={root}" in service
    assert "/var/lib/planetGen" not in path_unit.replace("(`PLANETGEN_STATE_DIR`", "")
    calls = log.read_text().splitlines()
    assert "daemon-reload" in calls
    assert "enable --now planetgen-update.path" in calls
    assert oct(stat.S_IMODE((units / "planetgen-update.path").stat().st_mode)) == "0o644"


def test_a_second_run_changes_nothing_and_does_not_reload(deployed):
    checkout, _, log, units = deployed
    _run(deployed, "/srv/planetgen/state", checkout_name="deployed")
    before = {p.name: p.read_text() for p in units.iterdir()}
    log.unlink()
    _run(deployed, "/srv/planetgen/state", checkout_name="deployed")
    assert {p.name: p.read_text() for p in units.iterdir()} == before
    assert not log.exists()                      # nothing to reload or enable


def test_the_deploy_scripts_set_the_units_up_and_the_maintenance_installer_has_the_opt_out():
    for name in ("install.sh", "update.sh"):
        with open(os.path.join(REPO_DIR, name), encoding="utf-8") as f:
            assert re.search(r"^install_maintenance_units( \|\|.*)?$", f.read(), re.M), name
    with open(os.path.join(MAINT, "install-maintenance-timer.sh"), encoding="utf-8") as f:
        text = f.read()
    assert "--skip-update-button" in text and "planetgen-update.path" in text
