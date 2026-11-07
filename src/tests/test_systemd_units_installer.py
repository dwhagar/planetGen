# tests/test_systemd_units_installer.py

"""
OPS.110: `install.sh` and `update.sh` render every maintenance systemd unit for the checkout and virtual environment they run from
(`install_maintenance_units`, scripts/deploy-common.sh), and `examples/maintenance/install-maintenance-timer.sh` is a thin caller of
the same function. The tests run it against a throwaway unit folder with a fake `systemctl`; no systemd and no root are needed.
"""

import os
import stat
import subprocess

import pytest

REPO_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
INSTALLER = os.path.join(REPO_DIR, "examples", "maintenance", "install-maintenance-timer.sh")
DEFAULTS = ("/var/lib/planetGen", "/usr/local/planetgen/venv")
UNITS = [
    "planetgen-orbits@.service", "planetgen-orbits@.timer", "planetgen-update.service", "planetgen-update.timer",
    "planetgen-update.path", "planetgen-tiles@.service", "planetgen-tiles@.timer",
]
SMOKE = ["planetgen-smoke.service", "planetgen-smoke.timer"]


@pytest.fixture
def host(tmp_path):
    """A checkout at a path of its own, a fake systemctl that logs its calls, and an empty unit folder."""
    checkout = tmp_path / "opt" / "pg"
    checkout.mkdir(parents=True)
    for name in ("examples", "scripts"):
        (checkout / name).symlink_to(os.path.join(REPO_DIR, name))
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    log = tmp_path / "systemctl.log"
    fake = bin_dir / "systemctl"
    enabled = tmp_path / "enabled.txt"
    fake.write_text(f'''#!/bin/bash
echo "$@" >> "{log}"
case "$1" in
    enable) shift; [[ "$1" == --now ]] && shift; echo "$1" >> "{enabled}" ;;
    is-enabled) shift; [[ "$1" == --quiet ]] && shift; grep -qx "$1" "{enabled}" 2>/dev/null || exit 1 ;;
esac
exit 0
''')
    fake.chmod(fake.stat().st_mode | stat.S_IXUSR)
    units = tmp_path / "units"
    units.mkdir()
    venv = tmp_path / "venv"
    return {"checkout": checkout, "bin": bin_dir, "log": log, "units": units, "venv": venv, "tmp": tmp_path, "state": tmp_path / "state"}


def _function(host, arguments="", checkout=None):
    """Runs `install_maintenance_units` from deploy-common.sh as a deploy would."""
    checkout = checkout or host["checkout"]
    script = (f'SCRIPT_DIR="{checkout}"; PYTHON=/nonexistent; source "{REPO_DIR}/scripts/deploy-common.sh"; '
              f'is_macos() {{ return 1; }}; install_maintenance_units {arguments}')
    return subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=_env(host), timeout=60)


def _env(host):
    return dict(os.environ, PATH=f"{host['bin']}:{os.environ['PATH']}", SYSTEMD_UNIT_DIR=str(host["units"]),
                PLANETGEN_VENV=str(host["venv"]), UNIT_STATE_DIR=str(host["state"]),
                MAINTENANCE_ENV_FILE=str(host["tmp"] / "etc" / "maintenance.env"), PLANETGEN_MYSQL_DATABASE="")


def _calls(host):
    return host["log"].read_text().splitlines() if host["log"].exists() else []


def test_every_unit_is_rendered_for_the_checkout_and_the_venv_with_no_default_path_left(host):
    result = _function(host)
    assert result.returncode == 0, result.stdout + result.stderr
    assert sorted(path.name for path in host["units"].iterdir()) == sorted(UNITS)
    checkout, venv = str(host["checkout"]), str(host["venv"])
    for path in host["units"].iterdir():
        text = path.read_text()
        for default in DEFAULTS:
            assert default not in text, f"{path.name} still names {default}"
        assert stat.S_IMODE(path.stat().st_mode) == 0o644
    service = (host["units"] / "planetgen-orbits@.service").read_text()
    assert f"WorkingDirectory={checkout}/src" in service
    assert f"ExecStart={venv}/bin/python -m planetgen.cli.orbits" in service
    assert f"ExecStart={checkout}/update.sh" in (host["units"] / "planetgen-update.service").read_text()
    assert f"PathExists={host['state']}/update.request" in (host["units"] / "planetgen-update.path").read_text()
    assert f"ExecStartPost=-{venv}/bin/python -m planetgen.cli.generate warm --progress" in (host["units"] / "planetgen-tiles@.service").read_text()


def test_the_smoke_units_render_with_their_option_and_the_output_path_follows_the_checkout(host):
    assert _function(host, "--with-smoke-timer").returncode == 0
    for name in SMOKE:
        text = (host["units"] / name).read_text()
        for default in DEFAULTS:
            assert default not in text, f"{name} still names {default}"
    assert f"--output {host['checkout']}/smoke-last.json" in (host["units"] / "planetgen-smoke.service").read_text()


def test_a_deploy_enables_what_it_installs_but_never_the_unattended_update_and_says_how(host):
    result = _function(host)
    calls = _calls(host)
    for name in ("planetgen-update.path", "planetgen-orbits@planetgen.timer", "planetgen-tiles@planetgen.timer"):
        assert f"enable --now {name}" in calls, name
    assert "enable --now planetgen-update.timer" not in calls
    assert (host["units"] / "planetgen-update.timer").is_file()              # installed, off
    text = result.stdout + result.stderr
    assert "sudo systemctl enable --now planetgen-update.timer" in text and "no review" in text
    assert "daemon-reload" in calls


def test_a_second_deploy_changes_nothing_and_says_again_while_the_update_timer_is_off(host):
    _function(host)
    before = {path.name: path.read_text() for path in host["units"].iterdir()}
    host["log"].unlink()
    result = _function(host)
    assert result.returncode == 0
    assert {path.name: path.read_text() for path in host["units"].iterdir()} == before
    assert not any(call.startswith(("enable", "daemon-reload")) for call in _calls(host))    # nothing reloaded, nothing switched on
    assert "sudo systemctl enable --now planetgen-update.timer" in result.stdout + result.stderr


def test_once_the_operator_turns_the_update_timer_on_the_deploy_stops_saying_so(host):
    _function(host)
    subprocess.run([str(host["bin"] / "systemctl"), "enable", "--now", "planetgen-update.timer"], check=True)
    result = _function(host)
    assert "sudo systemctl enable --now planetgen-update.timer" not in result.stdout + result.stderr


def test_the_installer_script_enables_the_update_timer_when_asked(host):
    result = _function(host, "--enable-all")
    assert "enable --now planetgen-update.timer" in _calls(host)
    assert "To turn it on" not in result.stdout + result.stderr


def test_a_unit_that_differs_is_replaced_and_the_replacement_is_logged(host):
    _function(host)
    stale = host["units"] / "planetgen-orbits@.service"
    stale.write_text(stale.read_text().replace(str(host["checkout"]), "/var/lib/planetGen"))
    host["log"].unlink()
    result = _function(host)
    assert result.returncode == 0
    assert "replaced" in result.stdout + result.stderr and "planetgen-orbits@.service" in result.stdout + result.stderr
    assert str(host["checkout"]) in stale.read_text() and "/var/lib/planetGen" not in stale.read_text()
    calls = _calls(host)
    assert "daemon-reload" in calls
    assert not any(call.startswith("enable") for call in calls)     # a replaced unit keeps the state the host gave it


def test_an_installed_optional_unit_is_kept_right_even_when_not_asked_for(host):
    _function(host, "--with-smoke-timer")
    smoke = host["units"] / "planetgen-smoke.service"
    smoke.write_text(smoke.read_text().replace(str(host["checkout"]), "/elsewhere"))
    assert _function(host).returncode == 0              # a later deploy: no option, the unit is still there
    assert "/elsewhere" not in smoke.read_text() and str(host["checkout"]) in smoke.read_text()


def test_a_path_with_a_space_or_special_character_is_refused_and_nothing_is_written(host):
    spaced = host["tmp"] / "my checkout"
    spaced.mkdir()
    for name in ("examples", "scripts"):
        (spaced / name).symlink_to(os.path.join(REPO_DIR, name))
    result = _function(host, checkout=spaced)
    assert result.returncode == 1
    assert "special character" in result.stdout + result.stderr
    assert not list(host["units"].iterdir()) and _calls(host) == []
    host["venv"] = host["tmp"] / "ven$v"                                    # the venv path is checked the same way
    result = _function(host)
    assert result.returncode == 1 and not list(host["units"].iterdir())


def test_the_options_leave_units_out_and_name_the_databases(host):
    result = _function(host, "--skip-update-timer --skip-update-button --skip-tiles-timer --enable-all alpha beta")
    assert result.returncode == 0, result.stderr
    assert sorted(path.name for path in host["units"].iterdir()) == ["planetgen-orbits@.service", "planetgen-orbits@.timer"]
    calls = _calls(host)
    assert "enable --now planetgen-orbits@alpha.timer" in calls and "enable --now planetgen-orbits@beta.timer" in calls
    assert not any("update" in call or "tiles" in call for call in calls if call.startswith("enable"))
    assert "To turn it on" not in result.stdout + result.stderr          # the update timer was left out on purpose: nothing to offer


def test_the_gunicorn_unit_is_never_installed(host):
    _function(host, "--enable-all --with-smoke-timer")
    assert not (host["units"] / "planetgen-gunicorn.service").exists()
    for path in (INSTALLER, os.path.join(REPO_DIR, "scripts", "deploy-common.sh")):
        with open(path, encoding="utf-8") as f:
            assert "render_systemd_unit \"$REPO_DIR/examples/systemd" not in f.read()


def test_the_installer_script_is_a_thin_caller_of_the_same_function(host):
    env = dict(_env(host), PATH=f"{host['bin']}:{os.environ['PATH']}")
    result = subprocess.run(["bash", INSTALLER, "--with-smoke-timer", "orbitdb"], capture_output=True, text=True, env=env, timeout=60)
    if os.geteuid() != 0:
        assert result.returncode == 1 and "must be run as root" in result.stderr
        pytest.skip("the installer script insists on root")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "enable --now planetgen-orbits@orbitdb.timer" in _calls(host)
    assert (host["tmp"] / "etc" / "maintenance.env").is_file()
    assert stat.S_IMODE((host["tmp"] / "etc" / "maintenance.env").stat().st_mode) == 0o600
    for name in UNITS + SMOKE:
        assert (host["units"] / name).is_file(), name
    with open(INSTALLER, encoding="utf-8") as f:
        text = f.read()
    assert "install_maintenance_units" in text and "sed -e" not in text    # no copy of the rendering lives in the script


def test_the_deploy_scripts_call_the_function_on_every_run():
    for name in ("install.sh", "update.sh"):
        with open(os.path.join(REPO_DIR, name), encoding="utf-8") as f:
            assert "install_maintenance_units" in f.read(), name
