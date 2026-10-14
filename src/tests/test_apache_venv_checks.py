# tests/test_apache_venv_checks.py

"""
OPS.46: `deploy-common.sh`'s checks of Apache's vhost against mod_wsgi's rules for a virtual environment. They only
warn, with the line to add, and never edit the vhost. Each case runs the real functions in bash against a copy of
the vhost (`APACHE_CONF_DIRS`) and a venv folder.
"""

import os
import shutil
import subprocess
import tempfile

import pytest

ROOT = os.path.realpath(os.path.join(os.path.dirname(__file__), "..", ".."))
GOOD = """
<VirtualHost *:80>
    WSGIDaemonProcess planetgen-api processes=1 threads=5 request-timeout=60 python-home={venv}
    WSGIScriptAlias / {root}/src/html/wsgi.py process-group=planetgen-api application-group=%{{GLOBAL}}
</VirtualHost>
"""


def _run(tmp_path, vhost, function="check_apache_venv_rules", user=None, venv=None):
    sites = tmp_path / "sites"
    sites.mkdir(exist_ok=True)
    (sites / "planetgen.conf").write_text(vhost.format(venv=venv or "/usr/local/planetgen/venv", root=ROOT))
    envvars = tmp_path / "envvars"
    envvars.write_text(f"export APACHE_RUN_USER={user or 'nobody'}\nexport APACHE_RUN_GROUP=nogroup\n")
    env = dict(os.environ, SCRIPT_DIR=ROOT, APACHE_CONF_DIRS=str(sites), APACHE_MAIN_CONFS="",
               APACHE_ENVVARS=str(envvars), PLANETGEN_VENV=venv or "/usr/local/planetgen/venv")
    done = subprocess.run(["bash", "-c", f'source "{ROOT}/scripts/deploy-common.sh"; {function}'], env=env,
                          capture_output=True, text=True)
    return done.stdout, done.stderr


@pytest.fixture
def base():
    """A folder anyone can enter (pytest's own temporary folder is private to its owner)."""
    folder = tempfile.mkdtemp(prefix="pg-venv-check-", dir="/tmp")
    os.chmod(folder, 0o755)
    yield folder
    shutil.rmtree(folder, ignore_errors=True)


@pytest.fixture
def venv(base):
    """A readable venv-shaped folder, so the access check has something to pass."""
    folder = os.path.join(base, "planetgen", "venv")
    for part in ("bin", "lib"):
        os.makedirs(os.path.join(folder, part))
    with open(os.path.join(folder, "pyvenv.cfg"), "w") as f:
        f.write("home = /usr/bin\n")
    for path in (os.path.join(base, "planetgen"), folder, os.path.join(folder, "bin"), os.path.join(folder, "lib")):
        os.chmod(path, 0o755)
    os.chmod(os.path.join(folder, "pyvenv.cfg"), 0o644)
    return folder


def test_a_vhost_that_follows_the_rules_says_so(tmp_path, venv):
    out, err = _run(tmp_path, GOOD, venv=venv)
    assert f"Apache runs the site under {venv}" in out
    assert "warning" not in err


def test_python_home_must_be_the_venv_root_not_the_interpreter(tmp_path, venv):
    vhost = GOOD.replace("python-home={venv}", "python-home={venv}/bin/python")
    out, err = _run(tmp_path, vhost, venv=venv)
    assert "a path into the venv's bin/" in err
    assert f"python-home={venv}" in err and "Apache runs the site under" not in out


@pytest.mark.parametrize("home", ["/usr/local/planetgen", "/opt/other/venv", "/usr/local/planetgen/venv/lib"])
def test_python_home_pointing_anywhere_else_is_warned(tmp_path, venv, home):
    out, err = _run(tmp_path, GOOD.replace("python-home={venv}", f"python-home={home}"), venv=venv)
    assert f"python-home={home}, not planetGen's virtual environment" in err
    assert "Apache runs the site under" not in out


def test_a_trailing_slash_on_the_venv_root_is_fine(tmp_path, venv):
    out, err = _run(tmp_path, GOOD.replace("python-home={venv}", "python-home={venv}/"), venv=venv)
    assert "Apache runs the site under" in out and "warning" not in err


def test_a_missing_python_home_is_warned_with_the_line_to_add(tmp_path, venv):
    _out, err = _run(tmp_path, GOOD.replace(" python-home={venv}", ""), venv=venv)
    assert "does not run the site under planetGen's virtual environment" in err
    assert f"Add python-home={venv}" in err


def test_embedded_mode_is_warned(tmp_path, venv):
    vhost = """
<VirtualHost *:80>
    WSGIScriptAlias / {root}/src/html/wsgi.py
</VirtualHost>
"""
    _out, err = _run(tmp_path, vhost, venv=venv)
    assert "embedded mode" in err and "process-group=planetgen-api application-group=%{GLOBAL}" in err


def test_a_process_group_directive_counts_as_the_option(tmp_path, venv):
    vhost = """
WSGIDaemonProcess other processes=1 threads=5 python-home={venv}
<VirtualHost *:80>
    WSGIProcessGroup other
    WSGIApplicationGroup %{{GLOBAL}}
    WSGIScriptAlias / {root}/src/html/wsgi.py
</VirtualHost>
"""
    out, err = _run(tmp_path, vhost, venv=venv)
    assert "process group other" in out and "warning" not in err


def test_a_group_with_no_daemon_process_is_warned(tmp_path, venv):
    vhost = GOOD.replace("    WSGIDaemonProcess planetgen-api processes=1 threads=5 request-timeout=60 python-home={venv}\n", "")
    _out, err = _run(tmp_path, vhost, venv=venv)
    assert "has no WSGIDaemonProcess planetgen-api line" in err


def test_the_main_interpreter_is_required(tmp_path, venv):
    _out, err = _run(tmp_path, GOOD.replace(" application-group=%{{GLOBAL}}", ""), venv=venv)
    assert "main Python interpreter" in err and "application-group=%{GLOBAL}" in err


def test_commented_lines_do_not_count(tmp_path, venv):
    vhost = "# " + GOOD.strip().replace("\n", "\n# ") + "\n"
    out, err = _run(tmp_path, vhost, venv=venv)
    assert out == "" and err == ""


def test_no_vhost_for_planetgen_says_nothing(tmp_path, venv):
    assert _run(tmp_path, "<VirtualHost *:80>\n</VirtualHost>\n", venv=venv) == ("", "")


needs_root = pytest.mark.skipif(os.geteuid() != 0 or not (shutil.which("runuser") or shutil.which("sudo")),
                                reason="runs a command as Apache's user: needs root and runuser or sudo")


@needs_root
def test_apache_user_who_can_read_the_venv_is_confirmed(tmp_path, venv):
    out, err = _run(tmp_path, GOOD, venv=venv)
    assert "Apache's user nobody can read" in out and "warning" not in err


@needs_root
def test_a_directory_apache_cannot_enter_is_named(tmp_path, base, venv):
    blocked = os.path.join(base, "planetgen")
    os.chmod(blocked, 0o700)
    _out, err = _run(tmp_path, GOOD, venv=venv)
    assert f"Apache's user nobody can't read {blocked}" in err
    assert f"sudo chmod o+rx {blocked}" in err


@needs_root
def test_an_unreadable_pyvenv_cfg_is_named(tmp_path, venv):
    os.chmod(os.path.join(venv, "pyvenv.cfg"), 0o600)
    _out, err = _run(tmp_path, GOOD, venv=venv)
    assert f"can't read {venv}/pyvenv.cfg" in err


def test_the_checks_never_edit_the_vhost(tmp_path, venv):
    _run(tmp_path, GOOD.replace(" python-home={venv}", ""), venv=venv)
    assert "python-home" not in (tmp_path / "sites" / "planetgen.conf").read_text()
