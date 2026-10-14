"""
scripts/install-python-deps.sh keeps its own list of runtime requirements
to check the venv it makes (OPS.45). This checks that list matches
setup.py's install_requires plus its 'api' extra, so a dependency added or
bumped in setup.py can't be silently missed, and that every part of the
deployment runs from that venv.

Run with: pytest src/tests/test_install_python_deps.py
"""
import ast
import os
import re
import subprocess
import sys

import pytest

ROOT = os.path.join(os.path.dirname(__file__), "..", "..")


def _setup_requirements():
    with open(os.path.join(ROOT, "setup.py"), encoding="utf-8") as f:
        tree = ast.parse(f.read())
    call = next(
        node for node in ast.walk(tree)
        if isinstance(node, ast.Call) and getattr(node.func, "id", None) == "setup"
    )
    kwargs = {kw.arg: kw.value for kw in call.keywords}
    reqs = ast.literal_eval(kwargs["install_requires"])
    reqs += ast.literal_eval(kwargs["extras_require"])["api"]
    return {r.lower() for r in reqs}


def _script_requirements():
    with open(os.path.join(ROOT, "scripts", "install-python-deps.sh"), encoding="utf-8") as f:
        text = f.read()
    block = re.search(r"^REQUIREMENTS=\((.*?)^\)", text, re.S | re.M).group(1)
    return re.findall(r'"(\S+)"', block)


def test_script_requirements_match_setup_py():
    assert {r.lower() for r in _script_requirements()} == _setup_requirements()


def test_every_requirement_has_a_floor():
    for spec in _script_requirements():
        assert re.fullmatch(r"[a-z0-9-]+>=[0-9.]+", spec), spec


def _read(*parts):
    with open(os.path.join(ROOT, *parts), encoding="utf-8") as f:
        return f.read()


def _code(text):
    """Shell source without its comment lines."""
    return "\n".join(line for line in text.splitlines() if not line.lstrip().startswith("#"))


def test_update_checks_instead_of_reinstalling():
    update = _code(_read("update.sh"))
    assert "install-python-deps.sh\" --check" in update
    assert "install.sh\"" not in update.replace("install-python-deps.sh\"", "")
    assert "--force-reinstall" not in update
    assert "--clear" not in update


def test_install_and_update_share_the_deploy_checks():
    for script in ("install.sh", "update.sh"):
        code = _code(_read(script))
        assert 'source "$SCRIPT_DIR/scripts/deploy-common.sh"' in code, script
        for check in ("ensure_nltk_words", "ensure_apache_modules", "check_app_imports", "ensure_redis"):
            assert check in code, (script, check)


def _function(script, name):
    body = script[script.index(f"{name}() {{"):]
    return body[:body.index("\n}\n")]


def test_everything_runs_from_the_venv():
    """OPS.45: libraries go into planetGen's own venv on every platform, never into the system Python."""
    code = _code(_read("scripts", "install-python-deps.sh"))
    assert "-m venv" in _function(code, "make_venv")
    assert "-m venv" not in code.replace(_function(code, "make_venv"), "")
    for word in ("--break-system-packages", "--ignore-installed", "apt-get install -y --no-install-recommends ${"):
        assert word not in code, word
    # The only package apt is asked for is the one that makes a venv possible.
    assert re.findall(r"apt-get install[^\n]*", code) == ["apt-get install -y --no-install-recommends python3-venv || true"]
    assert "sys.path.insert" not in code  # no .pth written any more
    assert "(3, 12)" in code  # setup.py's python_requires


def test_the_apache_example_runs_the_app_in_the_venv():
    conf = _read("examples", "apache", "planetgen.conf.example")
    line = next(l for l in conf.splitlines() if l.strip().startswith("WSGIDaemonProcess planetgen-api"))
    assert "python-home=/usr/local/planetgen/venv" in line
    common = _code(_read("scripts", "deploy-common.sh"))
    assert 'python-home=$VENV' in _function(common, "check_apache_python_home")
    assert "check_apache_python_home" in _function(common, "ensure_apache_modules")


def test_the_services_and_scripts_use_the_venv_python():
    assert "ExecStart=/usr/local/planetgen/venv/bin/gunicorn" in _read("examples", "systemd", "planetgen-gunicorn.service")
    assert "ExecStart=/usr/local/planetgen/venv/bin/python -m planetgen.cli.orbits" in _read(
        "examples", "maintenance", "planetgen-orbits@.service")
    for script in ("install.sh", "update.sh"):
        code = _code(_read(script))
        assert code.index("install-python-deps.sh") < code.index("use_site_python"), script
    assert 'PYTHON="$VENV/bin/python"' in _function(_code(_read("scripts", "deploy-common.sh")), "use_site_python")


def _lock_pins():
    sys.path.insert(0, os.path.join(ROOT, "scripts"))
    try:
        import lock_pins
    finally:
        sys.path.pop(0)
    return lock_pins


def _key(version):
    return tuple(int(p) for p in re.findall(r"\d+", version)[:3])


def test_lock_satisfies_setup_py():
    """
    requirements.lock pins every runtime requirement at or above
    setup.py's floor, for every Python range it splits into. Rerun
    scripts/lock-requirements.sh when this fails.
    """
    lock_pins = _lock_pins()
    pins = lock_pins.read_lock(os.path.join(ROOT, "requirements.lock"))
    for spec in _setup_requirements():
        name, floor = spec.split(">=")
        versions = [v for n, v, _, _ in pins if n == lock_pins.normalize(name)]
        assert versions, f"{name} is not in requirements.lock"
        for version in versions:
            assert _key(version) >= _key(floor), f"{name}=={version} is below {floor}"


def test_every_locked_pin_has_hashes():
    lock_pins = _lock_pins()
    pins = lock_pins.read_lock(os.path.join(ROOT, "requirements.lock"))
    assert len(pins) > len(_setup_requirements())  # dependencies are locked too
    for name, version, _, hashes in pins:
        assert hashes and all(h.startswith("sha256:") for h in hashes), f"{name}=={version}"


def test_every_pip_install_uses_the_lock():
    """
    Every pip install of a library checks the lock's hashes; the only
    unhashed ones are pip itself and planetGen's own checkout (--no-deps).
    """
    code = _code(_read("scripts", "install-python-deps.sh"))
    installs = [line.strip() for line in code.splitlines() if "-m pip install" in line]
    assert installs
    for line in installs:
        assert ("--require-hashes" in line or "--upgrade pip" in line
                or ("--no-deps" in line and '"$SCRIPT_DIR"' in line)), line


MACOS_SCRIPTS = [
    ("install.sh",), ("update.sh",), ("scripts", "deploy-common.sh"),
    ("scripts", "install-python-deps.sh"), ("scripts", "lock-requirements.sh"),
    ("examples", "apache", "apache-identity.sh"), ("examples", "apache", "create-cache-dir.sh"),
    ("examples", "apache", "set-permissions.sh"), ("examples", "apache", "setup-debug-log.sh"),
    ("examples", "maintenance", "install-maintenance-timer.sh"),
]


@pytest.mark.parametrize("parts", MACOS_SCRIPTS, ids=lambda p: "/".join(p))
def test_bash_scripts_avoid_bash4_only_features(parts):
    """
    macOS ships bash 3.2: no mapfile/readarray, associative arrays, case
    modification or nameref, and "${a[@]}" of an empty array is an
    unbound-variable error under set -u before bash 4.4. Every array
    expansion is either guarded (${a[@]+"${a[@]}"}) or of an array that
    is never empty (listed here).
    """
    code = _code(_read(*parts))
    for word in ("mapfile", "readarray", "declare -A", "local -A", "declare -n", "local -n", ",,}", "^^}"):
        assert word not in code, (parts, word)
    never_empty = {"REQUIREMENTS", "SPECS", "enabled", "need", "databases"}
    for name in re.findall(r'(?<!\+)"\$\{(\w+)\[@\]', code):
        if name not in never_empty:
            unguarded = re.search(r'(?<!\+)"\$\{%s\[@\]' % name, code.replace('+"${%s[@]' % name, ""))
            assert not unguarded, (parts, name)


def test_server_lock_has_the_app_servers():
    lock_pins = _lock_pins()
    pins = lock_pins.read_lock(os.path.join(ROOT, "requirements-server.lock"))
    names = {name for name, _, _, _ in pins}
    assert {"gunicorn"} <= names
    for name, version, _, hashes in pins:
        assert hashes, f"{name}=={version}"


def _run_install_checkout(tmp_path, script_dir):
    """Runs the script's install_checkout with the venv's python a stand-in that
    records `-m pip` calls and hands everything else to this Python."""
    script = _read("scripts", "install-python-deps.sh")
    calls = tmp_path / "pip-calls"
    (tmp_path / "venv" / "bin").mkdir(parents=True)
    python = tmp_path / "venv" / "bin" / "python"
    # PYTHONPATH pins which checkout `import planetgen` finds (TEST.143): the editable install of whatever
    # checkout this machine was set up from would otherwise decide, and a test run from another checkout
    # or a git worktree saw "planetgen is not imported from here" and reinstalled.
    python.write_text(f'#!/bin/sh\nif [ "$1" = -m ] && [ "$2" = pip ]; then echo "$*" >> "{calls}"; exit 0; fi\n'
                      f'PYTHONPATH="{os.path.join(os.path.realpath(ROOT), "src")}" exec "{sys.executable}" "$@"\n')
    python.chmod(0o755)
    command = tmp_path / "bin" / "planetgen"
    harness = "\n".join([
        "set -euo pipefail", f'VENV="{tmp_path / "venv"}"', f'SCRIPT_DIR="{script_dir}"',
        f'PLANETGEN_COMMAND="{command}"', _function(script, "install_checkout") + "\n}\n", "install_checkout",
    ])
    result = subprocess.run(["bash", "-c", harness], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return (calls.read_text().splitlines() if calls.exists() else []), command


def test_install_checkout_leaves_a_current_editable_install_alone(tmp_path):
    """update.sh reinstalls nothing when planetgen already imports from
    this checkout; the `planetgen` command links to pip's console script."""
    import sysconfig
    if not os.access(os.path.join(sysconfig.get_path("scripts"), "planetgen"), os.X_OK):
        pytest.skip("no planetgen console script")
    calls, command = _run_install_checkout(tmp_path, os.path.realpath(ROOT))
    assert calls == []
    assert os.readlink(command) == os.path.join(sysconfig.get_path("scripts"), "planetgen")


def test_install_checkout_installs_another_checkout_editable(tmp_path):
    """A checkout planetgen doesn't import from is pip-installed editable,
    without its dependencies (they come from the lock)."""
    other = tmp_path / "checkout"
    (other / "src").mkdir(parents=True)
    calls, _ = _run_install_checkout(tmp_path, other)
    assert len(calls) == 1
    assert calls[0].endswith(f"--no-deps -e {other}")
    assert "--break-system-packages" not in calls[0]


def test_the_web_app_never_imports_numba():
    """OPS.42: numba's default cache needs a writable code tree or home folder, which Apache's web user has neither
    of, so importing a compiled function there fails ("cannot cache function: no locator available"). The web app
    must stay clear of numba: only the workers, the settle step and the orbit update import it."""
    code = ("import sys\n"
            "from planetgen.web.app import create_app\n"
            "from planetgen.api import routes  # noqa: F401\n"
            "import planetgen.web.views, planetgen.web.system_pages  # noqa: F401\n"
            "assert 'numba' not in sys.modules, 'numba was imported by the web app'\n")
    done = subprocess.run([sys.executable, "-c", code], cwd=os.path.join(ROOT, "src"), capture_output=True, text=True)
    assert done.returncode == 0, done.stderr[-1500:]


def test_cryptography_is_a_requirement_for_mysql_8_logins():
    """The first real update.sh on a venv failed at the migration with "'cryptography' package is required for
    sha256_password or caching_sha2_password auth methods": pymysql needs it to log in to MySQL 8's default accounts."""
    assert "cryptography>=42.0.0" in _setup_requirements()
    assert "cryptography>=42.0.0" in _script_requirements()
    lock = _lock_pins().read_lock(os.path.join(ROOT, "requirements-server.lock"))
    assert "cryptography" in {name for name, _, _, _ in lock}


def test_the_library_check_fails_a_pymysql_that_cannot_use_cryptography():
    """A venv with pymysql but no cryptography is "broken" at the library step, not at the migration."""
    code = ("import sys; sys.modules['cryptography'] = None; sys.argv = ['probe', 'pymysql>=1.1.1']\n"
            f"exec(open({os.path.join(ROOT, 'scripts', 'probe_requirements.py')!r}).read())\n")
    done = subprocess.run([sys.executable, "-I", "-c", code], capture_output=True, text=True)
    assert done.stdout.startswith("broken pymysql>=1.1.1 ImportError:pymysql_has_no_cryptography"), done.stdout + done.stderr
