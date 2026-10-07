"""
scripts/install-python-deps.sh keeps its own list of runtime requirements
(pip requirement plus the apt package that provides it) for installing on
an externally managed Python. This checks that list matches setup.py's
install_requires plus its 'api' extra, so a dependency added or bumped in
setup.py can't be silently missed on a managed-Python install.

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
    return dict(re.findall(r'"(\S+) (python3-\S+)"', block))


def test_script_requirements_match_setup_py():
    assert {r.lower() for r in _script_requirements()} == _setup_requirements()


def test_every_requirement_has_a_floor_and_apt_package():
    for spec, package in _script_requirements().items():
        assert re.fullmatch(r"[a-z0-9-]+>=[0-9.]+", spec), spec
        assert package.startswith("python3-"), package


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


def test_no_virtual_environment():
    """Boss's call: on Linux, libraries go into the system Python, apt
    first. Only macOS's path (install_venv, picked on Darwin) makes one."""
    code = _code(_read("scripts", "install-python-deps.sh"))
    venv = _function(code, "install_venv")
    assert "-m venv" in venv
    assert "-m venv" not in code.replace(venv, "")
    assert 'if [[ "$(uname -s)" == Darwin ]]; then MODE=venv' in code
    assert "sys.path.insert" not in code  # no .pth written any more


def test_pip_never_removes_apt_files():
    """
    pip must never uninstall or upgrade in place on a managed Python: for
    an apt package with egg-info metadata (python3-pymysql) that deletes
    apt's own files. It resolves with --dry-run --report and installs the
    pinned result alongside with --ignore-installed --no-deps.
    """
    pip = _function(_code(_read("scripts", "install-python-deps.sh")), "pip_install_system")
    resolver = _read("scripts", "lock_pins.py")
    assert '"--dry-run"' in resolver and '"--report"' in resolver
    assert "--ignore-installed --no-deps" in pip
    assert "--upgrade" not in pip
    # The only uninstall is of copies outside apt's directory.
    assert '"/usr/lib/python3/"' in pip


def test_managed_paths_try_apt_first():
    script = _code(_read("scripts", "install-python-deps.sh"))
    for name in ("install_managed", "check_requirements"):
        body = _function(script, name)
        assert body.index("apt_install") < body.index("pip_install_needed"), name
        assert "remove_legacy_venv" in body, name


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
    unhashed ones are pip itself, planetGen's own checkout (--no-deps) and
    the fallback for a pip too old to report, which is held to the lock's
    versions.
    """
    code = _code(_read("scripts", "install-python-deps.sh"))
    installs = [line.strip() for line in code.splitlines() if "-m pip install" in line]
    assert installs
    for line in installs:
        assert ("--require-hashes" in line or "--upgrade pip" in line
                or ("--no-deps" in line and '"$SCRIPT_DIR"' in line)
                or '-c "$hashed"' in line), line


def test_resolve_writes_hashed_pins(tmp_path, monkeypatch):
    lock_pins = _lock_pins()
    lock = tmp_path / "requirements.lock"
    lock.write_text(
        "# header\n"
        "flask==3.1.3 \\\n    --hash=sha256:aa \\\n    --hash=sha256:bb\n"
        "    # via planetgen\n"
        "click==8.1.8 ; python_full_version < '3.10' \\\n    --hash=sha256:cc\n"
        "click==8.5.0 ; python_full_version >= '3.10' \\\n    --hash=sha256:dd\n",
        encoding="utf-8",
    )
    seen = []

    def fake_dry_run(flags, requirements, lines):
        seen.append(lines)
        return [("flask", "3.1.3"), ("click", "8.5.0")]

    monkeypatch.setattr(lock_pins, "_dry_run", fake_dry_run)
    monkeypatch.setattr(lock_pins, "_installed", lambda name: name == "click")
    out = tmp_path / "out.txt"
    assert lock_pins.resolve(str(lock), str(out), [], ["flask>=3.0.3"]) == 0
    # click was installed, so it isn't held on the first pass; pip chose
    # to replace it anyway, so the second pass holds it to the lock.
    assert seen[0] == ["flask==3.1.3"]
    assert "click==8.5.0 ; python_full_version >= '3.10'" in seen[1]
    assert out.read_text(encoding="utf-8").splitlines() == [
        "flask==3.1.3 --hash=sha256:aa --hash=sha256:bb",
        "click==8.5.0 --hash=sha256:dd",
    ]


def test_resolve_refuses_unlocked_versions(tmp_path, monkeypatch):
    lock_pins = _lock_pins()
    lock = tmp_path / "requirements.lock"
    lock.write_text("flask==3.1.3 \\\n    --hash=sha256:aa\n", encoding="utf-8")
    monkeypatch.setattr(lock_pins, "_dry_run", lambda *a: [("flask", "3.2.0")])
    monkeypatch.setattr(lock_pins, "_installed", lambda name: False)
    assert lock_pins.resolve(str(lock), str(tmp_path / "out"), [], ["flask"]) == 1


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
    never_empty = {"REQUIREMENTS", "SPECS", "NUMPY_STACK", "enabled", "need", "pins", "owned", "databases", "wanted", "lacking", "before"}
    for name in re.findall(r'(?<!\+)"\$\{(\w+)\[@\]', code):
        if name not in never_empty:
            unguarded = re.search(r'(?<!\+)"\$\{%s\[@\]' % name, code.replace('+"${%s[@]' % name, ""))
            assert not unguarded, (parts, name)


def test_powershell_requirements_match_setup_py():
    """scripts/deploy-common.ps1 checks the same requirements on Windows."""
    text = _read("scripts", "deploy-common.ps1")
    block = re.search(r"\$script:Requirements = @\((.*?)\)", text, re.S).group(1)
    assert set(re.findall(r'"([^"]+)"', block)) == _setup_requirements()
    server = re.search(r'\$script:ServerRequirement = "([^"]+)"', text).group(1)
    assert server + '; sys_platform == "win32"' in _read("setup.py")


POWERSHELL_SCRIPTS = [("install.ps1",), ("update.ps1",), ("scripts", "deploy-common.ps1"),
                      ("examples", "maintenance", "install-maintenance-task.ps1")]


@pytest.mark.parametrize("parts", POWERSHELL_SCRIPTS, ids=lambda p: "/".join(p))
def test_powershell_scripts_are_ascii(parts):
    """Windows PowerShell 5.1 reads a BOM-less script as the ANSI code
    page, so anything outside ASCII would come out garbled."""
    _read(*parts).encode("ascii")


def test_powershell_installers_share_the_steps():
    """install.ps1 and update.ps1 do what install.sh and update.sh do,
    through the shared deploy-common.ps1."""
    for script in ("install.ps1", "update.ps1"):
        text = _read(script)
        assert '. (Join-Path $Root "scripts\\deploy-common.ps1")' in text, script
        for step in ("Install-PythonDeps", "Install-NltkWords", "Invoke-MigrateOrReset",
                     "New-RuntimeDirs", "Set-PlanetGenPermissions", "Test-AppImports"):
            assert step in text, (script, step)
    assert "Install-PythonDeps -Check" in _read("update.ps1")
    assert "--force-reinstall" not in _read("update.ps1")
    common = _read("scripts", "deploy-common.ps1")
    assert "--require-hashes -r (Join-Path $Root \"requirements-server.lock\")" in common
    # The same migrate-or-delete question as install.sh: y/N, 30 seconds.
    assert "[y/N] (default N in 30s): \" 30" in common


def test_server_lock_has_the_app_servers():
    lock_pins = _lock_pins()
    pins = lock_pins.read_lock(os.path.join(ROOT, "requirements-server.lock"))
    names = {name for name, _, _, _ in pins}
    assert {"gunicorn", "waitress"} <= names
    for name, version, _, hashes in pins:
        assert hashes, f"{name}=={version}"


# --- OPS.26: the NumPy stack comes from one place ------------------------------

NUMPY_LINKED = {"numpy", "scipy", "astropy", "scikit-image"}
"""Runtime requirements with compiled code built against NumPy's C API."""


def _numpy_stack():
    text = _read("scripts", "install-python-deps.sh")
    block = re.search(r"^NUMPY_STACK=\((.*?)\)", text, re.M).group(1)
    names = re.search(r'^NUMPY_STACK_NAMES=" (.*) "', text, re.M).group(1).split()
    return re.findall(r'"(\S+)"', block), names


def test_the_numpy_stack_lists_every_numpy_built_requirement():
    specs, names = _numpy_stack()
    assert NUMPY_LINKED <= {spec.split(">=")[0] for spec in _script_requirements()}
    assert {spec.split(">=")[0] for spec in specs} == NUMPY_LINKED | {"pyerfa"} == set(names)
    # Each at setup.py's own floor.
    for spec in specs:
        if spec.split(">=")[0] in NUMPY_LINKED:
            assert spec in _script_requirements(), spec
    lock_pins = _lock_pins()
    locked = {name for name, _, _, _ in lock_pins.read_lock(os.path.join(ROOT, "requirements.lock"))}
    assert "pyerfa" in locked


def test_managed_paths_install_the_numpy_stack_whole():
    script = _code(_read("scripts", "install-python-deps.sh"))
    for name in ("install_managed", "check_requirements"):
        body = _function(script, name)
        assert "pip_install_needed" in body and "pip_install_system" not in body, name
    pip = _function(script, "pip_install_system")
    assert "--ignore-installed" in pip and "return 4" in pip


def _run_needed(tmp_path, specs, rest_status):
    """Runs the script's pip_install_needed with a stub pip_install_system
    that records each call and returns `rest_status` for one without
    --whole."""
    script = _read("scripts", "install-python-deps.sh")
    stack = re.search(r"^NUMPY_STACK=.*$", script, re.M).group(0)
    names = re.search(r"^NUMPY_STACK_NAMES=.*$", script, re.M).group(0)
    needed = _function(script, "pip_install_needed") + "\n}\n"
    calls = tmp_path / "calls"
    harness = "\n".join([
        "set -euo pipefail", stack, names,
        f'pip_install_system() {{ echo "$*" >> "{calls}"; [[ "$1" == --whole ]] || return {rest_status}; }}',
        needed, "pip_install_needed " + " ".join(f'"{spec}"' for spec in specs),
    ])
    subprocess.run(["bash", "-c", harness], check=True, capture_output=True, text=True)
    return calls.read_text().splitlines() if calls.exists() else []


@pytest.mark.skipif(sys.platform == "win32", reason="bash")
def test_a_stack_library_brings_the_whole_stack(tmp_path):
    specs, _names = _numpy_stack()
    calls = _run_needed(tmp_path, ["redis>=5.0.0", "scipy>=1.13.0"], 0)
    assert calls == ["redis>=5.0.0", "--whole " + " ".join(specs)]


@pytest.mark.skipif(sys.platform == "win32", reason="bash")
def test_a_library_that_would_pull_in_numpy_brings_the_whole_stack(tmp_path):
    specs, _names = _numpy_stack()
    calls = _run_needed(tmp_path, ["redis>=5.0.0"], 4)
    assert calls == ["redis>=5.0.0", "--whole " + " ".join(specs)]


@pytest.mark.skipif(sys.platform == "win32", reason="bash")
def test_other_libraries_leave_the_stack_alone(tmp_path):
    assert _run_needed(tmp_path, ["redis>=5.0.0"], 0) == ["redis>=5.0.0"]


def _run_install_checkout(tmp_path, script_dir, mode="unmanaged"):
    """Runs the script's install_checkout with PYTHON a stand-in that
    records `-m pip` calls and hands everything else to this Python."""
    script = _read("scripts", "install-python-deps.sh")
    calls = tmp_path / "pip-calls"
    python = tmp_path / "python"
    python.write_text(f'#!/bin/sh\nif [ "$1" = -m ] && [ "$2" = pip ]; then echo "$*" >> "{calls}"; exit 0; fi\n'
                      f'exec "{sys.executable}" "$@"\n')
    python.chmod(0o755)
    command = tmp_path / "bin" / "planetgen"
    harness = "\n".join([
        "set -euo pipefail", f'PYTHON="{python}"', f'SCRIPT_DIR="{script_dir}"', f'MODE={mode}',
        f'PLANETGEN_COMMAND="{command}"', _function(script, "install_checkout") + "\n}\n", "install_checkout",
    ])
    result = subprocess.run(["bash", "-c", harness], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return (calls.read_text().splitlines() if calls.exists() else []), command


def _installed_from_this_checkout():
    try:
        import planetgen
    except ImportError:
        return False
    src = os.path.realpath(os.path.join(ROOT, "src"))
    return os.path.dirname(os.path.dirname(os.path.realpath(planetgen.__file__))) == src


@pytest.mark.skipif(sys.platform == "win32", reason="bash")
@pytest.mark.skipif(not _installed_from_this_checkout(), reason="needs `pip install -e .` of this checkout")
def test_install_checkout_leaves_a_current_editable_install_alone(tmp_path):
    """update.sh reinstalls nothing when planetgen already imports from
    this checkout; the `planetgen` command links to pip's console script."""
    import sysconfig
    if not os.access(os.path.join(sysconfig.get_path("scripts"), "planetgen"), os.X_OK):
        pytest.skip("no planetgen console script")
    calls, command = _run_install_checkout(tmp_path, os.path.realpath(ROOT))
    assert calls == []
    assert os.readlink(command) == os.path.join(sysconfig.get_path("scripts"), "planetgen")


@pytest.mark.skipif(sys.platform == "win32", reason="bash")
@pytest.mark.parametrize("mode, flag", [("unmanaged", False), ("managed", True)])
def test_install_checkout_installs_another_checkout_editable(tmp_path, mode, flag):
    """A checkout planetgen doesn't import from is pip-installed editable,
    without its dependencies (they come from the lock); a PEP 668 Python
    needs --break-system-packages."""
    other = tmp_path / "checkout"
    (other / "src").mkdir(parents=True)
    calls, _ = _run_install_checkout(tmp_path, other, mode)
    assert len(calls) == 1
    assert calls[0].endswith(f"--no-deps -e {other}")
    assert ("--break-system-packages" in calls[0]) == flag
