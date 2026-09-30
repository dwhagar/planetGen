"""
scripts/install-python-deps.sh keeps its own list of runtime requirements
for planetGen's venv. This checks that list matches setup.py's
install_requires plus its 'api' extra, so a dependency added or bumped in
setup.py can't be silently missed on a server. It also pins the
deployment rules: libraries live in the dedicated venv, never in the
system Python, and update.sh checks rather than reinstalls.

Run with: pytest src/tests/test_install_python_deps.py
"""
import ast
import os
import re

ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
VENV = "/opt/planetgen/venv"


def _read(*parts):
    with open(os.path.join(ROOT, *parts), encoding="utf-8") as f:
        return f.read()


def _code(text):
    """Shell source without its comment lines."""
    return "\n".join(line for line in text.splitlines() if not line.lstrip().startswith("#"))


def _setup_requirements():
    tree = ast.parse(_read("setup.py"))
    call = next(
        node for node in ast.walk(tree)
        if isinstance(node, ast.Call) and getattr(node.func, "id", None) == "setup"
    )
    kwargs = {kw.arg: kw.value for kw in call.keywords}
    reqs = ast.literal_eval(kwargs["install_requires"])
    reqs += ast.literal_eval(kwargs["extras_require"])["api"]
    return {r.lower() for r in reqs}


def _script_requirements():
    text = _read("scripts", "install-python-deps.sh")
    block = re.search(r"^REQUIREMENTS=\((.*?)^\)", text, re.S | re.M).group(1)
    return re.findall(r'"(\S+)"', block)


def test_script_requirements_match_setup_py():
    assert {r.lower() for r in _script_requirements()} == _setup_requirements()


def test_every_requirement_has_a_floor():
    for spec in _script_requirements():
        assert re.fullmatch(r"[a-z0-9-]+>=[0-9.]+", spec), spec


def test_system_python_is_never_changed():
    """PEP 668: nothing pip-installs into the system Python."""
    for parts in (("install.sh",), ("update.sh",), ("scripts", "install-python-deps.sh"),
                  ("scripts", "deploy-common.sh")):
        code = _code(_read(*parts))
        assert "--break-system-packages" not in code, parts
        assert not re.search(r'"\$PYTHON" -m pip install', code), parts


def test_libraries_go_into_the_dedicated_venv():
    code = _code(_read("scripts", "install-python-deps.sh"))
    assert f"/opt/planetgen/venv" in code
    assert '"$PYTHON" -m venv' in code
    assert "--system-site-packages" not in code
    assert '"$VENV_PYTHON" -m pip install' in code


def test_update_checks_instead_of_reinstalling():
    update = _code(_read("update.sh"))
    assert 'install-python-deps.sh" --check' in update
    assert 'install.sh"' not in update.replace('install-python-deps.sh"', "")
    assert "--force-reinstall" not in update


def test_install_and_update_share_the_deploy_checks():
    for script in ("install.sh", "update.sh"):
        code = _code(_read(script))
        assert 'source "$SCRIPT_DIR/scripts/deploy-common.sh"' in code, script
        for check in ("ensure_nltk_words", "ensure_apache_modules", "check_app_imports"):
            assert check in code, (script, check)
        # Everything after the deps step runs with the venv's Python.
        assert 'PYTHON="${PLANETGEN_VENV_DIR:-/opt/planetgen/venv}/bin/python"' in code, script


def test_apache_and_services_run_the_venv():
    conf = _read("examples", "apache", "planetgen.conf.example")
    assert re.search(rf"WSGIDaemonProcess planetgen-api python-home={VENV} ", conf)
    orbits = _read("examples", "maintenance", "planetgen-orbits@.service")
    assert f"ExecStart={VENV}/bin/python src/updateOrbits.py" in orbits
