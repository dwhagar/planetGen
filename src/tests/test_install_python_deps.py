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
        for check in ("ensure_nltk_words", "ensure_apache_modules", "check_app_imports"):
            assert check in code, (script, check)


def test_check_mode_keeps_the_venv():
    script = _read("scripts", "install-python-deps.sh")
    check = script[script.index("check_requirements() {"):]
    check = check[:check.index("\n}\n")]
    assert "venv_install keep" in check
    assert "venv_install clear" not in check
    assert "--force-reinstall" not in check
