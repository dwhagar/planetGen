# tests/test_subprocess_checkout.py

"""A program a test starts imports `planetgen` from the checkout under test (TEST.145).

The editable install points at one checkout. Run the tests from another (a git worktree) and every child process
(the generate CLI, a worker, a script run with `sys.executable`) would import the install's code instead, so a
test would pass or fail on code it is not testing. `conftest.py` puts this checkout's `src` first on PYTHONPATH
for the whole run; these tests fail if a child still finds planetgen elsewhere.
"""

import os
import subprocess
import sys

SRC = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))


def _imported_from(env=None, args=("-c", "import planetgen; print(planetgen.__file__)")):
    result = subprocess.run([sys.executable, *args], capture_output=True, text=True, timeout=60, env=env, check=True)
    return os.path.realpath(result.stdout.strip())


def test_a_child_process_imports_planetgen_from_this_checkout():
    assert _imported_from().startswith(SRC + os.sep)


def test_a_child_given_a_copy_of_the_environment_does_too():
    assert _imported_from(env=dict(os.environ)).startswith(SRC + os.sep)


def test_a_child_started_as_a_module_does_too():
    result = subprocess.run([sys.executable, "-m", "planetgen.cli.generate", "--help"], capture_output=True,
                            text=True, timeout=120, check=False)
    assert "planetgen" in (result.stdout + result.stderr).lower()


def test_the_checkout_comes_first_on_the_path():
    assert os.environ["PYTHONPATH"].split(os.pathsep)[0] == SRC
