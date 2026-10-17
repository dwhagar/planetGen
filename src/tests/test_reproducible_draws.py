# tests/test_reproducible_draws.py

"""
GEN.56: every random draw in generation comes from the derived seeds,
through `planetgen.util.draw`'s helpers built on `random()` alone (see
docs/design/reproducible-galaxies.md, section 5).

- The helpers' own recipes are pinned, so a change to one shows up here
  rather than as a quietly different galaxy.
- A scan of the generation packages fails on any other random source:
  the stdlib `random` module, `secrets`, `uuid`, `os.urandom`,
  `SystemRandom`, `numpy.random`, or a stream seeded from the clock.
- One small sector generated under two `PYTHONHASHSEED` values and two
  locales comes out the same, so nothing depends on set or string-keyed
  dict iteration order or on locale sorting.
"""

import ast
import locale
import math
import os
import pathlib
import subprocess
import sys
import threading

import pytest

from planetgen.util import draw

SRC = pathlib.Path(__file__).resolve().parents[1]
PACKAGE = SRC / "planetgen"

GENERATION_PACKAGES = ("galaxy", "generation", "names", "physics", "population", "util")
"""The packages whose code draws for generation. `api/`, `web/`, `admin/`,
`queue/` and `db/` hold login, CSRF, API keys, job ids and retry pauses,
which are not generation draws (and call generators only through
`draw`)."""

ALLOWED = {
    "util/draw.py": {"random", "secrets", "numpy.random"},  # the one wrapper, the run stream's first seed, and the
    # hand-over of the stream's state to NumPy's MT19937 for a bulk draw (PERF.91)
    "galaxy/seed.py": {"secrets"},  # a brand-new galaxy seed (`new_seed`)
}
"""Modules allowed one of the forbidden imports, and which."""

FORBIDDEN_MODULES = {"random", "secrets", "uuid"}
SEEDING_CALLS = {"bound", "set_run_seed", "Stream"}


def _generation_modules():
    for name in GENERATION_PACKAGES:
        yield from sorted((PACKAGE / name).rglob("*.py"))


def _problems(path, relative=None):
    relative = relative or path.relative_to(PACKAGE).as_posix()
    allowed = ALLOWED.get(relative, set())
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    found = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                root = alias.name.split(".")[0]
                if (root in FORBIDDEN_MODULES and root not in allowed) or (
                        alias.name.startswith("numpy.random") and "numpy.random" not in allowed):
                    found.append(f"{relative}:{node.lineno} imports {alias.name}")
        elif isinstance(node, ast.ImportFrom):
            module = node.module or ""
            root = module.split(".")[0]
            if node.level == 0 and ((root in FORBIDDEN_MODULES and root not in allowed)
                                    or (module.startswith("numpy.random") and "numpy.random" not in allowed)):
                found.append(f"{relative}:{node.lineno} imports from {module}")
        elif isinstance(node, ast.Attribute):
            if node.attr in ("urandom", "SystemRandom", "uuid4", "token_bytes", "randbits") \
                    and not (node.attr in ("token_bytes", "randbits") and "secrets" in allowed):
                found.append(f"{relative}:{node.lineno} uses .{node.attr}")
            if node.attr == "random" and isinstance(node.value, ast.Attribute) and node.value.attr == "random" \
                    and isinstance(node.value.value, ast.Name) and node.value.value.id in ("np", "numpy"):
                found.append(f"{relative}:{node.lineno} uses numpy.random")
        elif isinstance(node, ast.Call):
            func = node.func
            name = func.attr if isinstance(func, ast.Attribute) else getattr(func, "id", None)
            if name in SEEDING_CALLS:
                for arg in ast.walk(ast.Tuple(elts=list(node.args) + [k.value for k in node.keywords])):
                    if isinstance(arg, ast.Name) and arg.id == "time" \
                            or isinstance(arg, ast.Attribute) and isinstance(arg.value, ast.Name) and arg.value.id == "time":
                        found.append(f"{relative}:{node.lineno} seeds a stream from the clock")
                        break
    return found


def test_generation_draws_only_through_the_draw_module():
    problems = [problem for path in _generation_modules() for problem in _problems(path)]
    assert not problems, "\n".join(problems)


def test_the_scan_catches_a_forbidden_source(tmp_path):
    bad = tmp_path / "probe.py"
    bad.write_text("import random\nimport os\nfrom planetgen.util import draw\nimport time\n"
                   "x = os.urandom(4)\ns = draw.Stream(time.time())\n", encoding="utf-8")
    problems = _problems(bad, "generation/probe.py")
    assert any("imports random" in p for p in problems)
    assert any(".urandom" in p for p in problems)
    assert any("clock" in p for p in problems)


# --- Maths that depends on the Python or the BLAS library (GEN.135) ---------------

NOT_DETERMINISTIC_EXEMPT = ("web/", "db/migrations/")
"""Where `math.hypot` and `math.dist` are allowed: drawing code, where the last bit does not matter, and migration
revisions, which are never edited once they have run."""


def _maths_problems(path, relative):
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    found = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Attribute) and node.attr in ("hypot", "dist") \
                and isinstance(node.value, ast.Name) and node.value.id == "math":
            found.append(f"{relative}:{node.lineno} uses math.{node.attr} (use planetgen.util.detmath.{node.attr})")
        elif isinstance(node, ast.ImportFrom) and node.module == "math" and node.level == 0 \
                and any(alias.name in ("hypot", "dist") for alias in node.names):
            found.append(f"{relative}:{node.lineno} imports hypot or dist from math (use planetgen.util.detmath)")
        elif isinstance(node, ast.BinOp) and isinstance(node.op, ast.MatMult) \
                or isinstance(node, ast.AugAssign) and isinstance(node.op, ast.MatMult):
            found.append(f"{relative}:{node.lineno} uses @, which runs through BLAS (use planetgen.util.detmath.matmul)")
        elif isinstance(node, ast.Attribute) and node.attr in ("dot", "matmul", "vdot", "inner", "tensordot") \
                and isinstance(node.value, ast.Name) and node.value.id in ("np", "numpy"):
            found.append(f"{relative}:{node.lineno} uses numpy.{node.attr}, which runs through BLAS")
    return found


def test_no_version_or_library_dependent_maths_in_the_package():
    problems = []
    for path in sorted(PACKAGE.rglob("*.py")):
        relative = path.relative_to(PACKAGE).as_posix()
        if relative.startswith(NOT_DETERMINISTIC_EXEMPT) or relative == "util/detmath.py":
            continue
        problems.extend(_maths_problems(path, relative))
    assert not problems, "\n".join(problems)


def test_the_maths_scan_catches_what_it_is_for(tmp_path):
    bad = tmp_path / "probe.py"
    bad.write_text("import math\nimport numpy as np\nfrom math import hypot\n"
                   "a = math.hypot(1.0, 2.0)\nb = math.dist((0, 0), (1, 1))\nc = m @ n\nd = np.dot(m, n)\n",
                   encoding="utf-8")
    problems = _maths_problems(bad, "galaxy/probe.py")
    assert len(problems) == 5, problems
    ok = tmp_path / "fine.py"
    ok.write_text("import math\nx = math.sqrt(2.0) + math.fsum([1.0])\n", encoding="utf-8")
    assert _maths_problems(ok, "galaxy/fine.py") == []


def test_the_helpers_are_the_plain_formulas():
    from planetgen.util import detmath
    assert detmath.hypot(3.0, 4.0) == 5.0
    assert detmath.hypot(0.1, 0.2, 0.3) == math.sqrt(0.1 * 0.1 + 0.2 * 0.2 + 0.3 * 0.3)
    assert detmath.hypot(1.0, 2.0, 3.0, 4.0) == math.sqrt(((1.0 + 4.0) + 9.0) + 16.0)
    assert detmath.hypot() == 0.0
    assert detmath.dist((1.0, 2.0, 3.0), (4.0, 6.0, 3.0)) == 5.0
    assert detmath.dist((1.0, 2.0), (4.0, 6.0)) == 5.0
    assert detmath.dist((1.0,), (4.0,)) == 3.0
    assert detmath.dist((1.0, 2.0, 3.0, 4.0), (1.0, 2.0, 3.0, 8.0)) == 4.0
    with pytest.raises(ValueError):
        detmath.dist((1.0, 2.0), (1.0, 2.0, 3.0))
    assert detmath.fsum([0.1] * 10) == 1.0


def test_matmul_is_the_left_to_right_sum_of_products():
    import numpy as np
    from planetgen.util import detmath
    a = np.array([[0.1, 0.2, 0.3], [1.5, -2.5, 3.5]])
    b = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
    expected = np.array([[(0.1 * 1.0 + 0.2 * 3.0) + 0.3 * 5.0, (0.1 * 2.0 + 0.2 * 4.0) + 0.3 * 6.0],
                         [(1.5 * 1.0 + -2.5 * 3.0) + 3.5 * 5.0, (1.5 * 2.0 + -2.5 * 4.0) + 3.5 * 6.0]])
    assert (detmath.matmul(a, b) == expected).all()
    assert detmath.matmul(np.array([1.0, 2.0, 3.0]), np.eye(3)).tolist() == [1.0, 2.0, 3.0]
    assert np.allclose(detmath.matmul(a, b), a @ b)


# --- The helpers' recipes ------------------------------------------------------

def test_the_recipes_are_pinned():
    """Built on `random()` alone, so these values hold on every Python
    release; a change here changes every galaxy and needs a note."""
    stream = draw.Stream(12345)
    first = stream.random()
    stream = draw.Stream(12345)
    assert stream.uniform(10.0, 20.0) == 10.0 + 10.0 * first
    stream = draw.Stream(12345)
    assert stream.randrange(10) == int(first * 10)
    stream = draw.Stream(12345)
    assert stream.choice("abcdefghij") == "abcdefghij"[int(first * 10)]


def test_the_helpers_stay_in_range():
    stream = draw.Stream(7)
    for _ in range(2000):
        assert 3 <= stream.randint(3, 5) <= 5
        assert 0 <= stream.randrange(4) < 4
        assert -2.0 <= stream.uniform(-2.0, 2.0) <= 2.0
        assert 0 <= stream.getrandbits(63) < 2 ** 63
    assert {stream.randint(0, 2) for _ in range(500)} == {0, 1, 2}
    with pytest.raises(ValueError):
        stream.randrange(0)
    with pytest.raises(IndexError):
        stream.choice([])
    with pytest.raises(IndexError):
        stream.choices([], weights=[])


def test_shuffle_and_sample_are_permutations():
    stream = draw.Stream(3)
    items = list(range(20))
    stream.shuffle(items)
    assert sorted(items) == list(range(20)) and items != list(range(20))
    picked = stream.sample(range(50), 10)
    assert len(set(picked)) == 10 and all(0 <= p < 50 for p in picked)


def test_choices_follow_their_weights_and_gauss_its_moments():
    stream = draw.Stream(11)
    picks = stream.choices("ab", weights=[1, 3], k=20000)
    assert picks.count("b") / len(picks) == pytest.approx(0.75, abs=0.02)
    values = [stream.gauss(5.0, 2.0) for _ in range(20000)]
    mean = sum(values) / len(values)
    spread = (sum((v - mean) ** 2 for v in values) / len(values)) ** 0.5
    assert mean == pytest.approx(5.0, abs=0.06) and spread == pytest.approx(2.0, abs=0.06)


def test_a_bound_stream_is_this_threads_only():
    seen = []
    with draw.bound(99) as stream:
        assert draw.current() is stream
        worker = threading.Thread(target=lambda: seen.append(draw.current()))
        worker.start()
        worker.join()
    assert seen[0] is not stream
    assert draw.current() is not stream


def test_a_unit_draws_the_same_whatever_ran_before():
    with draw.bound(5):
        first = [draw.random() for _ in range(3)]
    draw.random()
    with draw.bound(5):
        assert [draw.random() for _ in range(3)] == first


# --- Hash seeds and locales ----------------------------------------------------------

def _available_locales():
    found = []
    for name in ("C", "C.UTF-8", "C.utf8", "en_US.UTF-8"):
        try:
            locale.setlocale(locale.LC_ALL, name)
        except locale.Error:
            continue
        finally:
            locale.setlocale(locale.LC_ALL, "")
        found.append(name)
    return found


def _probe(hash_seed, locale_name):
    env = dict(os.environ, PYTHONHASHSEED=str(hash_seed), LC_ALL=locale_name,
               PYTHONPATH=os.pathsep.join([str(SRC), os.environ.get("PYTHONPATH", "")]))
    done = subprocess.run([sys.executable, str(SRC / "tests" / "reproducible_sector_probe.py"), "5"],
                          env=env, capture_output=True, text=True, timeout=600)
    assert done.returncode == 0, done.stderr[-2000:]
    return done.stdout.strip().splitlines()[-1]


@pytest.mark.slow
def test_a_sector_does_not_depend_on_hash_seed_or_locale():
    locales = _available_locales()
    names = [locales[0], locales[-1]] if locales else ["C"]
    digests = {(seed, name): _probe(seed, name) for seed in (0, 1) for name in dict.fromkeys(names)}
    assert len(set(digests.values())) == 1, digests
