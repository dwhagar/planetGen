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
    "util/draw.py": {"random", "secrets"},  # the one wrapper, and the run stream's first seed
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
                if (root in FORBIDDEN_MODULES and root not in allowed) or alias.name.startswith("numpy.random"):
                    found.append(f"{relative}:{node.lineno} imports {alias.name}")
        elif isinstance(node, ast.ImportFrom):
            module = node.module or ""
            root = module.split(".")[0]
            if node.level == 0 and ((root in FORBIDDEN_MODULES and root not in allowed)
                                    or module.startswith("numpy.random")):
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
