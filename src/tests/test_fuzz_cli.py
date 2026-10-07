# tests/test_fuzz_cli.py

"""
Property-based fuzzing of `generate.py`'s command line: `build_parser()`/
`process_args()` on random argv lists assembled from the parser's *real*
option names (read off the built parser, so a new option is fuzzed the
day it's added) mixed with hostile values and junk tokens, plus full
`main()` runs of every subcommand -- `system`, `sector`, `phenomenon`,
`plan`, and `galaxy` against a tiny, flat galaxy -- with tiny, huge,
negative, zero and non-finite numeric options.

The contract under test: every invocation ends in a clean `SystemExit`
(argparse's own code 2 for a rejected command line, 0 for --help/
--version, or the CLI's own documented code 1 for a runtime refusal such
as "outside the galaxy" or "no room") or a successful run -- never a
traceback, and never a hang.

`test_bughunt_generation_fuzz.py`/`test_bughunt_cli_edges.py` already
run `system`/`sector`/`phenomenon` with random *valid* flag combinations;
this file goes after the parser itself, the values of every numeric
option (including NaN/inf, which a plain argparse `type=float` happily
accepted), file-path options, database options, and `plan`/`galaxy`.

Every way found to get a traceback or hang is pinned as a regression
test with the exact argv that reproduced it; the property tests draw
from everything else, so they guard the rest of the surface.
"""

import argparse
import contextlib
import os
import json
import math
import random
import secrets
import signal
import sys
from unittest import mock

import pytest
from hypothesis import HealthCheck, assume, example, given, note, settings
from hypothesis import strategies as st

import generate
from stellarObjects.starData import STAR_TYPE_PATTERN
from stellarObjects import _db, generationLimits, log, program_constants
from stellarObjects.galaxyDensity import build_galaxy_shape
from stellarObjects.galaxySkeleton import expected_system_count_at_density_1
from stellarObjects.utils import ly_to_pc

from tests import worker_patches
from tests.bughunt_support import mysql_argv
from tests.fuzz_support import hostile_text, scaled

EDGE_PC = ly_to_pc(program_constants.DEFAULT_SECTOR_EDGE_LY)
DB_SETTINGS = settings(suppress_health_check=[HealthCheck.function_scoped_fixture, HealthCheck.too_slow])


# ---------------------------------------------------------------------------
# Harness
# ---------------------------------------------------------------------------

class _Timeout(BaseException):
    """BaseException so no `except Exception` inside the CLI swallows it."""


@contextlib.contextmanager
def _time_limit(seconds):
    def _raise(signum, frame):
        raise _Timeout(f"did not finish within {seconds}s")

    previous = signal.signal(signal.SIGALRM, _raise)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)


@contextlib.contextmanager
def _deterministic(seed):
    """`main()` seeds `random` from `secrets.randbits`, a new galaxy's seed
    comes from `secrets.token_bytes`, and generation draws from `random`
    or the galaxy seed (GEN.39): pin both."""
    rng = random.Random(seed)
    with mock.patch.object(secrets, "randbits", rng.getrandbits), \
            mock.patch.object(secrets, "token_bytes", lambda n: rng.getrandbits(8 * n).to_bytes(n, "big")):
        yield


@pytest.fixture(autouse=True)
def _restore_cli_globals():
    """`build_system_config` writes --flavor-chance-*/--max-planet-flavor
    straight into `program_constants`, and `main()` reconfigures the
    shared logger -- undo both after every test."""
    saved = {name: getattr(program_constants, name)
             for name in ("FLAVOR_CHANCE_SYSTEM", "FLAVOR_CHANCE_PLANET", "MAX_FLAVOR_TOTAL")}
    yield
    for name, value in saved.items():
        setattr(program_constants, name, value)
    log.reset_console()
    log.configure(log.NORMAL)


def _parse(argv):
    """`process_args()` on `argv`: the Namespace, or the SystemExit code."""
    old = sys.argv
    sys.argv = ["generate.py"] + list(argv)
    try:
        return generate.process_args()
    except SystemExit as exc:
        return exc.code
    finally:
        sys.argv = old


def _main(argv, *, seed=0, limit=60):
    """Runs `generate.main()`: returns 0 for a normal return, else the
    SystemExit code. Any other exception (a traceback) or a hang
    propagates, naming the argv."""
    old = sys.argv
    sys.argv = ["generate.py"] + list(argv)
    try:
        with _deterministic(seed), _time_limit(limit):
            generate.main()
        return 0
    except SystemExit as exc:
        return exc.code if isinstance(exc.code, int) else 1
    finally:
        sys.argv = old


def assert_clean(argv, *, seed=0, limit=60, allowed=(0, 1, 2)):
    try:
        code = _main(argv, seed=seed, limit=limit)
    except _Timeout as exc:
        pytest.fail(f"hang: generate.py {' '.join(argv)!s}: {exc}")
    except Exception as exc:  # noqa: BLE001 - a traceback is exactly the failure
        pytest.fail(f"traceback: generate.py {' '.join(map(str, argv))}: {type(exc).__name__}: {exc}")
    assert code in allowed, f"generate.py {' '.join(argv)} exited {code!r}"
    return code


# ---------------------------------------------------------------------------
# Option table, read off the real parser
# ---------------------------------------------------------------------------

_PARSER, _SUBPARSERS = generate.build_parser()
COMMANDS = sorted(_SUBPARSERS)
_SKIP_DESTS = {"help", "version", "mysql_host", "mysql_port", "mysql_user", "mysql_password",
               "mysql_database", "system_file", "output", "debug"}


def _options(command):
    """[(option_string, kind, choices)] for every fuzzable option of a
    subcommand -- `kind` is 'flag' (no value), 'int', 'float', 'str'."""
    table = []
    for action in _SUBPARSERS[command]._actions:
        if not action.option_strings or action.dest in _SKIP_DESTS:
            continue
        if action.nargs == 0:
            kind = "flag"
        elif action.type is int:
            kind = "int"
        elif action.type in (float, generate.finite_float):
            kind = "float"
        else:
            kind = "str"
        for option in action.option_strings:
            table.append((option, kind, action.choices))
    return table


OPTION_TABLE = {command: _options(command) for command in COMMANDS}

INT_TOKENS = ["0", "-0", "1", "2", "-1", "-7", "3", "10", "2147483648", "9" * 30, "-" + "9" * 30]
FINITE_FLOAT_TOKENS = ["0", "-0", "0.0", "1", "0.5", "1e-300", "5e-324", "1e300", "-1", "-1e300", "1e308", "2.5"]
NON_FINITE_TOKENS = ["nan", "NaN", "inf", "-inf", "Infinity", "1e309"]
JUNK_TOKENS = ["", " ", "x", "0x10", "1_000", "--", "-", "+", "--nope", "+nope", "-nope", "=", "1,5", "١٢"]


def _value_tokens(kind, choices, allow_non_finite):
    if choices:
        return st.one_of(st.sampled_from(list(choices)), st.sampled_from(JUNK_TOKENS))
    if kind == "int":
        return st.one_of(st.sampled_from(INT_TOKENS), st.integers(-10**6, 10**6).map(str),
                         st.sampled_from(JUNK_TOKENS + FINITE_FLOAT_TOKENS[:3]))
    if kind == "float":
        pool = FINITE_FLOAT_TOKENS + (NON_FINITE_TOKENS if allow_non_finite else [])
        return st.one_of(st.sampled_from(pool), st.floats(-1e6, 1e6).map(repr), st.sampled_from(JUNK_TOKENS))
    return st.one_of(hostile_text, st.sampled_from(["G2V", "M5V", "O5IA+", "K0VII", "young", "old"]))


@st.composite
def argv_for(draw, command, allow_non_finite=True, junk=True):
    table = OPTION_TABLE[command]
    argv = [command]
    for _ in range(draw(st.integers(0, 8))):
        option, kind, choices = draw(st.sampled_from(table))
        argv.append(option)
        if kind != "flag":
            value = draw(_value_tokens(kind, choices, allow_non_finite))
            if draw(st.integers(0, 5)) == 0 and option.startswith("--"):
                argv[-1] = f"{option}={value}"  # the --opt=value spelling
            else:
                argv.append(value)
        if junk and draw(st.integers(0, 9)) == 0:
            argv.append(draw(st.sampled_from(JUNK_TOKENS)))
    return argv


# ---------------------------------------------------------------------------
# Parse/validate layer: fast, no database
# ---------------------------------------------------------------------------

@settings(max_examples=scaled(150))
@given(data=st.data(), command=st.sampled_from(COMMANDS))
def test_process_args_never_tracebacks(data, command):
    argv = data.draw(argv_for(command))
    note(f"argv={argv}")
    result = _parse(argv)
    if isinstance(result, argparse.Namespace):
        assert result.command == command
    else:
        assert result in (0, 2), f"argv={argv} exited {result!r}"


@settings(max_examples=scaled(60))
@given(tokens=st.lists(st.one_of(hostile_text, st.sampled_from(JUNK_TOKENS + COMMANDS + ["--version", "-h"])),
                       max_size=6))
def test_arbitrary_token_soup_never_tracebacks(tokens):
    result = _parse(tokens)
    assert isinstance(result, argparse.Namespace) or result in (0, 2), (tokens, result)


@settings(max_examples=scaled(120))
@given(data=st.data(), command=st.sampled_from(COMMANDS))
def test_accepted_values_satisfy_every_documented_bound(data, command):
    """Whatever validation lets through (finite tokens only -- see the
    non-finite test below) must satisfy the bounds the validators and
    help text promise."""
    argv = data.draw(argv_for(command, allow_non_finite=False, junk=False))
    args = _parse(argv)
    if not isinstance(args, argparse.Namespace):
        return
    note(f"argv={argv}")
    get = lambda name: getattr(args, name, None)  # noqa: E731
    if get("num_orbits") is not None:
        assert 0 <= args.num_orbits <= generationLimits.MAX_NUM_ORBITS
    for name in ("flavor_chance_system", "flavor_chance_planet"):
        if get(name) is not None:
            assert 0.0 <= get(name) <= 1.0
    if command in ("sector", "galaxy"):
        assert get("density") is None or get("density") > 0
        assert get("num_systems") is None or get("num_systems") >= 1
        assert get("min_habitable") >= 0
        assert get("density") is None or get("num_systems") is None
        if get("num_systems") is not None:
            assert get("min_habitable") <= get("num_systems")
    if command == "sector":
        assert args.num_sectors >= 1
    if command == "galaxy":
        for name in ("ring", "slot", "max_ring"):
            assert get(name) is None or get(name) >= 0
        for name in ("ring", "max_ring"):
            assert get(name) is None or get(name) <= generationLimits.MAX_GENERATE_RING
        assert get("limit") is None or 1 <= get("limit") <= generationLimits.MAX_GENERATE_LIMIT
        assert get("radius_pc") is None or 0 < get("radius_pc") <= generationLimits.MAX_GENERATE_RADIUS_PC
        assert get("min_start_density") is None or get("min_start_density") > 0
    if command == "plan":
        assert 0 <= args.arm_amplitude < 1 and args.max_ring >= 1


# Float options whose validators use `x <= 0` / `x < 0 or x >= 1` style
# checks, which NaN (every comparison False) used to slip straight through
# (now rejected at parse time by `generate.finite_float`).
_FLOAT_OPTIONS = sorted({(c, o) for c in COMMANDS for o, kind, _ in OPTION_TABLE[c] if kind == "float"})
_BASE_ARGV = {"galaxy": ["--ring", "0"]}


@pytest.mark.parametrize("argv", [
    ["galaxy", "--workers=--"],
    ["galaxy", "--flavor-chance-planet=--"],
    ["system", "--flavor-chance-system=--"],
    ["plan", "--disk-scale-length-pc=--"],
])
def test_a_lone_double_dash_value_is_a_usage_error(argv):
    # PERF.27: Python 3.9's argparse stored [] for these instead of an
    # error, and the validators crashed comparing a list with a number.
    assert _parse(argv) == 2


@settings(max_examples=scaled(60))
@given(target=st.sampled_from(_FLOAT_OPTIONS), token=st.sampled_from(NON_FINITE_TOKENS))
@example(target=("sector", "--density"), token="nan")                   # hung in _sample_poisson_count
@example(target=("galaxy", "--density"), token="nan")                   # hung too
@example(target=("galaxy", "--radius-pc"), token="nan")                 # math.ceil(nan) ValueError
@example(target=("galaxy", "--radius-pc"), token="inf")                 # math.ceil(-inf) OverflowError
@example(target=("galaxy", "--min-start-density"), token="nan")         # filter silently off
@example(target=("plan", "--arm-amplitude"), token="nan")               # NaN reached MySQL
@example(target=("plan", "--disk-scale-length-pc"), token="nan")        # NaN reached MySQL
@example(target=("plan", "--calibration-radius-pc"), token="inf")       # math domain error
def test_non_finite_numbers_are_rejected_by_validation(target, token):
    command, option = target
    assert _parse([command] + _BASE_ARGV.get(command, []) + [option, token]) == 2


# ---------------------------------------------------------------------------
# Full runs: system (no database needed with --output)
# ---------------------------------------------------------------------------

SYSTEM_STAR_TYPES = ["G2V", "M9V", "O0V", "B0IA+", "M5VII", "K5VI", "A0III", "O5D", "G0IAB", "F5II", "g2v"]


@st.composite
def system_run_argv(draw):
    argv = ["system", "--quiet"]
    for name, _attr, _desc in generate.TRISTATE_OPTIONS:
        sign = draw(st.sampled_from(["+", "-", None, None]))
        if sign:
            argv.append(f"{sign}{name}")
    if draw(st.booleans()):
        argv += ["--star-type", draw(st.sampled_from(SYSTEM_STAR_TYPES))]
    if draw(st.booleans()):
        argv += ["--num-orbits", draw(st.sampled_from(["0", "1", "-1", "500", "9" * 25]))]
    if draw(st.booleans()):
        argv += ["--name", draw(hostile_text.filter(lambda s: not s.startswith(("-", "+"))))]
    if draw(st.booleans()):
        argv += ["--age", draw(st.sampled_from(["young", "old", "ancient"]))]
    for option in ("--flavor-chance-system", "--flavor-chance-planet"):
        if draw(st.booleans()):
            argv += [option, draw(st.sampled_from(["0", "1", "0.5", "1.0000001", "-0", "5e-324", "nan", "inf"]))]
    if draw(st.booleans()):
        argv.append("--max-planet-flavor")
    if draw(st.booleans()):
        argv.append("--markdown")
    return argv


@settings(max_examples=scaled(40))
@given(argv=system_run_argv(), seed=st.integers(0, 2**32))
def test_system_runs_clean_to_a_file(tmp_path_factory, argv, seed):
    out = tmp_path_factory.mktemp("sys") / "page.txt"
    code = assert_clean(argv + ["--output", str(out)], seed=seed)
    if code == 0:
        text = out.read_text(encoding="utf-8")
        assert text.strip()
        assert "Traceback" not in text


@settings(DB_SETTINGS, max_examples=scaled(10))
@given(argv=system_run_argv(), seed=st.integers(0, 2**32))
def test_system_runs_clean_into_the_database(mysql_config, argv, seed):
    assert_clean(argv + mysql_argv(mysql_config), seed=seed)


def test_system_output_to_stdout(capsys):
    assert assert_clean(["system", "--quiet", "--star-type", "G2V", "--output", "-"]) == 0
    assert capsys.readouterr().out.strip()


# ---- known tracebacks on the system subcommand ----------------------------

@settings(max_examples=scaled(30))
@given(star_type=st.one_of(hostile_text, st.sampled_from(SYSTEM_STAR_TYPES).map(lambda t: t + "x")),
       command=st.sampled_from(["system", "sector"]))
@example(star_type="bogus", command="system")   # was a raw ValueError traceback
@example(star_type="X9V", command="system")
@example(star_type="G2", command="system")
@example(star_type="G2I", command="system")
@example(star_type="G2Vjunk", command="sector")
def test_malformed_star_type_is_a_usage_error(star_type, command):
    assume(star_type and not star_type.startswith(("-", "+"))
           and not STAR_TYPE_PATTERN.fullmatch(star_type.upper()))
    assert _parse([command, "--quiet", "--star-type", star_type]) == 2


@pytest.mark.parametrize("content", [None, "{not json", "42", "null", "[]", '"G2V"', '{"star_type": "G2Vjunk"}',
                                     '{"star_type": 7}', b"\xff\xfe{"])
def test_bad_system_file_is_a_clean_error(tmp_path, content):
    path = tmp_path / "spec.json"
    if isinstance(content, bytes):
        path.write_bytes(content)  # not UTF-8
    elif content is not None:
        path.write_text(content, encoding="utf-8")
    assert_clean(["system", "--quiet", "--system-file", str(path), "--output", "-"], allowed=(1, 2))


def test_unwritable_output_path_is_a_clean_error(tmp_path):
    assert_clean(["system", "--quiet", "--output", str(tmp_path / "missing" / "page.txt")], allowed=(1, 2))


def test_unwritable_debug_file_is_a_clean_error(tmp_path):
    assert_clean(["system", "--debug", str(tmp_path / "missing" / "log.txt"), "--output", "-"],
                 allowed=(1, 2))


@pytest.mark.parametrize("port", ["1", "-1", "0", "99999999"])
def test_unreachable_database_is_a_clean_error(mysql_config, port):
    argv = mysql_argv(mysql_config)
    argv[argv.index("--mysql-port") + 1] = port
    assert_clean(["system", "--quiet"] + argv, allowed=(1, 2))


# ---------------------------------------------------------------------------
# Full runs: sector and phenomenon (database)
# ---------------------------------------------------------------------------

@st.composite
def sector_run_argv(draw):
    argv = ["sector", "--quiet"]
    for name, _attr, _desc in generate.TRISTATE_OPTIONS:
        sign = draw(st.sampled_from(["+", "-", None, None, None]))
        if sign:
            argv.append(f"{sign}{name}")
    count = draw(st.sampled_from(["num", "density", "none"]))
    if count == "num":
        argv += ["--num-systems", draw(st.sampled_from(["0", "1", "2", "3", "-1"]))]
    elif count == "density":
        argv += ["--density", draw(st.sampled_from(["5e-324", "1e-300", "0.01", "0.5", "1", "-1", "0"]))]
    if draw(st.booleans()):
        argv += ["--min-habitable", draw(st.sampled_from(["0", "1", "2", "-1", "99"]))]
    if draw(st.booleans()):
        argv += ["--num-sectors", draw(st.sampled_from(["0", "1", "2", "-3"]))]
    if draw(st.booleans()):
        argv += ["-n", draw(hostile_text.filter(lambda s: not s.startswith(("-", "+"))))]
    if draw(st.booleans()):
        argv += ["--star-type", draw(st.sampled_from(SYSTEM_STAR_TYPES))]
    return argv


@settings(DB_SETTINGS, max_examples=scaled(20))
@given(argv=sector_run_argv(), seed=st.integers(0, 2**32))
def test_sector_runs_clean(mysql_config, argv, seed):
    if "--num-systems" not in argv and "--density" not in argv:
        argv = argv + ["--num-systems", "1"]  # the default of 10 is slow, not interesting
    assert_clean(argv + mysql_argv(mysql_config), seed=seed)


@settings(DB_SETTINGS, max_examples=scaled(3))
@given(density=st.sampled_from(["inf", "1e300", "1e308"]))
def test_sector_with_absurd_density_stops_when_full(mysql_config, density):
    """An infinite/astronomical --density must not try to generate
    ~infinity systems: placement stops once the cube is full."""
    assert_clean(["sector", "--quiet", "--density", density] + mysql_argv(mysql_config), limit=90)


def test_sector_density_nan_does_not_hang(mysql_config):
    assert_clean(["sector", "--quiet", "--density", "nan"] + mysql_argv(mysql_config), limit=5)


class _TooManyConfigs(Exception):
    pass


def test_huge_num_systems_does_not_build_every_config_up_front(mysql_config, monkeypatch):
    real = generate.build_system_config
    calls = []

    def counting(args):
        calls.append(1)
        if len(calls) > 10_000:
            raise _TooManyConfigs(f"{len(calls)} configs built for one sector")
        return real(args)

    monkeypatch.setattr(generate, "build_system_config", counting)
    assert_clean(["sector", "--quiet", "--num-systems", "1000000000"] + mysql_argv(mysql_config), limit=60)


@settings(DB_SETTINGS, max_examples=scaled(15))
@given(
    ptype=st.one_of(st.none(), st.sampled_from(list(program_constants.PHENOMENON_TYPE_CHOICES) + ["wormhole"])),
    anchor=st.booleans(),
    orbits=st.one_of(st.none(), st.sampled_from(["0", "1", "-1", "40", "9" * 20])),
    name=st.one_of(st.none(), hostile_text.filter(lambda s: not s.startswith(("-", "+")))),
    seed=st.integers(0, 2**32),
)
def test_phenomenon_runs_clean(mysql_config, ptype, anchor, orbits, name, seed):
    argv = ["phenomenon", "--quiet"]
    if ptype:
        argv += ["--type", ptype]
    if anchor:
        argv.append("--anchor-system")
    if orbits is not None:
        argv += ["--num-orbits", orbits]
    if name is not None:
        argv += ["--name", name]
    assert_clean(argv + mysql_argv(mysql_config), seed=seed)


@pytest.mark.parametrize("sector_id", ["999999", "-1", "0"])
def test_phenomenon_unknown_sector_id_is_a_clean_error(mysql_config, sector_id):
    assert_clean(["phenomenon", "--quiet", "--type", "nebula", "--sector-id", sector_id]
                 + mysql_argv(mysql_config), allowed=(1, 2))


# ---------------------------------------------------------------------------
# plan
# ---------------------------------------------------------------------------

@st.composite
def plan_argv(draw):
    """A galaxy shape anywhere in the region `build_galaxy_shape` can
    evaluate -- positive finite scales, a non-zero arm count, a pitch
    strictly inside (0, 90] degrees, a calibration radius within a few
    dozen scale lengths (see the degenerate-shape test for what's rejected outside) --
    at scales from toy (1 pc) to absurd (1e5 pc), with a small --max-ring
    so a huge galaxy's scan stays fast."""
    length = draw(st.sampled_from([1.0, 12.0, 40.0, 2800.0, 1e5]))
    argv = [
        # The bright-star scatter walks every cell of the outline; its own
        # tests (test_bright_star_scatter.py) cover it on a small galaxy.
        "plan", "--quiet", "--no-bright-stars",
        "--disk-scale-length-pc", repr(length),
        "--disk-scale-height-pc", repr(draw(st.sampled_from([1e-3, 1.0, 350.0, 1e5]))),
        "--bulge-scale-radius-pc", repr(draw(st.sampled_from([1e-3, 10.0, 200.0, 1e5]))),
        "--bulge-amplitude", repr(draw(st.sampled_from([0.0, 1.0, 1e3, -0.5]))),
        "--arm-count", str(draw(st.sampled_from([1, 2, 4, 1000, -2]))),
        "--pitch-angle-deg", repr(draw(st.sampled_from([1e-3, 15.0, 89.999, 90.0, -15.0]))),
        # 1.0/-0.1 and --max-ring 0 are the validator's own rejections (exit 2),
        # drawn less often than the accepted values.
        "--arm-amplitude", repr(draw(st.sampled_from([0.0, 0.4, 0.999999, 0.0, 0.4, 1.0, -0.1]))),
        "--max-ring", str(draw(st.sampled_from([1, 5, 200, 1, 5, 200, 0]))),
    ]
    if draw(st.booleans()):
        argv += ["--calibration-radius-pc", repr(length * draw(st.sampled_from([1e-3, 1.0, 2.82, 30.0])))]
    return argv


@settings(DB_SETTINGS, max_examples=scaled(25))
@given(argv=plan_argv())
def test_plan_runs_clean_across_galaxy_shapes(mysql_config, argv):
    code = assert_clean(argv + mysql_argv(mysql_config), limit=60)
    if code == 0:
        conn = _db.get_connection(mysql_config)
        try:
            bounds = _db.get_galaxy_bounds(conn)
        finally:
            conn.close()
        assert bounds is not None


_PLAN_TRACEBACKS = [
    (["--arm-count", "0"], "ZeroDivisionError (math.pi / arm_count)"),
    (["--pitch-angle-deg", "0"], "ZeroDivisionError (log(...) / tan(0))"),
    (["--disk-scale-length-pc", "0"], "ZeroDivisionError"),
    (["--disk-scale-length-pc", "-1"], "ValueError: math domain error"),
    (["--disk-scale-height-pc", "0"], "ZeroDivisionError (z / height)"),
    (["--bulge-scale-radius-pc", "0"], "ZeroDivisionError"),
    (["--bulge-scale-radius-pc", "-1"], "OverflowError: math range error"),
    (["--calibration-radius-pc", "0"], "ValueError: math domain error (log(0))"),
    (["--calibration-radius-pc", "1e300"], "ZeroDivisionError (k_norm = 1 / underflowed density)"),
    (["--disk-scale-height-pc", "inf"], "pymysql ProgrammingError: inf can not be used with MySQL"),
]


@pytest.mark.parametrize("extra,consequence", [
    pytest.param(extra, why, id=" ".join(extra)) for extra, why in _PLAN_TRACEBACKS
])
def test_plan_degenerate_shape_is_a_usage_error(extra, consequence):
    """Regression: each of these used to die with the listed raw
    traceback; validate_plan_args only checked --arm-amplitude and
    --max-ring. Now every one is a usage error before touching MySQL."""
    assert _parse(["plan", "--quiet"] + extra) == 2


# ---------------------------------------------------------------------------
# galaxy, against a tiny flat galaxy
# ---------------------------------------------------------------------------

# Relative density ~1 everywhere near the center (no bulge, no arms, huge
# scale lengths, calibrated at 1 pc) so every sector expects the realistic
# ~4 systems -- fast -- and an outline of rings 0..3, layers -1..1.
_FLAT_SHAPE = build_galaxy_shape(
    disk_scale_length_pc=1e6, disk_scale_height_pc=1e6, bulge_scale_radius_pc=1.0,
    bulge_amplitude=0.0, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.0,
    calibration_radius_pc=1.0,
)
_OUTER_RING = 3


def _plan_flat_galaxy(mysql_config):
    _db.save_galaxy_shape(_FLAT_SHAPE, edge_pc=EDGE_PC, outer_ring_index=_OUTER_RING,
                          expected_system_count_at_density_1=expected_system_count_at_density_1(),
                          config=mysql_config)
    _db.replace_galaxy_layers([(1, _OUTER_RING), (0, _OUTER_RING), (-1, _OUTER_RING)], config=mysql_config)


def _placed_sector_id(mysql_config):
    """Generates ring 0 slot 0 and returns its sector id -- a real,
    galaxy-placed sector for --center-sector."""
    assert_clean(["galaxy", "--quiet", "--ring", "0", "--limit", "1", "--num-systems", "1"]
                 + mysql_argv(mysql_config))
    conn = _db.get_connection(mysql_config)
    try:
        return conn.execute("SELECT MIN(id) AS id FROM sectors").fetchone()["id"]
    finally:
        conn.close()


_FLAT_CENTERS = {}  # database name -> the placed center sector's id
GALAXY_INTS = ["0", "1", "2", "3", "4", "-1", "1000000000", "-1000000000"]
GALAXY_FLOATS = ["5e-324", "1e-9", "0.5", "4", "12", "-1", "0", "1e300"]


@st.composite
def galaxy_argv(draw, center_id):
    mode = draw(st.sampled_from(["ring", "slot", "center", "random"]))
    argv = ["galaxy", "--quiet"]
    if mode in ("ring", "slot"):
        argv += ["--ring", draw(st.sampled_from(GALAXY_INTS))]
        if draw(st.booleans()):
            argv += ["--layer", draw(st.sampled_from(GALAXY_INTS))]
        if mode == "slot":
            argv += ["--slot", draw(st.sampled_from(GALAXY_INTS))]
        else:
            argv += ["--limit", draw(st.sampled_from(["1", "2", "0", "-1"]))]
            if draw(st.booleans()):
                argv.append("--yes")
    elif mode == "center":
        argv += ["--center-sector", draw(st.sampled_from([str(center_id), "999999", "-1", "0"])),
                 "--radius-pc", draw(st.sampled_from(["5e-324", "1e-9", "0.5", "4", "-1", "0"]))]
    else:
        argv += ["--radius-pc", draw(st.sampled_from(["1e-9", "0.5", "-1", "0"]))]
        if draw(st.booleans()):
            argv += ["--max-ring", draw(st.sampled_from(["0", "1", "-1", "1000000000000"]))]
        if draw(st.booleans()):
            argv += ["--min-start-density", draw(st.sampled_from(["1e-300", "0.5", "1e300", "0", "-1"]))]
    if mode != "slot" and draw(st.booleans()):
        argv += ["--num-systems", draw(st.sampled_from(["1", "0", "-1"]))]
    elif mode != "slot" and draw(st.booleans()):
        argv += ["--density", draw(st.sampled_from(["1e-300", "0.2", "1", "0", "-1"]))]
    return argv


@settings(DB_SETTINGS, max_examples=scaled(30))
@given(data=st.data(), seed=st.integers(0, 2**32))
def test_galaxy_runs_clean_with_hostile_numbers(mysql_config, data, seed):
    if mysql_config.database not in _FLAT_CENTERS:
        _plan_flat_galaxy(mysql_config)
        _FLAT_CENTERS[mysql_config.database] = _placed_sector_id(mysql_config)
    argv = data.draw(galaxy_argv(_FLAT_CENTERS[mysql_config.database]))
    note(f"argv={argv}")
    assert_clean(argv + mysql_argv(mysql_config), seed=seed, limit=60)
    # Nothing generated may ever land outside the stored outline.
    conn = _db.get_connection(mysql_config)
    try:
        rows = conn.execute("SELECT ring_index, layer_index FROM sectors").fetchall()
    finally:
        conn.close()
    assert all(0 <= r["ring_index"] <= _OUTER_RING and -1 <= r["layer_index"] <= 1 for r in rows), rows


def test_galaxy_before_plan_is_a_clean_refusal(mysql_config):
    assert assert_clean(["galaxy", "--quiet", "--ring", "0", "--limit", "1", "--num-systems", "1"]
                        + mysql_argv(mysql_config)) == 1


_ASKED_ENV = "PLANETGEN_TEST_ASKED_ADDRESSES"


def _recording_fill_task(payload):
    """A `_fill_sector_task` that only records the address it was asked
    for, in the file `$PLANETGEN_TEST_ASKED_ADDRESSES` names."""
    worker_patches.append_json(os.environ[_ASKED_ENV], list(payload["address"]))
    return {"sector_id": 0, "name": "stub", "systems": 0, "phenomena": 0, "summary": ""}


@pytest.mark.parametrize("radius,consequence", [
    ("nan", "was a raw ValueError: cannot convert float NaN to integer"),
    ("inf", "was a raw OverflowError from math.ceil(-inf)"),
    ("1e300", "was a hang enumerating ~1e300 rings; now past the radius bound"),
    ("1e308", "the largest finite radius; now past the radius bound"),
    (repr(generationLimits.MAX_GENERATE_RADIUS_PC), "the largest allowed radius, wider than the test galaxy"),
])
def test_galaxy_absurd_radius_is_clean(mysql_config, monkeypatch, tmp_path, radius, consequence):
    """A radius beyond the whole galaxy just means "every sector in it":
    the enumeration must be trimmed to the galaxy up front, not walk the
    whole sphere. A radius past `generationLimits.MAX_GENERATE_RADIUS_PC`
    is a usage error before anything is generated. Generating the ~150
    sectors of the test galaxy for real would take a while and prove
    nothing extra, so each one is stubbed and only the addresses asked
    for are checked."""
    _plan_flat_galaxy(mysql_config)
    center = _placed_sector_id(mysql_config)
    # The stub is the task itself, so it reaches the workers (PERF.21).
    monkeypatch.setenv(_ASKED_ENV, str(tmp_path / "asked.jsonl"))
    monkeypatch.setattr(generate, "_fill_sector_task", _recording_fill_task)
    monkeypatch.setattr(generate, "_log_saved", lambda *a, **k: None)
    assert_clean(["galaxy", "--quiet", "--center-sector", str(center), "--radius-pc", radius,
                  "--num-systems", "1"] + mysql_argv(mysql_config), limit=10)
    asked = [tuple(address) for address in worker_patches.read_json_lines(str(tmp_path / "asked.jsonl"))]
    cells = sum(generate.ring_sector_count(ring) for ring in range(_OUTER_RING + 1)) * 3
    assert len(asked) == len(set(asked)) <= cells
    assert all(0 <= ring <= _OUTER_RING and -1 <= layer <= 1 for ring, layer, _slot in asked)
    if radius in ("1e300", "1e308"):
        assert asked == []
    elif radius not in ("nan", "inf"):
        assert len(asked) == cells - 1  # every sector but the (already generated) center


def test_galaxy_density_nan_does_not_hang(mysql_config):
    _plan_flat_galaxy(mysql_config)
    assert_clean(["galaxy", "--quiet", "--ring", "0", "--limit", "1", "--density", "nan"]
                 + mysql_argv(mysql_config), limit=5)
