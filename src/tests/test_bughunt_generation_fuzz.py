# tests/test_bughunt_generation_fuzz.py

"""
Tier 1 bug-hunt coverage: `generate.py system`/`sector`/`phenomenon` driven
through many randomized, potentially-conflicting flag combinations (via
the real CLI entry point, same convention as `test_bughunt_cli_edges.py`)
-- unlike that file (one specific documented-incompatible combo per test),
this one draws combinations at random every seed, covering corners no
single hand-picked test enumerates. A `SystemExit` from `parser.error`/
`build_system_config`'s own validation is an acceptable outcome (an
invalid combination being rejected is correct behavior); anything else
escaping (a raw `AttributeError`/`KeyError`/`ZeroDivisionError`/etc. deep
in generation) is a real bug.
"""

import random

import pytest

from stellarObjects import program_constants
from tests.bughunt_support import FUZZ_SEEDS, mysql_argv, run_cli

_SEEDS = FUZZ_SEEDS[:60]

_TRISTATE_FLAGS = [
    "habitable_world", "asteroid_belt", "comets", "large_star", "moons",
    "max_planets", "intelligent_life", "binary_system", "wide_binary", "planets",
]
_STAR_TYPES = [None, "M5V", "K2V", "G2V", "F5V", "A1V", "B3V", "O5V", "G2III", "M5I"]
_AGES = [None, "young", "old"]


def _random_system_argv():
    argv = []
    for flag in _TRISTATE_FLAGS:
        choice = random.choice(["+", "-", None])
        if choice:
            argv.append(f"{choice}{flag}")
    star_type = random.choice(_STAR_TYPES)
    if star_type:
        argv += ["--star-type", star_type]
    age = random.choice(_AGES)
    if age:
        argv += ["--age", age]
    if random.random() < 0.3:
        argv += ["--num-orbits", str(random.choice([0, 1, 3, 8, 20]))]
    if random.random() < 0.3:
        argv += ["--flavor-chance-system", str(round(random.uniform(0.0, 1.0), 2))]
    if random.random() < 0.3:
        argv += ["--flavor-chance-planet", str(round(random.uniform(0.0, 1.0), 2))]
    return argv


def test_random_system_flag_combinations_never_crash_uncleanly(mysql_config):
    def check(seed):
        argv = _random_system_argv()
        try:
            run_cli("system", argv + mysql_argv(mysql_config))
        except SystemExit:
            pass  # a clean, validated rejection is an acceptable outcome
        except Exception as exc:  # noqa: BLE001
            pytest.fail(f"seed={seed} argv={argv}: uncaught {type(exc).__name__}: {exc}")

    from tests.bughunt_support import run_seeded
    run_seeded(check, seeds=_SEEDS)


def test_random_sector_flag_combinations_never_crash_uncleanly(mysql_config):
    def check(seed):
        argv = _random_system_argv()  # sector shares the same shared-option surface
        argv += ["--num-systems", str(random.choice([1, 2, 5]))]
        try:
            run_cli("sector", argv + mysql_argv(mysql_config))
        except SystemExit:
            pass
        except Exception as exc:  # noqa: BLE001
            pytest.fail(f"seed={seed} argv={argv}: uncaught {type(exc).__name__}: {exc}")

    from tests.bughunt_support import run_seeded
    run_seeded(check, seeds=_SEEDS[:30])


_PHENOMENON_TYPES = [
    "black-hole", "neutron-star", "nebula", "supernova-remnant",
    "rogue-planet", "comet", "asteroid-field",
]


def test_every_phenomenon_type_generates_without_crashing(mysql_config):
    def check(seed):
        ptype = random.choice(_PHENOMENON_TYPES)
        argv = ["--type", ptype]
        try:
            run_cli("phenomenon", argv + mysql_argv(mysql_config))
        except SystemExit:
            pass
        except Exception as exc:  # noqa: BLE001
            pytest.fail(f"seed={seed} type={ptype}: uncaught {type(exc).__name__}: {exc}")

    from tests.bughunt_support import run_seeded
    run_seeded(check, seeds=_SEEDS[:20])
