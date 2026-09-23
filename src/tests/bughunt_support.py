# tests/bughunt_support.py

"""
Shared infrastructure for the `test_bughunt_*.py` files -- a deliberate,
exhaustive bug-hunting pass over the generator, the database layer, the
API, and the web interface (see `docs/TODO.md` history and CHANGELOG.md
for the recurring bug categories this targets: overlap/placement bugs,
unit-conversion bugs, NaN/degenerate-physics bugs, infinite loops in name
generation, markdown-rendering breakage).

This is deliberately *not* a new `hypothesis`-based property-testing
layer -- it extends this repo's own existing convention instead:
`stellarObjects/plausibility.py`/`phenomenaPlausibility.py` already split
checks into hard invariants (always enforced) and statistical outliers
(Tukey's fences, reported for human review). Every `test_bughunt_*.py`
file follows the same split, formalized here as two tiers:

- **Tier 1 (hard invariants)** -- plain `assert`/`pytest.raises`. A
  failure here is a genuine bug: an uncaught exception, a NaN/Inf
  escaping a physics calculation, a Hill-sphere overlap, a lost database
  field, a page that 500s. These run unconditionally in every `pytest`
  invocation, same as any other test in this suite.
- **Tier 2 (secondary/statistical)** -- marked `@pytest.mark.tier2`
  (registered in `pytest.ini`). Still runs every time, but checks
  something softer: a statistical outlier, a rare-but-plausible shape, a
  near-boundary condition that isn't itself wrong. A Tier 2 test may still
  use a hard `assert` for genuine invariants discovered along the way --
  the marker is about the *kind* of check a file is built around, not a
  license to skip real bugs.

Since generation draws from the module-level `random` (every
`stellarObjects/*.py` generator uses `import random`, not its own
`random.Random` instance -- confirmed by grep), reproducibility here means
seeding that global RNG before each fuzz iteration, exactly like the
existing ad hoc `random.seed(...)` calls in `test_space_sector.py`/
`test_galaxy_density.py` -- this module just gives that pattern a shared,
documented home and a consistent failure-reporting format.
"""

from __future__ import annotations

import random
import sys
from dataclasses import dataclass, field
from typing import Callable, Sequence

import pytest

# A fixed, reproducible seed list every fuzz test iterates over -- small
# enough to stay fast (a few hundred iterations, matching the existing
# plausibility CLIs' own `--n 150`/`--n 200` defaults), but large enough to
# turn up rare edge cases across many runs. 0/1/-1-adjacent values are
# included deliberately since off-by-one seed-dependent branches are a
# classic source of "works on my seed" bugs.
FUZZ_SEEDS: tuple[int, ...] = tuple(range(0, 150)) + (
    2**31 - 1,
    2**32 - 1,
    999_999_937,  # a large prime, away from any small-seed correlation
)


@dataclass
class FuzzFailure(Exception):
    """
    Raised (or attached to a re-raised exception's `__cause__` chain via
    `report_seed_failure`) so a fuzz-loop failure always names the exact
    seed and call that produced it -- without this, a seeded-loop
    assertion failure just reports "assert False" with no way to
    reproduce it.
    """

    seed: int
    detail: str = ""

    def __str__(self) -> str:  # pragma: no cover - trivial
        return f"failed at seed={self.seed}{': ' + self.detail if self.detail else ''}"


def run_seeded(fn: Callable[[int], None], seeds: Sequence[int] = FUZZ_SEEDS) -> None:
    """
    Calls `fn(seed)` once per seed in `seeds`, seeding the module-level
    `random` RNG first so `fn` (and anything it calls) draws a
    reproducible-but-varied sequence. On the first failure, re-raises with
    the failing seed attached so a human can reproduce it with
    `random.seed(<seed>); fn(<seed>)` alone -- matching the
    `physical_plausibility_cli.py`/`phenomena_plausibility_cli.py`
    convention of reporting exactly what to re-run, not just that
    something failed.

    Args:
        fn: Called once per seed; receives the seed itself as its only
            argument (most callers ignore it -- it's there for logging).
        seeds: The seed list to iterate. Defaults to `FUZZ_SEEDS`.
    """
    for seed in seeds:
        random.seed(seed)
        try:
            fn(seed)
        except Exception as exc:  # noqa: BLE001 - re-raising with context is the point
            raise FuzzFailure(seed=seed, detail=f"{type(exc).__name__}: {exc}") from exc


def run_cli(subcommand: str, argv: Sequence[str]) -> None:
    """
    Runs `generate.py`'s real `main()` for the given subcommand
    (`system`/`sector`/`galaxy`/`plan`/`phenomenon`) via `sys.argv`,
    restoring it afterward -- generalizes `test_galaxy_gen.py`'s own
    `_run_cli`/`_run_sector_gen_cli` helpers (previously duplicated
    per-subcommand, `galaxy`-only) to every subcommand, since the fuzz
    tests here need all five. A clean, expected `parser.error()` rejection
    surfaces as `SystemExit` (argparse's own behavior) -- callers that
    want to fuzz *valid* argv should let `SystemExit` propagate and assert
    on it explicitly; callers checking "no crash on this exact argv"
    should catch it themselves if a `SystemExit` is an acceptable outcome
    for that argv.

    Args:
        subcommand: One of generate.py's subcommands.
        argv: The rest of the command line (excluding "generate.py" and
            the subcommand itself).
    """
    import generate as _generate

    old_argv = sys.argv
    try:
        sys.argv = ["generate.py", subcommand] + list(argv)
        _generate.main()
    finally:
        sys.argv = old_argv


def mysql_argv(mysql_config) -> list[str]:
    """The `--mysql-*` argv fragment pointing a CLI invocation at a
    `mysql_config` fixture's throwaway database -- same helper as
    `test_galaxy_gen.py`'s private `_mysql_argv`, shared here so every
    bug-hunt file that drives the CLI against a real database doesn't
    reimplement it."""
    return [
        "--mysql-host", mysql_config.host,
        "--mysql-port", str(mysql_config.port),
        "--mysql-user", mysql_config.user,
        "--mysql-password", mysql_config.password,
        "--mysql-database", mysql_config.database,
    ]


@dataclass
class Tier2Report:
    """
    Accumulates soft findings for a Tier 2 test so they're reported
    together at the end (via `pytest.fail` if `strict` is requested, or
    just printed to stdout under `-s` otherwise) instead of stopping at
    the first one -- mirrors `plausibility.format_report`'s "collect
    everything, report once" shape.
    """

    findings: list[str] = field(default_factory=list)

    def add(self, message: str) -> None:
        self.findings.append(message)

    def flush(self, label: str) -> None:
        """Prints every accumulated finding, prefixed with a Tier 2 banner
        so it's easy to grep out of a full-suite run's output."""
        if not self.findings:
            return
        print(f"\n[tier2:{label}] {len(self.findings)} soft finding(s):")
        for msg in self.findings:
            print(f"  - {msg}")


def tier2(fn):
    """Convenience decorator combining `@pytest.mark.tier2` with nothing
    else -- kept as a function (not just the raw marker) so call sites
    read `@bughunt_support.tier2` next to the Tier 1 tests in the same
    file without an extra `import pytest` just for the marker."""
    return pytest.mark.tier2(fn)
