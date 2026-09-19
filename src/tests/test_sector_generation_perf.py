# tests/test_sector_generation_perf.py

"""
Sector generation performance benchmark
=========================================

Not a correctness test -- profiles how long `generate.generate_sector`
(the `sector` subcommand's actual, DB-free, in-memory generation path:
`build_sector_configs` + one `StarSystem` per config + `SpaceSector`
placement + `generate_sector_phenomena`) takes, broken down by phase.

The phase breakdown comes from `--debug`'s own logging, not a separate
profiling API: `stellarObjects.log.timed_phase` (used by `generate_sector`
and `StarSystem.__init__` at each of their major phases) logs
`"{label}: {elapsed}ms"` at DEBUG severity, with a timestamp, exactly as
`generate.py sector --debug` would print to a terminal. This module just
attaches its own `logging.Handler` to capture those same records instead
of a human reading them off the console, and aggregates the elapsed
values already embedded in each message. A `cProfile` run cross-checks
the same conclusion at the function level, in case a phase boundary
placed by hand hides an internal hotspot.

Skipped by default (opt-in via `PLANETGEN_RUN_PERF_BENCHMARK=1`), the
same "opt-in real work" treatment `conftest.py`'s `wikijs_config`/
`mediawiki_config` give a real external service: this deliberately
generates dozens of sectors purely to gather timing statistics, which
would slow down every ordinary `pytest` run for no benefit to a
correctness check. No database is involved -- `generate_sector` itself
never touches one (see its own docstring); only `run_sector`'s later
`_db.save_sector` call would.

Run with:
    PLANETGEN_RUN_PERF_BENCHMARK=1 pytest -s src/tests/test_sector_generation_perf.py

`-s` is required to actually see the printed report; pytest otherwise
captures stdout and only shows it on failure.
"""

import cProfile
import io
import logging
import os
import pstats
import re
import time
from collections import defaultdict

import pytest

import generate as gen
from stellarObjects import log

pytestmark = pytest.mark.skipif(
    not os.environ.get("PLANETGEN_RUN_PERF_BENCHMARK"),
    reason="opt-in perf benchmark -- set PLANETGEN_RUN_PERF_BENCHMARK=1 to run it "
           "(pytest -s src/tests/test_sector_generation_perf.py, add -s to see the report)",
)

# Matches exactly what `log.timed_phase` emits: "<label>: <elapsed>ms".
_PHASE_LINE = re.compile(r"^(?P<label>.+): (?P<ms>[\d.]+)ms$")
# Collapses per-system/per-attempt labels ("generate system 7/35",
# "planet generation attempt 2") down to one bucket per phase so the
# report aggregates across every system/attempt instead of listing each
# one as its own row.
_FRACTION = re.compile(r"\d+/\d+")
_ATTEMPT = re.compile(r"attempt \d+")


def _normalize_label(label):
    label = _FRACTION.sub("N/M", label)
    label = _ATTEMPT.sub("attempt N", label)
    return label


class _PhaseCapture(logging.Handler):
    """
    Collects every `log.timed_phase` DEBUG record emitted while attached
    to the `"planetgen"` logger, bucketed by normalized label.

    Attributes:
        samples (dict): Normalized phase label -> list of elapsed-ms
                        values, one per `timed_phase` call.
    """

    def __init__(self):
        super().__init__(level=logging.DEBUG)
        self.samples = defaultdict(list)

    def emit(self, record):
        match = _PHASE_LINE.match(record.getMessage())
        if match:
            self.samples[_normalize_label(match.group("label"))].append(float(match.group("ms")))


def _sector_args(num_systems, name):
    """Builds a real `sector` subcommand namespace via `generate.py`'s own
    parser/validators (not a hand-built `SimpleNamespace`), so this
    benchmark exercises exactly the option surface/defaults a real
    `generate.py sector` invocation would."""
    parser, command_parsers = gen.build_parser()
    args = parser.parse_args(["sector", "--num-systems", str(num_systems), "--name", name])
    command_parser = command_parsers["sector"]
    gen.validate_shared_generation_args(args, command_parser)
    gen.validate_sector_args(args, command_parser)
    return args


def _print_phase_report(samples, grand_total_ms, title):
    rows = sorted(
        ((label, sum(vals), len(vals), sum(vals) / len(vals), max(vals)) for label, vals in samples.items()),
        key=lambda row: -row[1],
    )
    print(f"\n{title}")
    print(f"{'phase':40s} {'total_ms':>10s} {'count':>7s} {'avg_ms':>9s} {'max_ms':>9s} {'% of total':>11s}")
    for label, total, count, avg, worst in rows:
        pct = 100 * total / grand_total_ms if grand_total_ms else 0.0
        print(f"{label:40s} {total:10.2f} {count:7d} {avg:9.3f} {worst:9.3f} {pct:10.2f}%")


def test_sector_generation_phase_breakdown():
    """
    Generates a handful of sectors under `--debug`-equivalent logging,
    captures every `log.timed_phase` record, and prints (and sanity-checks)
    a per-phase timing report -- which phase of `generate_sector` (and, for
    each system it builds, which phase of `StarSystem.__init__`) the time
    actually goes to. See this module's docstring for how to view the
    report (`-s`).
    """
    sector_sizes = [12, 20, 28]

    log.configure(log.DEBUG)
    # The console handler would otherwise also print every debug line
    # (correct behavior for a real --debug run, but pure noise for a
    # benchmark whose only interesting output is the aggregated report
    # below) -- raise its own threshold above CRITICAL rather than
    # touching the logger's level, so records still reach _capture.
    previous_console_level = log._console_handler.level
    log._console_handler.setLevel(logging.CRITICAL + 1)
    logger = logging.getLogger("planetgen")
    capture = _PhaseCapture()
    logger.addHandler(capture)

    try:
        wall_times = []
        for i, size in enumerate(sector_sizes):
            args = _sector_args(size, f"BenchSector{i}")
            start = time.perf_counter()
            _name, sector = gen.generate_sector(args)
            elapsed_ms = (time.perf_counter() - start) * 1000
            wall_times.append((size, len(sector.entries), elapsed_ms))
    finally:
        logger.removeHandler(capture)
        log._console_handler.setLevel(previous_console_level)
        log.configure(log.NORMAL)

    grand_total_ms = sum(elapsed for _size, _placed, elapsed in wall_times)

    print("\nper-sector wall time (generate_sector, no DB):")
    for size, placed, elapsed in wall_times:
        print(f"  requested {size:3d} systems, placed {placed:3d}: {elapsed:8.2f}ms")

    _print_phase_report(capture.samples, grand_total_ms, "aggregated phase breakdown (all sectors combined):")

    assert grand_total_ms > 0
    assert capture.samples, "expected at least one log.timed_phase record to have been captured"
    # generate_sector always wraps each system's construction in its own
    # "generate system N/M" phase -- if that bucket is missing, the
    # capture/normalization above is broken, not that generation was
    # somehow instantaneous.
    assert "generate system N/M" in capture.samples


def test_sector_generation_cprofile_hotspots():
    """
    Cross-checks the phase breakdown above at the function level: profiles
    one `generate_sector` call with `cProfile` and prints the hottest
    functions by both cumulative and internal (self) time. Catches a
    hotspot hiding *inside* a phase `log.timed_phase` treats as one
    opaque block (e.g. one specific helper deep inside "planet generation
    attempt" actually dominating that phase), which the timestamp-based
    report above can't distinguish on its own.
    """
    args = _sector_args(30, "ProfiledSector")

    profiler = cProfile.Profile()
    profiler.enable()
    _name, sector = gen.generate_sector(args)
    profiler.disable()

    assert len(sector.entries) > 0

    for sort_key, title in (("cumulative", "top 25 by cumulative time"), ("tottime", "top 20 by internal (self) time")):
        buffer = io.StringIO()
        stats = pstats.Stats(profiler, stream=buffer).sort_stats(sort_key)
        stats.print_stats(25 if sort_key == "cumulative" else 20)
        print(f"\ncProfile {title}:\n{buffer.getvalue()}")
