# tests/audit/prevalence_disk_age.py

"""
The disk-and-planet census by star age (GEN.1553, report reports/young-stars-disks-no-bodies-2026-11-11.md): for stars of a chosen mass whose age is pinned in each of
the bins 0 to 1, 1 to 3, 3 to 6, 6 to 10 and 10 to 30 Myr, how many have a disk, bodies, both, or neither, set beside the survey curve a young star's disk follows. It measures; it asserts nothing. It is the check
GEN.1550 to GEN.1553 are built against (a disk wholly gone under 3 Myr and no bodies until 10 Myr is what it was made to find).

The star is pinned by building the system with a `Star` that takes the mass and age given (`initial_mass_sol`, `age_gy`, the same keywords `StarSystem(primary_star_params=...)` passes), so no module is patched
and the age is the draw's input, not its outcome; the record keeps both the age asked for and the age the finished star holds (the age floor, GEN.1546, may raise it until it is removed). The disk is read twice:
`disk_at_birth` straight after the star is made (`Star._finish_init`), and `disk_final` on the finished system (`StarSystem._set_debris_disk` clears a young disk and keeps only a debris one).

    PYTHONPATH=src python3 -m tests.audit.prevalence_disk_age --count 40 --out /tmp/disk-age-census.md
"""

import argparse
import json
import math
import multiprocessing
import os
import signal
import unittest.mock as mock

from tests.audit import prevalence as census

SYSTEM_TIME_LIMIT_S = 240
"""int: Seconds one system may take before it is abandoned and counted."""

MASSES_SOL = (("M2V", 0.4), ("K5V", 0.7), ("G2V", 1.0))
"""tuple: The census stars and the initial mass each is pinned to, solar masses."""

AGE_BINS_MYR = ((0.0, 1.0), (1.0, 3.0), (3.0, 6.0), (6.0, 10.0), (10.0, 30.0))
"""tuple: The age bins, Myr, a star's age is drawn uniform inside."""

SURVEY_DISK_FRACTION = {(0.0, 1.0): 0.85, (1.0, 3.0): 0.65, (3.0, 6.0): 0.30, (6.0, 10.0): 0.05, (10.0, 30.0): 0.0}
"""dict: The share of stars of that age that still hold a dusty disk (Haisch, Lada and Lada 2001 style curve: about 0.8 at 1 Myr, 0.5 at 3, 0.1 at 5 to 6, none past 10), taken at each bin's middle. Stand-in,
Confident Unverified, L0 (GEN.1551 registers the curve when it is built); it is the line to compare against, not a target."""


class _Abandoned(BaseException):
    """Raised by the alarm when one system ran past `SYSTEM_TIME_LIMIT_S`."""


def _pinned_star(captured, mass_sol, age_gy):
    """A `Star` subclass of that mass and age; it copies its disk fields into `captured` the moment it is built."""
    from planetgen.generation.star import Star

    class Pinned(Star):
        def __init__(self, system_config, **kwargs):
            super().__init__(system_config, initial_mass_sol=mass_sol, age_gy=age_gy, **kwargs)
            captured.update(birth_age_gy=self.age, disk_at_birth=getattr(self, "disk_kind", None), yerkes_at_birth=self.yerkes_class)
    return Pinned


def _record(system, captured, asked_myr):
    """One system's census row: the age asked for and the age held, the disk at birth and at the end, and the bodies it ended with."""
    star = system.primary_star
    planets = [body for body in system.all_planets if body.body_type != "a"]
    belts = [body for body in system.all_planets if body.body_type == "a"]
    return {"asked_age_myr": asked_myr, "star_age_myr": star.age * 1.0e3, "yerkes": star.yerkes_class, "yerkes_at_birth": captured.get("yerkes_at_birth"),
            "disk_at_birth": captured.get("disk_at_birth"), "disk_final": getattr(star, "disk_kind", None),
            "planets": len(planets), "belts": len(belts), "bodies": len(planets) + len(belts)}


def draw_cell(label, mass_sol, age_bin, count):
    """`count` systems of one star mass with the age drawn uniform in `age_bin` (Myr), as rows; one over the time limit is `{"abandoned": True}`, a crash `{"error": ...}`."""
    from planetgen.generation import system as system_module
    from planetgen.generation.config import SystemConfig
    from planetgen.util import draw

    stream = draw.Stream(f"disk-age-census:{label}:{age_bin[0]}:{age_bin[1]}")
    rows = []

    def _alarm(*_):
        raise _Abandoned()

    previous = signal.signal(signal.SIGALRM, _alarm)
    try:
        for index in range(count):
            asked_myr = stream.uniform(*age_bin)
            captured = {}
            config = SystemConfig()
            config.PLANETS = True
            config.COMETS = False
            config.BINARY_SYSTEM = False
            signal.alarm(SYSTEM_TIME_LIMIT_S)
            try:
                with mock.patch.object(system_module, "Star", _pinned_star(captured, mass_sol, asked_myr / 1.0e3)):
                    with draw.bound(census.seed_for("disk-age", label, age_bin[0], index)):
                        rows.append(_record(system_module.StarSystem(system_config=config), captured, asked_myr))
            except _Abandoned:
                rows.append({"abandoned": True, "asked_age_myr": asked_myr})
            except Exception as error:
                rows.append({"error": f"{type(error).__name__}: {error}", "asked_age_myr": asked_myr})
            finally:
                signal.alarm(0)
    finally:
        signal.signal(signal.SIGALRM, previous)
    return [{"star": label, "bin": list(age_bin), **row} for row in rows]


def _job(job):
    return draw_cell(*job)


def draw_cells(count, masses=MASSES_SOL, bins=AGE_BINS_MYR, workers=None):
    jobs = [(label, mass, age_bin, count) for label, mass in masses for age_bin in bins]
    with multiprocessing.get_context("fork").Pool(min(workers or os.cpu_count() or 1, len(jobs))) as pool:
        results = pool.map(_job, jobs, chunksize=1)
    return [row for result in results for row in result]


def summarize(rows):
    """`{(star, bin): {n, abandoned, errors, disk_at_birth, disk_final, bodies, both, neither, moved}}`, the shares being fractions of the finished systems."""
    out = {}
    for key in sorted({(row["star"], tuple(row["bin"])) for row in rows}):
        mine = [row for row in rows if (row["star"], tuple(row["bin"])) == key]
        good = [row for row in mine if "bodies" in row]
        n = len(good) or 1
        share = lambda test: sum(1 for row in good if test(row)) / n   # noqa: E731
        out[key] = {"n": len(good), "abandoned": sum(1 for row in mine if row.get("abandoned")), "errors": [row["error"] for row in mine if "error" in row][:3],
                    "disk_at_birth": share(lambda r: r["disk_at_birth"] in ("class-0", "class-i", "class-ii")),
                    "disk_final": share(lambda r: r["disk_final"] in ("class-0", "class-i", "class-ii")),
                    "bodies": share(lambda r: r["bodies"] > 0),
                    "both": share(lambda r: r["bodies"] > 0 and r["disk_at_birth"] is not None),
                    "neither": share(lambda r: r["bodies"] == 0 and r["disk_at_birth"] is None),
                    "age_raised": share(lambda r: r["star_age_myr"] > 1.05 * r["asked_age_myr"] + 0.01),
                    "mean_planets": sum(r["planets"] for r in good) / n}
    return out


def format_report(summary, count):
    lines = [f"# Disk and planet census by star age ({count} systems per cell requested)", "",
             "Shares are fractions of the finished systems. `disk at birth` is the star's disk straight after it is made; `disk final` is what the finished system holds (a debris disk only, today). "
             "`survey` is the dusty-disk share of real stars of that age (stand-in curve); `age raised` is the share whose star age ended above the age asked for (the age floor, GEN.1546).", "",
             "| star | age Myr | n | disk at birth | survey | disk final | with bodies | disk and bodies | neither | age raised | mean planets |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for (star, age_bin), v in summary.items():
        lines.append(f"| {star} | {age_bin[0]:g} to {age_bin[1]:g} | {v['n']} ({v['abandoned']} abandoned) | {v['disk_at_birth']:.2f} | {SURVEY_DISK_FRACTION.get(age_bin, math.nan):.2f} | "
                     f"{v['disk_final']:.2f} | {v['bodies']:.2f} | {v['both']:.2f} | {v['neither']:.2f} | {v['age_raised']:.2f} | {v['mean_planets']:.2f} |")
        if v["errors"]:
            lines.append(f"|  | errors: {v['errors']} | | | | | | | | | |")
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--count", type=int, default=40, help="systems per star mass and age bin")
    parser.add_argument("--out", default="/tmp/disk-age-census.md")
    parser.add_argument("--json", default=None, help="also write the raw rows here")
    parser.add_argument("--workers", type=int, default=None)
    args = parser.parse_args(argv)
    rows = draw_cells(args.count, workers=args.workers)
    if args.json:
        with open(args.json, "w", encoding="utf-8") as handle:
            json.dump(rows, handle)
    report = format_report(summarize(rows), args.count)
    with open(args.out, "w", encoding="utf-8") as handle:
        handle.write(report)
    print(report)


if __name__ == "__main__":
    main()
