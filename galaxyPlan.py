#!/usr/bin/env python
# galaxyPlan.py

"""
Galaxy-Wide Density Skeleton Builder
=======================================

Builds and persists the galaxy's density "skeleton" -- see
`docs/design/galaxy-coordinate-system.md` section 9's storage-analysis
addendum and `stellarObjects/galaxySkeleton.py`'s own module docstring for
the reasoning this implements. In short: a sector's position, density, and
vertices are pure deterministic functions of its `(shell_index,
shell_slot_index)` address and a handful of galaxy-wide shape parameters,
so none of that gets stored per sector here -- what this script actually
persists is:

- `galaxy_shape`: the galaxy's shape parameters, calibration constant,
  sector edge length, and outer edge (the last shell with any qualifying
  content) -- one singleton row.
- `galaxy_shell_band`: one row per contiguous *candidate* slot-index band
  per shell (almost always exactly one) -- a safe, cheap-to-compute
  superset of where a shell's qualifying sectors could be, not an exact
  per-sector list. See `galaxySkeleton.find_shell_bands`'s own docstring
  for exactly what "candidate" means and why that's the right thing to
  store.

No sector content, vertices, or even individual sector addresses are
generated or stored by this script -- that stays lazy, happening only when
a sector is actually visited (`galaxyGen.ensure_sector_generated`), per
this project's own "most sectors stay unvisited forever" design decision.
Re-running this script (e.g. after retuning the shape parameters) always
replaces the entire skeleton wholesale -- there is no incremental update,
since a full build is already well under a minute even at real Milky-Way
scale (see the module docstring above for why: per-shell work here is a
closed-form calculation, not a per-sector scan).

**Parallelization**: shells are scanned outward in order (needed to detect
the galaxy's true edge -- a run of consecutive empty shells), but each
shell's own band-finding is fully independent of every other shell's, so
shells are dispatched to a `multiprocessing.Pool` in fixed-size chunks
(`--chunk-size`, default scaled to `--workers`) -- one process-pool round
trip per chunk, main process only ever reading pool results and writing to
the database (workers never touch it -- there's no need for more than one
writer given how little data this script's own writes involve, a
singleton row plus a few thousand band rows at most).
"""

import argparse
import math
import multiprocessing
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))

from stellarObjects import _db, program_constants
from stellarObjects._version import VersionAction, version_banner
from stellarObjects.galaxyDensity import build_galaxy_shape
from stellarObjects.galaxySkeleton import expected_system_count_at_density_1, find_shell_bands
from stellarObjects.utils import ly_to_pc

DEFAULT_EMPTY_STREAK_TO_STOP = 50
"""int: How many consecutive empty shells confirms the galaxy's true edge
has been reached -- see `galaxySkeleton`'s own docstring on why the
qualifying region is expected to shrink monotonically outward, plus a
safety margin against a razor-thin band a coarse scan could miss for one
single shell."""

DEFAULT_MAX_SHELL = 100000
"""int: Hard cap on how many shells this script will ever scan, regardless
of `--empty-streak-to-stop` -- protects against a pathological shape
parameter choice (e.g. a bulge_amplitude/threshold combination with no
real outward decay) turning into an unbounded scan; real Milky-Way-scale
parameters reach their edge by shell ~4,100, so this cap is never expected
to bind in ordinary use."""


def _worker_find_bands(shape, edge_pc, shell_index, threshold_rho):
    """
    The one unit of work dispatched to each pool worker -- must be a
    plain, picklable module-level function (not a closure/lambda) for
    `multiprocessing.Pool` to use it under the `spawn` start method.

    Args:
        shape (galaxyDensity.GalaxyShape): The galaxy's shape parameters.
        edge_pc (float): The sector edge length, parsecs.
        shell_index (int): The shell index to find bands for.
        threshold_rho (float): The qualifying `relative_density` threshold.

    Returns:
        tuple: `(shell_index, [galaxySkeleton.ShellBand, ...])`.
    """
    return shell_index, find_shell_bands(shape, edge_pc, shell_index, threshold_rho)


def process_args():
    """
    Parses command-line arguments for `galaxyPlan.py`.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Galaxy Density Skeleton Builder",
        epilog="Builds and persists the galaxy_shape/galaxy_shell_band skeleton this galaxy's shape "
               "implies -- see this script's module docstring for what that means and why it's a small, "
               "fast build even at real Milky-Way scale.",
    )

    parser.add_argument('--version', action=VersionAction, banner=version_banner('galaxyPlan.py'))

    shape_group = parser.add_argument_group("galaxy shape (galaxyDensity.GalaxyShape)")
    shape_group.add_argument('--disk-scale-length-pc', type=float, default=2800.0,
                             help="Exponential disk radial scale length, parsecs. Default: 2800 "
                                  "(real Milky Way scale).")
    shape_group.add_argument('--disk-scale-height-pc', type=float, default=350.0,
                             help="Disk vertical scale height, parsecs. Default: 350.")
    shape_group.add_argument('--bulge-scale-radius-pc', type=float, default=200.0,
                             help="Bulge exponential scale radius, parsecs. Default: 200.")
    shape_group.add_argument('--bulge-amplitude', type=float, default=1.0,
                             help="Bulge amplitude, relative to the disk term. Default: 1.0.")
    shape_group.add_argument('--arm-count', type=int, default=2,
                             help="Number of spiral arms. Default: 2 (grand-design).")
    shape_group.add_argument('--pitch-angle-deg', type=float, default=15.0,
                             help="Spiral arm pitch angle, degrees. Default: 15.")
    shape_group.add_argument('--arm-amplitude', type=float, default=0.4,
                             help="Arm/inter-arm density contrast amplitude, in [0, 1). Default: 0.4.")
    shape_group.add_argument('--calibration-radius-pc', type=float, default=None,
                             help="In-plane radius the relative_density=1.0 calibration point sits at. "
                                  "Defaults to build_galaxy_shape's own default (2.82x disk scale length).")

    parser.add_argument('--edge-ly', type=float, default=program_constants.DEFAULT_SECTOR_EDGE_LY,
                        help=f"Sector edge length, light-years. Default: "
                             f"{program_constants.DEFAULT_SECTOR_EDGE_LY}.")
    parser.add_argument('--empty-streak-to-stop', type=int, default=DEFAULT_EMPTY_STREAK_TO_STOP,
                        help=f"Consecutive empty shells before concluding the galaxy's edge has been "
                             f"reached. Default: {DEFAULT_EMPTY_STREAK_TO_STOP}.")
    parser.add_argument('--max-shell', type=int, default=DEFAULT_MAX_SHELL,
                        help=f"Hard cap on shells scanned, regardless of --empty-streak-to-stop. "
                             f"Default: {DEFAULT_MAX_SHELL}.")
    parser.add_argument('--workers', type=int, default=None,
                        help="Worker process count for the parallel per-shell scan. Defaults to "
                             "os.cpu_count() (or 1 if that can't be determined).")
    parser.add_argument('--chunk-size', type=int, default=None,
                        help="Shells dispatched per parallel round. Defaults to 8x --workers.")
    _db.add_mysql_connection_args(parser)

    args = parser.parse_args()

    if args.arm_amplitude < 0 or args.arm_amplitude >= 1:
        parser.error("--arm-amplitude must be in [0, 1).")
    if args.empty_streak_to_stop < 1:
        parser.error("--empty-streak-to-stop must be a positive integer.")
    if args.max_shell < 1:
        parser.error("--max-shell must be a positive integer.")
    if args.workers is not None and args.workers < 1:
        parser.error("--workers must be a positive integer.")
    if args.chunk_size is not None and args.chunk_size < 1:
        parser.error("--chunk-size must be a positive integer.")

    if args.workers is None:
        args.workers = os.cpu_count() or 1
    if args.chunk_size is None:
        args.chunk_size = args.workers * 8

    return args


def build_skeleton(args):
    """
    Runs the full parallel skeleton build and persists it.

    Args:
        args (argparse.Namespace): Parsed arguments (see `process_args`).

    Returns:
        dict: Summary stats -- `shells_scanned`, `outer_shell_index`,
              `total_bands`, `total_candidate_slots`, `elapsed_s`.
    """
    edge_pc = ly_to_pc(args.edge_ly)
    e_value = expected_system_count_at_density_1(args.edge_ly)
    threshold_rho = 1.0 / e_value

    shape = build_galaxy_shape(
        disk_scale_length_pc=args.disk_scale_length_pc,
        disk_scale_height_pc=args.disk_scale_height_pc,
        bulge_scale_radius_pc=args.bulge_scale_radius_pc,
        bulge_amplitude=args.bulge_amplitude,
        arm_count=args.arm_count,
        pitch_angle_rad=math.radians(args.pitch_angle_deg),
        arm_amplitude=args.arm_amplitude,
        calibration_radius_pc=args.calibration_radius_pc,
    )

    print(
        f"Building skeleton: disk_scale_length_pc={shape.disk_scale_length_pc} "
        f"disk_scale_height_pc={shape.disk_scale_height_pc} "
        f"bulge_scale_radius_pc={shape.bulge_scale_radius_pc} "
        f"bulge_amplitude={shape.bulge_amplitude} arm_count={shape.arm_count} "
        f"k_norm={shape.k_norm:.4f} threshold_rho={threshold_rho:.6f} "
        f"workers={args.workers} chunk_size={args.chunk_size}"
    )

    t0 = time.perf_counter()
    all_bands = []  # (shell_index, band_index, slot_index_min, slot_index_max)
    empty_streak = 0
    outer_shell_index = -1
    shell = 0
    shells_scanned = 0

    with multiprocessing.Pool(processes=args.workers) as pool:
        while shell <= args.max_shell:
            chunk = range(shell, min(shell + args.chunk_size, args.max_shell + 1))
            results = pool.starmap(
                _worker_find_bands,
                [(shape, edge_pc, k, threshold_rho) for k in chunk],
            )
            stop = False
            for shell_index, bands in results:
                shells_scanned += 1
                if bands:
                    empty_streak = 0
                    outer_shell_index = shell_index
                    for band_index, band in enumerate(bands):
                        all_bands.append((shell_index, band_index, band.slot_index_min, band.slot_index_max))
                else:
                    empty_streak += 1
                    if empty_streak >= args.empty_streak_to_stop:
                        stop = True
                        break
            if stop:
                break
            shell = chunk.stop

    # True only if the loop actually broke on a confirmed run of
    # --empty-streak-to-stop consecutive empty shells, not because
    # --max-shell was reached first -- see DEFAULT_MAX_SHELL's docstring.
    edge_confirmed = empty_streak >= args.empty_streak_to_stop

    elapsed = time.perf_counter() - t0

    mysql_config = _db.mysql_config_from_args(args)
    _db.replace_galaxy_shell_bands(all_bands, config=mysql_config)
    _db.save_galaxy_shape(
        shape, edge_pc=edge_pc, outer_shell_index=outer_shell_index,
        expected_system_count_at_density_1=e_value, config=mysql_config,
    )

    total_candidate_slots = sum(b[3] - b[2] + 1 for b in all_bands)
    return {
        "shells_scanned": shells_scanned,
        "outer_shell_index": outer_shell_index,
        "total_bands": len(all_bands),
        "total_candidate_slots": total_candidate_slots,
        "elapsed_s": elapsed,
        "edge_confirmed": edge_confirmed,
    }


def main():
    args = process_args()
    summary = build_skeleton(args)
    print(
        f"Skeleton built in {summary['elapsed_s']:.2f}s: scanned {summary['shells_scanned']} shells, "
        f"outer edge = shell {summary['outer_shell_index']}, {summary['total_bands']} band(s) stored, "
        f"~{summary['total_candidate_slots']:,} candidate sector slots."
    )
    if not summary["edge_confirmed"]:
        print(
            f"WARNING: reached --max-shell ({args.max_shell}) without a run of "
            f"{args.empty_streak_to_stop} consecutive empty shells -- the galaxy's true edge was not "
            f"confirmed. Re-run with a larger --max-shell if these shape parameters really do produce "
            f"a galaxy this large."
        )


if __name__ == "__main__":
    main()
