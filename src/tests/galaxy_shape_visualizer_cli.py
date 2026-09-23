#!/usr/bin/env python3
# tests/galaxy_shape_visualizer_cli.py

"""
Diagnostic tool: renders this galaxy's real density model
(`stellarObjects.galaxyDensity.GalaxyShape`) as actual images, so a
person can eyeball whether a set of shape parameters (`generate.py
plan`'s own args) actually looks like a recognizable spiral galaxy --
face-on and edge-on -- rather than only ever judging it through
`relative_density` numbers or a live web page's own scatter of
illustrative points (`stellarObjects.galaxyViewport.density_sample_points`).

Not a pytest test (no `test_` prefix, so pytest never collects it --
matching this project's existing `climate_tuning_cli.py`/
`phenomena_plausibility_cli.py`/`physical_plausibility_cli.py` precedent
of a diagnostic CLI tool that happens to live under `tests/`, run by hand
rather than automatically).

Two panels are always rendered:

- **Face-on** (top-down, the x/y plane at `z=0`): a `relative_density`
  heatmap over a grid -- the same "how does the spiral look from directly
  above" view a real galaxy photograph shows, and the most direct way to
  confirm the spiral-arm math (`docs/design/galaxy-disk-density.md`)
  actually produces a recognizable spiral for a given set of parameters.
- **Edge-on** (the `r_cyl`/`z` plane, azimuthally averaged -- `arm_factor`
  averages to exactly `1` over a full revolution, see
  `galaxySkeleton._bound_raw_density_at_phi`'s own docstring for the same
  fact used there): a heatmap of the disk's pure bulge+disk vertical
  envelope, independent of any arm structure (which only exists in
  azimuth, so a single fixed-theta slice would show one arbitrary
  in-arm-or-not cross-section instead of the shape's real vertical
  profile).

With `--mysql-*` connection args (or `PLANETGEN_MYSQL_*` env vars, this
project's usual convention -- see `stellarObjects._db`), a third overlay
is added to the face-on panel: every real, already-generated sector's own
stored galaxy position (`sectors.center_x/y/z_pc`), so "does this look
like a spiral" and "is my already-generated ('known') space actually
sitting where the model predicts it should be dense" can both be checked
in the same image. Skipped (with a printed note, not an error) if no
connection was given, or the database has no galaxy-placed sectors yet --
the shape-only panels still render regardless.

Usage:
    python src/tests/galaxy_shape_visualizer_cli.py
    python src/tests/galaxy_shape_visualizer_cli.py --arm-count 4 --pitch-angle-deg 25
    python src/tests/galaxy_shape_visualizer_cli.py --mysql-database mygalaxy -o mygalaxy.png
"""

import argparse
import math
import os
import sys

import matplotlib

matplotlib.use("Agg")  # headless -- this tool only ever writes a PNG, never opens a window
import matplotlib.pyplot as plt
import numpy as np

_SRC_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _SRC_DIR)

from stellarObjects import _db  # noqa: E402
from stellarObjects.galaxyDensity import build_galaxy_shape, relative_density  # noqa: E402

DEFAULT_DISK_SCALE_LENGTH_PC = 2800.0
DEFAULT_DISK_SCALE_HEIGHT_PC = 350.0
DEFAULT_BULGE_SCALE_RADIUS_PC = 200.0
DEFAULT_BULGE_AMPLITUDE = 1.0
DEFAULT_ARM_COUNT = 2
DEFAULT_PITCH_ANGLE_DEG = 15.0
DEFAULT_ARM_AMPLITUDE = 0.4
"""These mirror `generate.py plan`'s own `add_plan_arguments` defaults
exactly, so running this tool with no shape flags at all visualizes the
same default galaxy a bare `generate.py plan` would build."""

FACE_ON_EXTENT_SCALE_LENGTHS = 4.0
"""float: Default face-on half-width, in disk scale lengths -- wide
enough to show a couple of full spiral windings plus genuinely empty
outer space around them, without the bulge swallowing the whole frame
the way a tighter crop would."""

EDGE_ON_R_EXTENT_SCALE_LENGTHS = 2.0
EDGE_ON_Z_EXTENT_SCALE_HEIGHTS = 8.0
"""Edge-on panel extents, in scale lengths/heights respectively -- the
vertical axis needs far more scale-height multiples than the radial axis
needs scale-length multiples to show the disk's own thin/thick-disk
falloff clearly, since a real disk's height is a small fraction of its
radius (see `docs/design/galaxy-coordinate-system.md`'s own "disk
scale height is meaningfully smaller than disk scale length" framing)."""


def _add_shape_arguments(parser):
    """Adds every `GalaxyShape` parameter as a CLI flag, mirroring
    `generate.py`'s own `add_plan_arguments` names/help text/defaults
    exactly, so the two tools' output can be compared apples-to-apples."""
    parser.add_argument("--disk-scale-length-pc", type=float, default=DEFAULT_DISK_SCALE_LENGTH_PC,
                        help="Disk radial exponential scale length, parsecs.")
    parser.add_argument("--disk-scale-height-pc", type=float, default=DEFAULT_DISK_SCALE_HEIGHT_PC,
                        help="Disk vertical scale height, parsecs.")
    parser.add_argument("--bulge-scale-radius-pc", type=float, default=DEFAULT_BULGE_SCALE_RADIUS_PC,
                        help="Central bulge exponential scale radius, parsecs.")
    parser.add_argument("--bulge-amplitude", type=float, default=DEFAULT_BULGE_AMPLITUDE,
                        help="Bulge peak density, relative to the disk's own calibration point.")
    parser.add_argument("--arm-count", type=int, default=DEFAULT_ARM_COUNT,
                        help="Number of spiral arms.")
    parser.add_argument("--pitch-angle-deg", type=float, default=DEFAULT_PITCH_ANGLE_DEG,
                        help="Spiral pitch angle, degrees (larger = more tightly wound... "
                             "counterintuitively named the opposite in some conventions -- see "
                             "galaxyDensity.GalaxyShape's own docstring).")
    parser.add_argument("--arm-amplitude", type=float, default=DEFAULT_ARM_AMPLITUDE,
                        help="Spiral arm/inter-arm density contrast, 0-1ish.")
    parser.add_argument("--calibration-radius-pc", type=float, default=None,
                        help="Where relative_density == 1.0 is calibrated -- defaults to "
                             "2.82 * disk-scale-length-pc (build_galaxy_shape's own default).")


def _shape_from_args(args):
    return build_galaxy_shape(
        disk_scale_length_pc=args.disk_scale_length_pc,
        disk_scale_height_pc=args.disk_scale_height_pc,
        bulge_scale_radius_pc=args.bulge_scale_radius_pc,
        bulge_amplitude=args.bulge_amplitude,
        arm_count=args.arm_count,
        pitch_angle_rad=math.radians(args.pitch_angle_deg),
        arm_amplitude=args.arm_amplitude,
        calibration_radius_pc=args.calibration_radius_pc,
    )


def _face_on_grid(shape, extent_pc, resolution):
    """`relative_density` over an `resolution` x `resolution` grid in the
    `z=0` plane, `[-extent_pc, extent_pc]` on both axes -- the face-on
    view."""
    xs = np.linspace(-extent_pc, extent_pc, resolution)
    ys = np.linspace(-extent_pc, extent_pc, resolution)
    grid = np.empty((resolution, resolution))
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            grid[j, i] = relative_density((float(x), float(y), 0.0), shape)
    return xs, ys, grid


def _edge_on_grid(shape, r_extent_pc, z_extent_pc, resolution, n_theta=24):
    """Azimuthally-averaged `relative_density` over an `resolution` x
    `resolution` grid in the `(r_cyl, z)` half-plane, `r_cyl` in
    `[0, r_extent_pc]` and `z` in `[-z_extent_pc, z_extent_pc]` -- the
    edge-on view. Averaging over `n_theta` evenly-spaced angles cancels
    the spiral-arm term's own azimuthal dependence (its mean over a full
    revolution is exactly `1`, the same fact
    `galaxySkeleton._bound_raw_density_at_phi` relies on), leaving the
    pure bulge+disk vertical envelope."""
    r_vals = np.linspace(0.0, r_extent_pc, resolution)
    z_vals = np.linspace(-z_extent_pc, z_extent_pc, resolution)
    grid = np.empty((resolution, resolution))
    thetas = [2.0 * math.pi * k / n_theta for k in range(n_theta)]
    for j, z in enumerate(z_vals):
        for i, r in enumerate(r_vals):
            total = 0.0
            for theta in thetas:
                total += relative_density((float(r) * math.cos(theta), float(r) * math.sin(theta), float(z)), shape)
            grid[j, i] = total / n_theta
    return r_vals, z_vals, grid


def _real_sector_positions(mysql_config):
    """Every galaxy-placed sector's own stored `(x, y, z)`, parsecs --
    `[]` if the connection fails, the schema has no `sectors` rows with a
    galaxy position yet, or `mysql_config` is `None` (no `--mysql-*`
    given). Never raises -- this overlay is optional/best-effort, not a
    reason to fail the whole visualization."""
    if mysql_config is None:
        return []
    try:
        conn = _db.get_connection(mysql_config)
    except Exception as exc:  # noqa: BLE001 -- this overlay is best-effort, never fatal
        print(f"Note: could not connect to a database ({exc}) -- rendering the shape-only panels "
              "without a 'known space' overlay.")
        return []
    try:
        rows = conn.execute(
            "SELECT center_x_pc, center_y_pc, center_z_pc FROM sectors WHERE center_x_pc IS NOT NULL"
        ).fetchall()
    finally:
        conn.close()
    return [(row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]) for row in rows]


def render(args):
    shape = _shape_from_args(args)

    face_extent_pc = args.extent_pc or (FACE_ON_EXTENT_SCALE_LENGTHS * args.disk_scale_length_pc)
    edge_r_extent_pc = EDGE_ON_R_EXTENT_SCALE_LENGTHS * args.disk_scale_length_pc
    edge_z_extent_pc = EDGE_ON_Z_EXTENT_SCALE_HEIGHTS * args.disk_scale_height_pc

    print(f"Rendering face-on panel ({args.resolution}x{args.resolution}, "
          f"+/-{face_extent_pc:,.0f} pc)...")
    face_xs, face_ys, face_grid = _face_on_grid(shape, face_extent_pc, args.resolution)

    print(f"Rendering edge-on panel ({args.resolution}x{args.resolution}, "
          f"r in [0, {edge_r_extent_pc:,.0f}] pc, z in +/-{edge_z_extent_pc:,.0f} pc)...")
    edge_rs, edge_zs, edge_grid = _edge_on_grid(shape, edge_r_extent_pc, edge_z_extent_pc, args.resolution)

    # Always attempted, not gated on --mysql-database being explicitly
    # given -- mysql_config_from_args already falls through to
    # $PLANETGEN_MYSQL_* / MySQLConfig's own defaults, so a deployment
    # that configures its database via env vars gets the overlay for
    # free; _real_sector_positions itself never raises if no server is
    # actually reachable there either.
    mysql_config = _db.mysql_config_from_args(args)
    sector_positions = _real_sector_positions(mysql_config)
    if sector_positions:
        print(f"Overlaying {len(sector_positions)} real, already-generated sector position(s).")
    else:
        print("No galaxy-placed sectors found (or no database reachable) -- face-on panel shows the "
              "shape only.")

    fig, (ax_face, ax_edge) = plt.subplots(1, 2, figsize=(14, 7), facecolor="#0b0e14")
    for ax in (ax_face, ax_edge):
        ax.set_facecolor("#0b0e14")
        ax.tick_params(colors="#c8ccd8")
        for spine in ax.spines.values():
            spine.set_color("#3a3f55")

    # log1p compresses the exponential falloff (bulge/disk both fall off
    # by orders of magnitude across the frame) into a range human vision
    # can actually distinguish -- the same "display-only contrast
    # stretch" reasoning html/lib/galaxymap.py's own (now-removed)
    # _DENSITY_RADIAL_GAMMA used for the flat map's density shading.
    face_display = np.log1p(face_grid)
    ax_face.pcolormesh(face_xs, face_ys, face_display, shading="auto", cmap="inferno")
    ax_face.set_aspect("equal")
    ax_face.set_title("Face-on (top-down, z=0)", color="#e8eaf0")
    ax_face.set_xlabel("x (pc)")
    ax_face.set_ylabel("y (pc)")

    if sector_positions:
        sector_x = [p[0] for p in sector_positions]
        sector_y = [p[1] for p in sector_positions]
        ax_face.scatter(sector_x, sector_y, s=14, c="#5ce1ff", marker="*",
                        edgecolors="#ffffff", linewidths=0.3, label="Known space (generated sectors)", zorder=5)
        ax_face.legend(loc="upper right", facecolor="#0b0e14", edgecolor="#3a3f55", labelcolor="#e8eaf0")

    edge_display = np.log1p(edge_grid)
    ax_edge.pcolormesh(edge_rs, edge_zs, edge_display, shading="auto", cmap="inferno")
    ax_edge.set_title("Edge-on (azimuthally averaged, r_cyl vs z)", color="#e8eaf0")
    ax_edge.set_xlabel("r_cyl (pc)")
    ax_edge.set_ylabel("z (pc)")

    shape_summary = (
        f"disk_scale_length={args.disk_scale_length_pc:,.0f}pc  disk_scale_height={args.disk_scale_height_pc:,.0f}pc  "
        f"bulge_scale_radius={args.bulge_scale_radius_pc:,.0f}pc  bulge_amplitude={args.bulge_amplitude:g}  "
        f"arms={args.arm_count}  pitch={args.pitch_angle_deg:g}deg  arm_amplitude={args.arm_amplitude:g}"
    )
    fig.suptitle(shape_summary, color="#c8ccd8", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(args.output, dpi=150, facecolor=fig.get_facecolor())
    print(f"Saved to {args.output}")


def main():
    parser = argparse.ArgumentParser(
        description="Renders this galaxy's real density model as face-on/edge-on images -- see this "
                     "file's own module docstring for what each panel shows.",
    )
    _add_shape_arguments(parser)
    parser.add_argument("--extent-pc", type=float, default=None,
                        help="Face-on panel half-width, parsecs -- defaults to "
                             f"{FACE_ON_EXTENT_SCALE_LENGTHS:g} * --disk-scale-length-pc.")
    parser.add_argument("--resolution", type=int, default=320,
                        help="Grid resolution per panel (resolution x resolution samples). Higher is "
                             "sharper but slower -- relative_density is cheap, so even 600+ finishes in "
                             "a few seconds.")
    parser.add_argument("-o", "--output", type=str, default="galaxy_shape.png",
                        help="Output PNG path.")
    _db.add_mysql_connection_args(parser)
    args = parser.parse_args()

    if args.resolution < 10:
        parser.error("--resolution must be at least 10.")

    render(args)


if __name__ == "__main__":
    main()
