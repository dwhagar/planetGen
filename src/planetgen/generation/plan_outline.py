# planetgen/generation/plan_outline.py

"""
The galaxy a `plan` would make, worked out without storing it (GEN.617): the shape drawn from the dials and the seed,
its components and its outline, and the line that says how dense the galaxy is and how much of it the outline holds.

`plan-skeleton` builds its outline here and stores it; the Generate page asks `preview` for the same numbers while the
admin sets the dials, so the line under the dials and the line in the plan's output cannot differ.
"""

import argparse
import dataclasses
import math
import time

from planetgen import tuning as program_constants
from planetgen.galaxy import component_field as componentField
from planetgen.galaxy import dials
from planetgen.galaxy import shape_draw as shapeDraw
from planetgen.galaxy.density import model_terms
from planetgen.galaxy.skeleton import (
    DEFAULT_REACH_SHELL, build_shell_extents, expected_system_count_at_density_1, shell_candidate_sector_count,
)
from planetgen.generation import clusters, satellites
from planetgen.generation import open_clusters as openClusters

PREVIEW_SEED = bytes(range(16))
"""bytes: The seed a preview draws the Auto levers from when the galaxy has none yet (the real run draws its own)."""


@dataclasses.dataclass
class Outline:
    """A planned galaxy: `shape` (with its components), `drawn` (`shape_draw.GalaxyDraw`), the components found, the
    outline's `extents`, `outer_shell_index` and `edge_confirmed` (`skeleton.build_shell_extents`), and the edge,
    one-star density and `shape_args` it was built for."""
    shape: object
    drawn: object
    clusters: list
    young: list
    satellites: list
    extents: list
    outer_shell_index: int
    edge_confirmed: bool
    edge_pc: float
    e_value: float
    shape_args: object
    galaxy_seed: bytes
    elapsed_s: float


def build(shape_args, galaxy_seed, edge_pc, reach_shell_index=DEFAULT_REACH_SHELL):
    """
    The galaxy `shape_args` (the lever table's parsed values, `None` for Auto) and `galaxy_seed` plan: the shape, its
    components and its outline. Nothing is stored.

    Raises:
        ValueError: A lever is outside its range, or the Edge dial is not a length the galaxy can use.
    """
    e_value = expected_system_count_at_density_1(program_constants.DEFAULT_SECTOR_EDGE_LY)
    shape, drawn = shapeDraw.plan_shape(shape_args, galaxy_seed)
    scale = shape.disk_scale_length_pc / program_constants.HALO_SCALE_LENGTH_PC
    epoch = program_constants.SCATTER_PHASE_ORIGIN_UNIX
    found = clusters.generate_clusters(galaxy_seed, epoch_unix=epoch, scale=scale)
    dwarfs = satellites.place_satellites(
        galaxy_seed, program_constants.REFERENCE_RADIUS_SCALE_LENGTHS * shape.disk_scale_length_pc,
        model_terms(shape)["reference_angle_rad"], epoch_unix=epoch, scale=scale)
    young = openClusters.generate_open_clusters(galaxy_seed, shape, edge_pc, e_value, epoch_unix=epoch)
    shape = shape._replace(components=componentField.ComponentField(found + young + dwarfs, edge_pc, e_value))
    edge_limit = shapeDraw.edge_limit_pc(dials.resolved(shape_args, "edge"), shape.disk_scale_length_pc)
    max_shell = reach_shell_index if edge_limit is None else min(math.ceil(edge_limit / edge_pc), reach_shell_index)
    started = time.perf_counter()
    extents, outer_shell_index, edge_confirmed = build_shell_extents(
        shape, edge_pc, 1.0 / e_value, max_shell=max_shell)
    elapsed = time.perf_counter() - started
    return Outline(shape, drawn, found, young, dwarfs, extents, outer_shell_index, edge_confirmed, edge_pc, e_value,
                   shape_args, galaxy_seed, elapsed)


EMPTY_OUTLINE_REASON = (
    "no sector of this galaxy reaches the one-star-per-sector density, so there is nothing to plan: raise the "
    "stellar mass or the central surface density, or lower the scale length or the scale height")


@dataclasses.dataclass
class Report:
    """What an outline says about its galaxy: `crest` (the thin disc's central density over the one-star density, the
    density-limits report's D), `share` (the share of the stars in planned sectors, approximate), `shells`, `sectors`
    (candidate sectors) and `radius_pc` of the outline, `empty` and the `warnings` (sentences)."""
    crest: float
    share: float
    shells: int
    sectors: int
    radius_pc: float
    empty: bool
    warnings: list

    def line(self):
        """The one line the plan prints and the Generate page shows under the dials."""
        if self.empty:
            return f"Outline: empty (the disc's crest is {self.crest:.2g} times the one-star density)"
        return (f"Outline: the disc's crest is {self.crest:.3g} times the one-star density, about "
                f"{100 * self.share:.0f} percent of the stars sit in planned sectors, {self.shells:,} shells out to "
                f"{self.radius_pc / 1000:,.1f} kpc, {self.sectors:,} candidate sectors")


def report(outline):
    """The `Report` of a built `outline`."""
    drawn = outline.drawn
    sigma0 = shapeDraw.drawn_surface_density(drawn)
    crest = dials.disc_crest(sigma0, drawn.disk_scale_height_pc)
    empty = not outline.extents
    radius = 0.0 if empty else (outline.outer_shell_index + 1) * outline.edge_pc
    return Report(
        crest=crest, share=0.0 if empty else dials.planned_share(crest), shells=len(outline.extents),
        sectors=shell_candidate_sector_count(outline.extents), radius_pc=radius, empty=empty,
        warnings=dials.shape_warnings(sigma0, drawn.disk_scale_height_pc, drawn.bulge_scale_radius_pc))


def _parser():
    parser = argparse.ArgumentParser(add_help=False, exit_on_error=False)
    dials.add_arguments(parser)
    parser.add_argument("--reach-shells", type=int, default=DEFAULT_REACH_SHELL)
    return parser


def preview(argv, galaxy_seed=None, edge_pc=program_constants.DEFAULT_SECTOR_EDGE_PC):
    """
    The report the Generate page shows for the dials in `argv` (the flags of the `plan` command, as the page builds
    them).

    Args:
        argv (list): The lever flags and values.
        galaxy_seed (bytes or None): The galaxy's seed; none draws the Auto levers from `PREVIEW_SEED`.
        edge_pc (float): The sector edge.

    Returns:
        tuple: `(report, refusal)`: the `Report` and `None`, or the report of an empty outline and why the plan refuses.

    Raises:
        ValueError: The flags are not valid.
    """
    try:
        args = _parser().parse_args(argv)
    except (argparse.ArgumentError, SystemExit) as exc:
        raise ValueError(str(exc)) from exc
    built = build(args, galaxy_seed or PREVIEW_SEED, edge_pc, args.reach_shells)
    made = report(built)
    return made, EMPTY_OUTLINE_REASON if made.empty else None
