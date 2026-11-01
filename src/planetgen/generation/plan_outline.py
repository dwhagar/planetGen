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
from planetgen.galaxy import plan_figures as planFigures
from planetgen.galaxy import shape_draw as shapeDraw
from planetgen.galaxy.skeleton import (
    DEFAULT_REACH_SHELL, build_shell_extents, expected_system_count_at_density_1, shell_candidate_sector_count,
)
from planetgen.generation import clusters, satellites
from planetgen.generation import open_clusters as openClusters
from planetgen.util import step_timer as stepTimer

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
    timer: object


def build(shape_args, galaxy_seed, edge_pc, reach_shell_index=DEFAULT_REACH_SHELL):
    """
    The galaxy `shape_args` (the lever table's parsed values, `None` for Auto) and `galaxy_seed` plan: the shape, its
    components and its outline. Nothing is stored.

    Raises:
        ValueError: A lever is outside its range, or the Edge dial is not a length the galaxy can use.
    """
    timer = stepTimer.StepTimer()
    e_value = expected_system_count_at_density_1(program_constants.DEFAULT_SECTOR_EDGE_LY)
    with timer.step("the galaxy's shape (the levers, the halo, the disc's history)"):
        shape, drawn = shapeDraw.plan_shape(shape_args, galaxy_seed)
    scale = shape.disk_scale_length_pc / program_constants.HALO_SCALE_LENGTH_PC
    epoch = program_constants.SCATTER_PHASE_ORIGIN_UNIX
    with timer.step("the globular clusters"):
        count = clusters.cluster_count(galaxy_seed, drawn.chain.halo_mass_solar)
        found = clusters.generate_clusters(galaxy_seed, count, epoch_unix=epoch, scale=scale, shape=shape)
    with timer.step("the satellite galaxies"):
        dwarfs = satellites.draw_satellites(galaxy_seed, drawn.chain, epoch_unix=epoch)
    with timer.step("the young open clusters"):
        young = openClusters.generate_open_clusters(galaxy_seed, shape, edge_pc, e_value, epoch_unix=epoch)
    with timer.step("the component field (each cluster's and satellite's systems)"):
        shape = shape._replace(components=componentField.ComponentField(found + young + dwarfs, edge_pc, e_value))
    edge_limit = shapeDraw.edge_limit_pc(dials.resolved(shape_args, "edge"), shape.disk_scale_length_pc)
    max_shell = reach_shell_index if edge_limit is None else min(math.ceil(edge_limit / edge_pc), reach_shell_index)
    with timer.step("the outline (the skeleton's shells)"):
        started = time.perf_counter()
        extents, outer_shell_index, edge_confirmed = build_shell_extents(
            shape, edge_pc, 1.0 / e_value, max_shell=max_shell)
        elapsed = time.perf_counter() - started
    return Outline(shape, drawn, found, young, dwarfs, extents, outer_shell_index, edge_confirmed, edge_pc, e_value,
                   shape_args, galaxy_seed, elapsed, timer)


EMPTY_OUTLINE_REASON = (
    "no sector of this galaxy reaches the one-star-per-sector density, so there is nothing to plan: raise the "
    "stellar mass or the central surface density, or lower the scale length or the scale height")


@dataclasses.dataclass
class Report:
    """What an outline says about its galaxy: `crest` (the thin disc's expected systems per sector on an arm crest at
    its centre, which is its central density over the one-star density, the density-limits report's D), `share` (the
    share of the galaxy's systems in sectors that expect one, integrated from its density; approximate), `shells`,
    `sectors` (candidate sectors, the outline's upper bound) and `reaching` (the sectors that expect at least one system,
    integrated from the density; approximate, GEN.627) and `radius_pc` of the outline, `empty` and the `warnings`
    (sentences)."""
    crest: float
    share: float
    shells: int
    sectors: int
    reaching: int
    radius_pc: float
    empty: bool
    warnings: list

    def line(self):
        """The one line the plan prints and the Generate page shows under the dials."""
        if self.empty:
            return f"Outline: empty (the disc's crest is {self.crest:.2g} times the one-star density)"
        return (f"Outline: the disc's crest is {self.crest:.3g} times the one-star density, about "
                f"{100 * self.share:.0f} percent of the stars sit in planned sectors, {self.shells:,} shells out to "
                f"{self.radius_pc / 1000:,.1f} kpc, {self.sectors:,} candidate sectors, about {self.reaching:,} of them expect "
                f"at least one system")


def report(outline):
    """The `Report` of a built `outline`."""
    drawn = outline.drawn
    crest = planFigures.crest(outline.shape, outline.e_value)
    empty = not outline.extents
    radius = 0.0 if empty else (outline.outer_shell_index + 1) * outline.edge_pc
    figures = None if empty else planFigures.planned(outline.shape, outline.e_value, radius)
    share = 0.0 if empty else figures.share
    reaching = 0 if empty else round(planFigures.reaching_sectors(
        outline.shape, outline.e_value, radius, outline.edge_pc, figures=figures))
    return Report(
        crest=crest, share=share, shells=len(outline.extents),
        sectors=shell_candidate_sector_count(outline.extents), reaching=reaching, radius_pc=radius, empty=empty,
        warnings=dials.shape_warnings(crest, share, drawn.disk_scale_height_pc, drawn.bulge_scale_radius_pc))


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
