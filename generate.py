#!/usr/bin/env python3
# generate.py

"""
Unified Generation CLI
=========================

The single entry point for every generator in this project -- star
systems, sectors, whole galaxies, the galaxy's density skeleton, and
exotic stellar phenomena all live in this one file now, dispatched by
subcommand:

    generate.py system [options]      -- one star system
    generate.py sector [options]      -- one or more independent sectors
    generate.py galaxy [options]      -- many sectors placed as one galaxy
    generate.py plan [options]        -- the galaxy's density skeleton
    generate.py phenomenon [options]  -- one exotic stellar phenomenon
    generate.py population [options]  -- species, civilizations and territories

Run `generate.py <command> --help` for that command's own full option
list. This replaces the five separate scripts this project used to ship
(`systemGen.py`, `sectorGen.py`, `galaxyGen.py`, `galaxyPlan.py`,
`phenomenonGen.py`) -- their logic now lives here directly rather than
being spread across five files that imported each other (`sectorGen.py`
called into `systemGen.py`, `galaxyGen.py` called into `sectorGen.py`,
and so on); every option and behavior those scripts had is preserved
exactly, just reachable through one program and one subcommand instead
of five separate ones.

Sections below, in the same layering the old scripts had (later
sections build on earlier ones):

1. System generation (`system`) -- `SystemConfig`/`--system-file`
   handling, tri-state `+name`/`-name` options, `build_system_config`.
2. Sector generation (`sector`) -- builds one `SystemConfig` per system
   via section 1, adds a sparse population of exotic phenomena, and
   places every system in a shared `SpaceSector`.
3. Galaxy generation (`galaxy`) -- calls section 2's `generate_sector`
   per sector, additionally supplying a real galaxy-frame position
   (ring/layer/slot address in the cylindrical grid,
   `galactic_center_dist_ly`) and persisting that position, with systems
   placed inside the sector's own cell. Four modes: `--ring` (batch),
   `--ring --slot` (one address), `--center-sector`/`--radius-pc` (local
   neighborhood), or neither (random start). `ensure_sector_generated`
   is this same per-sector logic exposed as a non-CLI, visit-triggered
   entry point against the galaxy skeleton section 4 builds.
4. Galaxy density skeleton (`plan`) -- the compact `galaxy_shape`/
   `galaxy_layer` summary (the galaxy's outline, one row per layer) `ensure_sector_generated` consults to
   decide, cheaply and exactly, whether a given address is worth
   generating at all, without ever enumerating the galaxy's ~10 billion
   candidate sector slots.
5. Exotic phenomena (`phenomenon`) -- black holes, neutron stars,
   nebulae, supernova remnants, rogue planets, interstellar comets, and
   standalone asteroid fields, generated on demand; section 2 also
   reuses this section's `generate_phenomenon` to seed every sector with
   its own sparse, science-based population of the same seven types.
6. Population and politics (`population`) -- species, civilizations
   and territories from what is already stored
   (`stellarObjects/population.py`); also run after a `sector` or
   `galaxy` run given `--population`.
7. The unified CLI itself (argument parsing/validation, dispatch,
   `main`).
"""

import argparse
import copy
import getpass
import logging
import math
import os
import random
import re
import secrets
import sys
import time
from collections import Counter

import pymysql
from rich.progress import (
    BarColumn, MofNCompleteColumn, Progress, TextColumn, TimeElapsedColumn, TimeRemainingColumn,
)

# stellarObjects lives at src/stellarObjects (src layout) -- add src/ to the
# import path so this keeps working without requiring `pip install .` first.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))

from stellarObjects import (
    _db, activitylog, brightStars, generationLimits, log, physical_constants, population, program_constants, progressFile,
    workQueue,
)
from stellarObjects._version import VersionAction, version_banner
from stellarObjects.asteroidFieldData import AsteroidField
from stellarObjects.compactRemnant import BlackHole, NeutronStar
from stellarObjects.config import SystemConfig
from stellarObjects.galaxyDensity import build_galaxy_shape, predicted_star_count, relative_density
from stellarObjects.galaxyDrill import (
    DRILL_LEVELS, DrillBlock, drill_block_sectors, drill_children, drill_parent, drill_slabs, format_drill_key,
    parse_drill_key,
)
from stellarObjects.galaxyGeometry import (
    SectorCell, enumerate_sectors_within_radius, galactic_radius_pc,
    provisional_sector_designation, ring_bounds_pc, ring_sector_count, sector_address_at, sector_position_pc,
)
from stellarObjects.galaxySkeleton import (
    DEFAULT_MAX_RING, build_layer_extents, candidate_sector_count, expected_system_count_at_density_1,
)
from stellarObjects.nebulaData import Nebula, choose_weighted_class
from stellarObjects.quasarData import Quasar
from stellarObjects.roguePlanetData import InterstellarComet, RoguePlanet
from stellarObjects.spaceSector import SpaceSector, _sample_poisson_count
from stellarObjects.starData import STAR_TYPE_PATTERN
from stellarObjects.stellarPopulation import bright_star_fraction
from stellarObjects.supernovaRemnantData import SupernovaRemnant
from stellarObjects.systemData import StarSystem
from stellarObjects.systemRender import render_star_system
from stellarObjects.utils import generate_sector_name, ly_to_pc, pc_to_ly

# Suppress transformers warnings
logging.getLogger("transformers").setLevel(logging.ERROR)


class _ReportingProgress(Progress):
    """`rich.progress.Progress` that also mirrors its most recently
    changed task to `stellarObjects.progressFile` (a no-op unless the web
    interface started this run)."""

    def _report(self, task_id, force=False):
        task = self._tasks.get(task_id)
        if task is not None:
            progressFile.report(task.completed, task.total, task.description, force=force)

    def add_task(self, description, *args, **kwargs):
        task_id = super().add_task(description, *args, **kwargs)
        self._report(task_id, force=True)
        return task_id

    def update(self, task_id, **kwargs):
        super().update(task_id, **kwargs)
        self._report(task_id, force=kwargs.get("total") is not None)

    def advance(self, task_id, advance=1):
        super().advance(task_id, advance)
        self._report(task_id)

    def stop(self):
        for task_id in list(self._tasks):
            self._report(task_id, force=True)
        super().stop()


def _generation_progress():
    """
    Builds the shared `rich.progress.Progress` used by `run_galaxy`'s three
    modes -- one "Sectors" task per run tracking how many sectors have been
    generated so far. Deliberately not used by `run_sector`'s own
    `--num-sectors` loop: that command's own run is normally short enough
    (and its sector count small enough) that a bar added more visual noise
    than it was worth, whereas a `galaxy` run (a whole ring, a
    neighborhood, or a random start's default 12 pc one) can mean
    thousands of sectors and legitimately benefit from a progress display.
    There also used to be a second, nested bar for "systems in the current
    sector," added/removed once per sector; that one is gone for good --
    most sectors hold anywhere from zero to a handful of systems, and
    system generation itself is fast, so a bar that flashed on and off
    again within a single frame for nearly every sector was pure noise,
    not something worth reviving alongside this one.

    Every task shows both elapsed time and an estimated time remaining
    (`TimeElapsedColumn`/`TimeRemainingColumn`) -- the remaining estimate
    only becomes accurate once a task has advanced enough for rich's own
    rate estimate to settle, same as any ETA.

    Callers use this as a context manager (`with _generation_progress() as
    progress:`); `rich.progress.Progress` is a `Live` display under the
    hood, so **every status line logged while it's active must go through
    `progress.console`, never a raw stdout write** -- printing directly to
    stdout fights with the `Live` region's own redraws (each plain `print`
    call forces the bar to erase itself, scroll up with the new text, and
    get redrawn at the bottom again), which is exactly what caused the
    flicker/scrolling a real terminal used to show once before. `run_galaxy`
    calls `log.set_console(progress.console)` right after opening this
    context manager (and `log.reset_console()` once it's closed), so every
    `log.normal(...)`/`log.debug(...)` call made anywhere during a galaxy
    run -- including deep inside `stellarObjects` modules -- is
    automatically routed through `progress.console.print(...)` instead of a
    raw stdout write for as long as the bar is live. Routed that way, rich
    prints each line safely *above* the live region and leaves the bar
    itself pinned at the bottom, redrawn in place with no flicker. This only
    matters when stdout is a real interactive terminal in the first place --
    `Console` auto-detects that (`Console.is_terminal`) and falls back to
    plain, periodic line-by-line bar output otherwise (piped to a file, a CI
    log, etc.), so no separate handling is needed for that case.

    The same counts also go to `$PLANETGEN_PROGRESS_FILE` when it is set
    (`stellarObjects.progressFile`), so the web interface's Generate page
    can show a run it started in the background.

    Returns:
        Progress: Not yet started.
    """
    return _ReportingProgress(
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TextColumn("[dim]elapsed"),
        TimeElapsedColumn(),
        TextColumn("[dim]remaining"),
        TimeRemainingColumn(),
    )


# ===========================================================================
# 1. System generation
# ===========================================================================

# Tri-state options settable with `+name` (force True) / `-name` (force False).
# Each maps directly onto a `SystemConfig` attribute of the same name.
TRISTATE_OPTIONS = [
    ("habitable_world", "HABITABLE_WORLD", "a habitable world"),
    ("asteroid_belt", "ASTEROID_BELT", "an asteroid belt"),
    ("comets", "COMETS", "a star-bound comet"),
    ("large_star", "LARGE_STAR", "a large, massive star"),
    ("moons", "MOONS", "moons on every planet"),
    ("max_planets", "MAX_PLANETS", "the maximum number of orbital objects"),
    ("intelligent_life", "INTELLIGENT_LIFE", "intelligent life on a planet"),
    ("binary_system", "BINARY_SYSTEM", "a binary star system"),
    ("wide_binary", "WIDE_BINARY", "an S-type (wide) binary instead of a P-type (close) one"),
    ("planets", "PLANETS", "at least one planet or asteroid belt"),
]


class TristateAction(argparse.Action):
    """
    Sets `namespace.dest` to True when invoked as `+name`, or False when
    invoked as `-name`. Leaving the option off the command line leaves the
    `default` (None) in place, meaning "let the generator decide".
    """

    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, option_string.startswith('+'))


def finite_float(text):
    """
    `argparse` type for every float option: a plain `float`, but NaN and
    +/-inf (`nan`, `inf`, `1e309`) are rejected at parse time. Every
    validator below uses `x <= 0`-style checks, which NaN (every
    comparison False) slips straight through -- `--density nan` used to
    hang forever in `_sample_poisson_count`, `--radius-pc inf` crashed in
    `math.ceil`, and NaN/inf shape parameters reached MySQL.
    """
    try:
        value = float(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"invalid float value: {text!r}") from None
    if not math.isfinite(value):
        raise argparse.ArgumentTypeError(f"must be a finite number, got {text!r}")
    return value


def _validate_star_type(args, parser):
    """--star-type must be a whole spectral type (e.g. G2V), checked here
    rather than surfacing as a traceback from `Star.generate_star`."""
    if args.star_type and not STAR_TYPE_PATTERN.fullmatch(args.star_type.upper()):
        parser.error(f"--star-type {args.star_type!r} is not a spectral type: expected a class letter "
                     f"(OBAFGKM), a subclass digit (0-9) and a Yerkes class (0, IA+, IA, IAB, IB, II, "
                     f"III, IV, V, VI, VII or D), e.g. G2V.")


def add_logging_arguments(parser):
    """
    Adds the logging options every subcommand shares to `parser`:
    `--debug` and `--quiet`/`--silent`.

    `--debug` alone enables debug-severity output to the console; giving it
    a filename (`--debug FILE`) additionally mirrors that same output to
    `FILE`. `--quiet`/`--silent` (two spellings of the same flag) restrict
    output to errors only. `validate_logging_args` enforces the one
    combination that doesn't make sense: `--quiet`/`--silent` together with
    a filename-less `--debug`, which would have nowhere to send debug
    output once the console is silenced.

    Args:
        parser (argparse.ArgumentParser): The parser to add options to.
    """
    parser.add_argument('--debug', nargs='?', const='', default=None, metavar='FILE',
                        help="Log every choice the generator makes, and why, with timestamps, to the "
                             "console. If FILE is given, also mirror that output to FILE.")
    parser.add_argument('--quiet', '--silent', dest='quiet', action='store_true',
                        help="Suppress all output except errors.")


def validate_logging_args(args, parser):
    """
    Validates the logging options `add_logging_arguments` added, calling
    `parser.error` (which exits) if `--quiet`/`--silent` is combined with a
    filename-less `--debug`.

    Args:
        args (argparse.Namespace): Parsed arguments.
        parser (argparse.ArgumentParser): The parser to raise errors
                                          through (so the caller's own
                                          `--help`/usage text is shown).
    """
    if args.quiet and args.debug is not None and not args.debug:
        parser.error("--debug requires a FILE argument when combined with --quiet/--silent "
                     "(there's nowhere else for debug output to go).")


def add_system_arguments(parser):
    """
    Adds every option the `system` subcommand accepts (besides
    `--version`, which only the top-level parser offers) to `parser` --
    the tri-state `+name`/`-name` flags (see `TRISTATE_OPTIONS`),
    `--system-file`, the MySQL connection args, `--markdown`,
    `--star-type`, `--num-orbits`, `--name`, `--age`, the flavor-text
    overrides, and the logging options (see `add_logging_arguments`).

    Args:
        parser (argparse.ArgumentParser): The parser to add options to.
                                          Must accept `+`/`-` prefix chars
                                          (`prefix_chars='-+'`) for the
                                          tri-state options.
    """
    for name, _attr, description in TRISTATE_OPTIONS:
        parser.add_argument(f'-{name}', f'+{name}', dest=name, action=TristateAction,
                            nargs=0, default=None,
                            help=f"+{name} forces the system to have {description}; "
                                 f"-{name} forces the system to not have {description}.")

    # Load system options from a JSON file
    parser.add_argument('--system-file', '-f', type=str,
                        help="Load system generation options from a JSON file. Command-line options "
                             "override the values it sets.")

    # Database persistence
    _db.add_mysql_connection_args(parser)

    # Output in Markdown format
    parser.add_argument('--markdown', '-m', action='store_true', help="Output in Markdown format.")

    # Write the page instead of saving (a one-off system)
    parser.add_argument('--output', '-o', type=str, metavar='FILE',
                        help="Write the system's page (wikitext, or Markdown with --markdown) to FILE "
                             "instead of saving the system to the database; '-' writes it to stdout.")

    # Logging (--debug, --quiet/--silent)
    add_logging_arguments(parser)

    # Star Type
    parser.add_argument('--star-type', type=str,
                        help="Force the generation of a specific star type (e.g., G2V).")

    # Number of Orbits
    parser.add_argument('--num-orbits', type=int,
                        help="Force an exact number of orbital slots (planets and asteroid belts "
                             "combined) to be generated.")

    # System Name
    parser.add_argument('--name', type=str, help="Force the name of the star system.")

    # System Age
    parser.add_argument('--age', type=str, choices=['young', 'old'],
                        help="Specify the age of the star system (young or old).")

    # Override Flavor Chance System
    parser.add_argument('--flavor-chance-system', type=finite_float,
                        help="Override the default FLAVOR_CHANCE_SYSTEM constant.")

    # Override Flavor Chance Planet
    parser.add_argument('--flavor-chance-planet', type=finite_float,
                        help="Override the default FLAVOR_CHANCE_PLANET constant.")

    # Max Planet Flavor
    parser.add_argument('--max-planet-flavor', action='store_true',
                        help="Sets the maximum flavor text total for planets to 99.")


def validate_system_args(args, parser):
    """
    Validates every option `add_system_arguments` added, calling
    `parser.error` (which exits) on the first problem found.

    Args:
        args (argparse.Namespace): Parsed arguments.
        parser (argparse.ArgumentParser): The parser to raise errors
                                          through (so the caller's own
                                          `--help`/usage text is shown).
    """
    if args.planets is False and (args.moons or args.max_planets or args.habitable_world):
        parser.error("-planets cannot be combined with +moons, +max_planets, or +habitable_world.")

    if args.star_type and args.large_star:
        parser.error("--star-type cannot be combined with +large_star.")

    _validate_star_type(args, parser)

    if args.intelligent_life is not None and args.habitable_world is False:
        parser.error("+intelligent_life/-intelligent_life cannot be combined with -habitable_world.")

    if args.system_file:
        # Read once here so a missing/unreadable/non-JSON/non-object file
        # is a usage error, not a traceback from build_system_config.
        try:
            file_data = load_system_file(args.system_file)
        except (OSError, ValueError) as exc:
            parser.error(f"--system-file {args.system_file!r} could not be read as JSON: {exc}")
        if not isinstance(file_data, dict):
            parser.error(f"--system-file {args.system_file!r} must contain a JSON object, "
                         f"not {type(file_data).__name__}.")
        file_star_type = file_data.get("star_type")
        if file_star_type is not None and not (
                isinstance(file_star_type, str) and STAR_TYPE_PATTERN.fullmatch(file_star_type.upper())):
            parser.error(f"--system-file {args.system_file!r}: star_type {file_star_type!r} is not a "
                         f"spectral type (e.g. G2V).")
        file_num_orbits = file_data.get("num_orbits")
        if file_num_orbits is not None and not (
                isinstance(file_num_orbits, int) and not isinstance(file_num_orbits, bool)
                and 0 <= file_num_orbits <= generationLimits.MAX_NUM_ORBITS):
            parser.error(f"--system-file {args.system_file!r}: num_orbits {file_num_orbits!r} must be a "
                         f"whole number from 0 to {generationLimits.MAX_NUM_ORBITS}.")

    if args.num_orbits is not None and args.num_orbits < 0:
        parser.error("--num-orbits must be zero or a positive integer.")
    if args.num_orbits is not None and args.num_orbits > generationLimits.MAX_NUM_ORBITS:
        parser.error(f"--num-orbits must be at most {generationLimits.MAX_NUM_ORBITS}.")

    if args.num_orbits is not None and args.planets is False:
        parser.error("--num-orbits cannot be combined with -planets.")

    if args.flavor_chance_system is not None and not (0.0 <= args.flavor_chance_system <= 1.0):
        parser.error("--flavor-chance-system must be a float between 0.0 and 1.0.")

    if args.flavor_chance_planet is not None and not (0.0 <= args.flavor_chance_planet <= 1.0):
        parser.error("--flavor-chance-planet must be a float between 0.0 and 1.0.")


def load_system_file(path):
    """
    Loads a system generation specification from a JSON file.

    The JSON file may contain any of the following keys, each corresponding
    to a `SystemConfig` attribute of the same name (see that class's
    docstrings for details): `star_type`, `name`, `age`, `num_orbits`,
    `slots`, `habitable_world`, `asteroid_belt`, `comets`, `large_star`, `moons`,
    `max_planets`, `intelligent_life`, `binary_system`, `wide_binary`,
    `planets`, `markdown`, `flavor_chance_system`, `flavor_chance_planet`,
    `max_planet_flavor`, `output`.

    `slots`, if present, is a list whose entries are either `null` or an
    object with a required `"type"` ("planet" or "asteroid_belt") and, for
    planets, optional `"planet_class"` (e.g. "M") and `"moons"` (an exact
    moon count) keys, e.g.:

        {
          "star_type": "G2V",
          "num_orbits": 5,
          "habitable_world": true,
          "slots": [
            {"type": "planet", "planet_class": "M", "moons": 1},
            {"type": "asteroid_belt"},
            null,
            {"type": "planet", "planet_class": "J", "moons": 4},
            null
          ]
        }

    Args:
        path (str): Path to the JSON system specification file.

    Returns:
        dict: The parsed JSON content.
    """
    import json
    with open(path, 'r') as f:
        return json.load(f)


def apply_system_file(system_config, data):
    """
    Applies a system specification (as loaded by `load_system_file`) onto a
    `SystemConfig` instance.

    Args:
        system_config (SystemConfig): The config object to update in place.
        data (dict): The parsed JSON system specification.
    """
    simple_keys = [
        "star_type", "name", "age", "num_orbits", "slots",
        "habitable_world", "asteroid_belt", "comets", "large_star", "moons",
        "max_planets", "intelligent_life", "binary_system", "wide_binary",
        "planets", "markdown",
    ]
    key_to_attr = {key: key.upper() for key in simple_keys}

    for key, attr in key_to_attr.items():
        if key in data:
            setattr(system_config, attr, data[key])

    if "flavor_chance_system" in data and data["flavor_chance_system"] is not None:
        program_constants.FLAVOR_CHANCE_SYSTEM = data["flavor_chance_system"]

    if "flavor_chance_planet" in data and data["flavor_chance_planet"] is not None:
        program_constants.FLAVOR_CHANCE_PLANET = data["flavor_chance_planet"]

    if data.get("max_planet_flavor"):
        program_constants.MAX_FLAVOR_TOTAL = 99
        program_constants.FLAVOR_CHANCE_PLANET = 1


def build_system_config(args):
    """
    Builds a fully-configured `SystemConfig` from parsed command-line
    arguments: a `--system-file` JSON spec (if given) sets the baseline,
    then any explicitly-given tri-state or value option on the command
    line overrides it, then the cross-option normalization rules
    (intelligent life implies a habitable world; a habitable world plus
    an asteroid belt implies a large star) are applied last.

    Shared by the `system` and `sector` subcommands (`build_sector_configs`
    calls this once per system in the sector), so the two can never
    silently drift apart on what a given set of options actually means.

    Args:
        args (argparse.Namespace): Parsed arguments (from the `system`
            subcommand, or a namespace shaped the same way, as
            `build_sector_configs` builds for each system in a sector).

    Returns:
        SystemConfig: The resolved configuration.

    Raises:
        SystemExit: If a habitable world and an asteroid belt are both
                   forced while a large star is explicitly forbidden (see
                   TRISTATE_OPTIONS's `large_star`) -- that combination has
                   no room to be satisfied.
    """
    system_config = SystemConfig()

    if args.system_file:
        file_data = load_system_file(args.system_file)
        apply_system_file(system_config, file_data)

    # Command-line tri-state options override anything set by --system-file.
    for name, attr, _description in TRISTATE_OPTIONS:
        value = getattr(args, name)
        if value is not None:
            setattr(system_config, attr, value)

    # Command-line value options override anything set by --system-file.
    if args.markdown:
        system_config.MARKDOWN = True
    if args.star_type is not None:
        system_config.STAR_TYPE = args.star_type
    if args.num_orbits is not None:
        system_config.NUM_ORBITS = args.num_orbits
    if args.name is not None:
        system_config.NAME = args.name
    if args.age is not None:
        system_config.AGE = args.age

    if args.flavor_chance_system is not None:
        program_constants.FLAVOR_CHANCE_SYSTEM = args.flavor_chance_system

    if args.flavor_chance_planet is not None:
        program_constants.FLAVOR_CHANCE_PLANET = args.flavor_chance_planet

    if args.max_planet_flavor:
        program_constants.MAX_FLAVOR_TOTAL = 99
        program_constants.FLAVOR_CHANCE_PLANET = 1

    # Intelligent life requires (and forbids the absence of) a habitable world.
    if system_config.INTELLIGENT_LIFE is not None:
        system_config.HABITABLE_WORLD = True

    if system_config.HABITABLE_WORLD is True and system_config.ASTEROID_BELT is True:
        if system_config.LARGE_STAR is False:
            log.error("Error: forcing both a habitable world and an asteroid belt requires a large "
                      "star; -large_star cannot be combined with +habitable_world and +asteroid_belt.")
            raise SystemExit(1)
        system_config.LARGE_STAR = True

    return system_config


def run_system(args):
    """
    Generates one star system and saves it to the database, or, with
    `--output`, writes its page to a file (or stdout) and touches no
    database at all -- the admin site's one-off system page runs it that
    way (`src/html/web/system_page.py`).

    Args:
        args (argparse.Namespace): Validated arguments (`command ==
            "system"`).
    """
    system_config = build_system_config(args)
    system = StarSystem(system_config=system_config)

    if args.output:
        text = render_star_system(system, "markdown" if system_config.MARKDOWN else "wikitext")
        if args.output == "-":
            sys.stdout.write(text if text.endswith("\n") else text + "\n")
        else:
            with open(args.output, "w", encoding="utf-8") as f:
                f.write(text)
            log.normal(f"Wrote system '{system.name}' to {args.output} (not saved to the database).")
        return

    mysql_config = _db.mysql_config_from_args(args)
    star_system_id = _db.save_system(system, system_config, config=mysql_config)
    RUN_COUNTS["systems"] += 1
    log.normal(f"Saved system '{system.name}' to the database (star_system_id={star_system_id}, "
               f"{mysql_config.database}@{mysql_config.host}:{mysql_config.port}).")


RUN_COUNTS = Counter()
"""Counter: What this run saved (`sectors`, `systems`, `phenomena`), for
the activity log's `generate.finish` line."""


def _count_sector(sector):
    RUN_COUNTS["sectors"] += 1
    RUN_COUNTS["systems"] += len(sector.entries)
    RUN_COUNTS["phenomena"] += len(sector.phenomena)


# ===========================================================================
# 2. Sector generation
# ===========================================================================

SECTOR_PHENOMENON_KINDS = (
    ("black-hole", "black-hole",
     lambda config, dist_ly: BlackHole(config, galactic_center_dist_ly=dist_ly)),
    ("neutron-star", "neutron-star",
     lambda config, dist_ly: NeutronStar(config, galactic_center_dist_ly=dist_ly)),
    ("planetary-nebula", "nebula", lambda config, dist_ly: Nebula(config, nebula_type="planetary")),
    ("molecular-cloud", "nebula", lambda config, dist_ly: Nebula(config, nebula_type="dark")),
    ("supernova-remnant", "supernova-remnant", lambda config, dist_ly: SupernovaRemnant(config)),
    ("rogue-planet", "rogue-planet", lambda config, dist_ly: RoguePlanet(config)),
    ("brown-dwarf", "rogue-planet", lambda config, dist_ly: RoguePlanet(config, mass_bin="brown-dwarf")),
    ("comet", "comet", lambda config, dist_ly: InterstellarComet(config)),
    ("asteroid-field", "asteroid-field", lambda config, dist_ly: AsteroidField(config)),
)
"""tuple: `(program_constants.PHENOMENON_DENSITY_PC3 key, phenomenon type,
factory)` for each kind `generate_sector_phenomena` rolls. A factory takes
`(config, galactic_center_dist_ly)`; only black holes and neutron stars
use the distance (their Hill sphere and galactic orbit). A planetary
nebula is generated around its own new hot white dwarf system
(`_add_planetary_nebula`). Runaway and hypervelocity stars are flags on
ordinary systems (`flag_fast_stars`), and emission and reflection nebulae
grow around a sector's own hot stars (`add_star_hosted_nebulae`). Diffuse
gas (classes A-B) fills about half the disk, so it is the background, not
generated."""


def add_shared_generation_options(parser):
    """
    Adds every per-system generation-tuning option the `sector` and
    `galaxy` subcommands share -- everything `build_sector_configs()`
    (and, via it, `build_system_config`) needs from parsed args -- to
    `parser`. Factored out so the two subcommands' option surfaces can
    never silently drift apart on what a given flag means.

    Also adds the logging options (see `add_logging_arguments`), shared
    identically by every subcommand.

    Deliberately excludes `--name`/`-n` (sector naming is the `sector`
    subcommand's own per-invocation concept; `galaxy` names each
    generated sector itself, once per grid address).

    Args:
        parser (argparse.ArgumentParser): The parser to add options to.
                                          Must accept `+`/`-` prefix chars
                                          (`prefix_chars='-+'`) for the
                                          tri-state options below.
    """
    for name, _attr, description in TRISTATE_OPTIONS:
        parser.add_argument(f'-{name}', f'+{name}', dest=name, action=TristateAction,
                            nargs=0, default=None,
                            help=f"+{name} forces every system in the sector to have {description}; "
                                 f"-{name} forces every system in the sector to not have {description}.")

    parser.add_argument('--num-systems', type=int, default=None,
                        help="The exact number of star systems to generate in the sector. Cannot be "
                             "combined with --density. Defaults to 10 if neither is given.")
    parser.add_argument('--density', type=finite_float, default=None,
                        help="Scale the number of systems generated per sector by this multiplier on "
                             "real local stellar density (1.0 = a realistic sector this size; 2.0 = "
                             "twice as dense; 0.5 = half). The actual count is randomly sampled per "
                             "sector (Poisson-distributed), so it varies run to run and sector to "
                             "sector even at the same density -- this is what lets one sector be made "
                             "meaningfully denser or sparser than another. Cannot be combined with "
                             "--num-systems.")
    parser.add_argument('--min-habitable', type=int, default=0,
                        help="Guarantee at least this many systems in the sector have a habitable world, "
                             "chosen randomly among them, without requiring every system to have one.")
    parser.add_argument('--markdown', '-m', action='store_true', help="Output in Markdown format.")
    parser.add_argument('--star-type', type=str,
                        help="Force every system's star to a specific type (e.g., G2V).")
    parser.add_argument('--age', type=str, choices=['young', 'old'],
                        help="Specify the age of every star in the sector (young or old).")
    parser.add_argument('--flavor-chance-system', type=finite_float,
                        help="Override the default FLAVOR_CHANCE_SYSTEM constant.")
    parser.add_argument('--flavor-chance-planet', type=finite_float,
                        help="Override the default FLAVOR_CHANCE_PLANET constant.")
    parser.add_argument('--max-planet-flavor', action='store_true',
                        help="Sets the maximum flavor text total for planets to 99.")
    parser.add_argument('--workers', type=int, default=None,
                        help="How many sectors to generate at once, each in its own low-priority worker "
                             "process. Default: 80%% of this machine's cores (one fewer when MySQL runs "
                             "here too), or PLANETGEN_WORKERS; 1 generates one sector at a time in this "
                             "process.")
    parser.add_argument('--population', action='store_true',
                        help="Also run the population pass (species, civilizations, territories) "
                             "after the sectors are saved. Off by default; 'generate.py population' "
                             "runs it any time.")

    add_logging_arguments(parser)


def validate_shared_generation_args(args, parser):
    """
    Validates every option `add_shared_generation_options` added, calling
    `parser.error` (which exits) on the first problem found. Shared by
    the `sector` and `galaxy` subcommands for the same reason
    `add_shared_generation_options` is.

    Args:
        args (argparse.Namespace): Parsed arguments.
        parser (argparse.ArgumentParser): The parser to raise errors
                                          through (so the caller's own
                                          `--help`/usage text is shown).
    """
    if args.density is not None and args.num_systems is not None:
        parser.error("--density cannot be combined with --num-systems.")
    if args.workers is not None and args.workers < 0:
        parser.error("--workers must be 0 (automatic) or more.")

    if args.density is not None and args.density <= 0:
        parser.error("--density must be a positive number.")

    if args.density is None and args.num_systems is None:
        # `galaxy` mode leaves both None rather than defaulting to a flat
        # count here -- run_ring_batch/run_local_neighborhood compute each
        # sector's own --density from the galaxy skeleton's real
        # position-based relative_density instead (see their own
        # docstrings and `_resolve_batch_density`), the same mechanism
        # `ensure_sector_generated` already uses for a single lazily-
        # generated sector, applied uniformly across a whole batch run.
        # `getattr` (not `args.command` directly): `_default_generation_args`
        # calls this same validator against a bare, subcommand-less
        # namespace with no `command` attribute at all -- falling through
        # to the flat default there is harmless, since
        # `ensure_sector_generated` immediately overwrites both fields
        # right after with its own already-computed density anyway.
        if getattr(args, 'command', None) != 'galaxy':
            args.num_systems = 10

    if args.num_systems is not None:
        if args.num_systems < 1:
            parser.error("--num-systems must be a positive integer.")

        if args.min_habitable > args.num_systems:
            parser.error("--min-habitable cannot exceed --num-systems.")

    if args.min_habitable < 0:
        parser.error("--min-habitable cannot be negative.")

    if args.min_habitable > 0 and args.habitable_world is False:
        parser.error("--min-habitable cannot be combined with -habitable_world.")

    if args.planets is False and (args.moons or args.max_planets or args.habitable_world):
        parser.error("-planets cannot be combined with +moons, +max_planets, or +habitable_world.")

    if args.star_type and args.large_star:
        parser.error("--star-type cannot be combined with +large_star.")

    _validate_star_type(args, parser)

    if args.intelligent_life is not None and args.habitable_world is False:
        parser.error("+intelligent_life/-intelligent_life cannot be combined with -habitable_world.")

    if args.flavor_chance_system is not None and not (0.0 <= args.flavor_chance_system <= 1.0):
        parser.error("--flavor-chance-system must be a float between 0.0 and 1.0.")

    if args.flavor_chance_planet is not None and not (0.0 <= args.flavor_chance_planet <= 1.0):
        parser.error("--flavor-chance-planet must be a float between 0.0 and 1.0.")


def add_sector_arguments(parser):
    """
    Adds the `sector` subcommand's own sector-specific options -- on top
    of whatever `add_shared_generation_options` already added -- to
    `parser`: `--name`/`-n` (sector name), `--num-sectors`, and the MySQL
    connection args.

    Args:
        parser (argparse.ArgumentParser): The parser to add options to.
    """
    parser.add_argument('--name', '-n', dest='sector_name', type=str,
                        help="Force the name of the sector, overriding the default random two-word name. "
                             "Cannot be combined with --num-sectors > 1.")
    parser.add_argument('--num-sectors', type=int, default=1,
                        help="Generate this many independent sectors, each with no galactic positioning, "
                             "saving all of them into the same database. Defaults to 1.")
    _db.add_mysql_connection_args(parser)


def validate_sector_args(args, parser):
    """
    Validates `add_sector_arguments`'s own options (`--num-sectors`/
    `--name` conflicts) and shapes `args` into the namespace
    `build_system_config()` expects, calling `parser.error` (which
    exits) on the first problem found. Doesn't include
    `add_shared_generation_options`'s own validation
    (`validate_shared_generation_args`) -- callers run both.

    Args:
        args (argparse.Namespace): Parsed arguments.
        parser (argparse.ArgumentParser): The parser to raise errors
                                          through (so the caller's own
                                          `--help`/usage text is shown).
    """
    if args.num_sectors < 1:
        parser.error("--num-sectors must be a positive integer.")

    if args.num_sectors > 1 and args.sector_name:
        parser.error("--name/-n cannot be combined with --num-sectors > 1 (every generated sector would "
                     "share the same forced name).")

    # build_system_config() expects a namespace shaped like the `system`
    # subcommand's own output, including these three -- deliberately not
    # exposed as sector-level flags (see add_sector_arguments's docstring),
    # so they get the same "not given" default the `system` subcommand's
    # own parser would.
    args.system_file = None
    args.num_orbits = None
    args.name = None


def build_sector_configs(args):
    """
    Builds one `SystemConfig` per system in the sector, sharing the same
    tri-state/value options across all of them (via `build_system_config`,
    so this can never silently drift from what those options mean for a
    single system), then -- if `--min-habitable` was given and not
    already guaranteed by a uniform `+habitable_world` -- forces
    `HABITABLE_WORLD = True` on that many randomly-chosen configs among
    the rest.

    Args:
        args (argparse.Namespace): Parsed arguments, with
            `args.num_systems` already resolved to a concrete count --
            either the user's explicit `--num-systems`, or (see
            `generate_sector`) a per-sector value sampled from `--density`.

    Returns:
        list: A list of `num_systems` `SystemConfig` instances, one per
              system the sector will contain.

    Raises:
        SystemExit: If `--min-habitable` exceeds `args.num_systems` -- always
            caught earlier by `validate_shared_generation_args` for an
            explicit `--num-systems`, but only knowable here for a
            `--density`-driven count, which isn't resolved until generation
            time.
    """
    return list(iter_sector_configs(args))


def iter_sector_configs(args):
    """
    `build_sector_configs`, lazily: yields the same `args.num_systems`
    configs one at a time, so a huge count (`--num-systems 1000000000`,
    or a huge `--density` draw) never builds them all up front --
    `generate_sector` stops pulling once the sector is full. Every config
    comes from the same `args`, so they share one `HABITABLE_WORLD`
    value; the `--min-habitable` indices are therefore chosen (and any
    conflict reported) before the first config is yielded.

    Raises:
        SystemExit: As `build_sector_configs`.
    """
    count = args.num_systems
    if args.min_habitable > count:
        log.error(
            f"Error: --min-habitable ({args.min_habitable}) exceeds this sector's generated system "
            f"count ({count}); with --density, the count is randomly sampled per sector and can "
            f"land below --min-habitable. Try a smaller --min-habitable, a higher --density, or "
            f"--num-systems for an exact count instead."
        )
        raise SystemExit(1)
    if count <= 0:
        return

    first = build_system_config(args)
    forced = set()
    if args.min_habitable > 0 and first.HABITABLE_WORLD is not True:
        forced = set(random.sample(range(count), k=args.min_habitable))
        # Mirrors build_system_config's own habitable-world + asteroid-belt
        # normalization, reapplied here since it ran before this override.
        if first.ASTEROID_BELT is True and first.LARGE_STAR is False:
            log.error(
                "Error: --min-habitable requires forcing a habitable world onto a system that also "
                "has +asteroid_belt forced sector-wide; that combination needs a large star, but "
                "-large_star was also forced sector-wide."
            )
            raise SystemExit(1)

    for i in range(count):
        config = first if i == 0 else build_system_config(args)
        if i in forced:
            config.HABITABLE_WORLD = True
            if config.ASTEROID_BELT is True:
                config.LARGE_STAR = True
        yield config


def sector_star_count(sector):
    """How many stars `sector`'s systems hold (a binary counts two) --
    what every per-star phenomenon rate multiplies."""
    return sum(len(getattr(entry.star_system, "stars", None) or [None]) for entry in sector.entries)


def generate_sector_phenomena(sector, args, galactic_center_dist_ly=None):
    """
    Populates an already-built `sector` with its exotic phenomena, each
    kind in `SECTOR_PHENOMENON_KINDS` sampled independently by a Poisson
    draw (`spaceSector._sample_poisson_count`, the same mechanism
    `--density`'s own system count uses) whose mean is the kind's rate
    per star (`program_constants.phenomenon_rate_per_star`, from Boss's
    research density at the solar neighborhood) times the sector's star
    count. A star count already tracks the local stellar density, so a
    denser sector gets proportionally more. A local-density sector gets
    about 58 rogue planets and 2 brown dwarfs; rarer kinds mostly none.

    Black holes and neutron stars are real, stellar-mass gravitating
    bodies, so they're added via `SpaceSector.add_phenomenon`'s
    Hill-sphere-aware placement (never within a neighboring star system's
    or another compact remnant's own Hill sphere); every other phenomenon
    type has no comparable gravitational footprint at this generator's
    scale and is simply placed at a random point in the sector's cube --
    see that method's own docstring.

    Args:
        sector (SpaceSector): An already-populated sector (every star
            system already added via `add_system`, so their Hill spheres
            are in place for a black hole's/neutron star's own placement
            check) to add phenomena to, in place.
        args (argparse.Namespace): Parsed arguments; only `args.markdown`
            is consulted (threaded into each phenomenon's own
            `SystemConfig`, matching every star system's own config).
        galactic_center_dist_ly (float, optional): As in `generate_sector`
            -- threaded into a standalone black hole's/neutron star's own
            Hill-sphere/galactic-orbit calculation, exactly like every star
            in the sector already gets.

    Returns:
        list: The newly created `SectorPhenomenonEntry` instances. May hold
             fewer than the Poisson draw sampled for a massive type (a
             black hole/neutron star) if the sector was too full/small to
             fit it without violating another massive object's Hill sphere
             -- that one draw is silently skipped (see `ValueError` below)
             rather than crashing the whole sector's generation over a
             single rare phenomenon that didn't fit.
    """
    star_count = sector_star_count(sector)
    new_entries = []

    for kind, phenomenon_type, factory in SECTOR_PHENOMENON_KINDS:
        count = _sample_poisson_count(program_constants.phenomenon_rate_per_star(kind) * star_count)
        for _ in range(count):
            phenomenon_config = SystemConfig()
            phenomenon_config.MARKDOWN = args.markdown
            phenomenon = factory(phenomenon_config, galactic_center_dist_ly)
            if kind == "planetary-nebula":
                entry = _add_planetary_nebula(sector, args, phenomenon, galactic_center_dist_ly)
                if entry is not None:
                    new_entries.append(entry)
                continue
            try:
                new_entries.append(sector.add_phenomenon(phenomenon, phenomenon_type))
            except ValueError:
                # Only a massive type (black-hole/neutron-star) can raise
                # here (see SpaceSector._random_position) -- no room left
                # to place it without violating another massive object's
                # Hill sphere. Skip just this one draw rather than this
                # phenomenon type, or the whole sector.
                continue

    new_entries.extend(add_star_hosted_nebulae(sector, args))
    return new_entries


def _spectral_code(star):
    """A star's spectral code (`"G2V"`), the first word of `Star.type`."""
    return (getattr(star, "type", "") or "").split(" ", 1)[0]


def _add_planetary_nebula(sector, args, nebula, galactic_center_dist_ly=None):
    """
    Adds `nebula` (a planetary nebula) around a new system whose star is
    its hot central white dwarf (`PLANETARY_NEBULA_CENTRAL_STAR_TYPES`),
    placed like any other system. Returns the nebula's entry, or `None`
    when the sector has no room left for the system.
    """
    config = SystemConfig()
    config.MARKDOWN = args.markdown
    config.STAR_TYPE = random.choice(program_constants.PLANETARY_NEBULA_CENTRAL_STAR_TYPES)
    config.BINARY_SYSTEM = False
    system = StarSystem(system_config=config, galactic_center_dist_ly=galactic_center_dist_ly)
    try:
        system_entry = sector.add_system(system, system_config=config)
    except ValueError:
        return None
    log.debug(f"Sector {sector.name!r}: planetary nebula {nebula.name!r} around new central star "
              f"{config.STAR_TYPE} {system.name!r}")
    return sector.add_phenomenon(nebula, "nebula", position=system_entry.position)


def add_star_hosted_nebulae(sector, args):
    """
    Grows emission and reflection nebulae around `sector`'s own hot
    stars (GEN.10): each system whose primary is a main-sequence
    star matching a `NEBULA_HOST_RULES` row rolls that row's chance, and
    on a hit gets a nebula of one of the row's classes centered on it.

    Returns:
        list: The new nebulae's `SectorPhenomenonEntry` instances.
    """
    entries = []
    for system_entry in list(sector.entries):
        stars = getattr(system_entry.star_system, "stars", None) or []
        if not stars:
            continue
        code = _spectral_code(stars[0])
        match = re.fullmatch(r"([OBAFGKM])([0-9])V", code)
        if match is None:
            continue
        letter, subclass = match.group(1), int(match.group(2))
        for rule_letter, low, high, classes, chance in program_constants.NEBULA_HOST_RULES:
            if letter != rule_letter or not (low <= subclass <= high):
                continue
            if random.random() < chance:
                config = SystemConfig()
                config.MARKDOWN = args.markdown
                nebula_class = choose_weighted_class(classes, "Star-hosted nebula class")
                nebula = Nebula(config, nebula_class=nebula_class)
                log.debug(f"Sector {sector.name!r}: {code} star {system_entry.star_system.name!r} "
                          f"lights class {nebula_class} nebula {nebula.name!r}")
                entries.append(sector.add_phenomenon(nebula, "nebula", position=system_entry.position))
            break
    return entries


def _log_uniform(low, high):
    return math.exp(random.uniform(math.log(low), math.log(high)))


def flag_fast_stars(sector, galactic_center_dist_ly=None):
    """
    Marks some of `sector`'s systems as runaway or hypervelocity stars
    (`StarSystem.runaway_class` and `runaway_speed_kms`, schema v37): an
    ordinary generated system moving unusually fast, not a separate
    phenomenon. Each system rolls a hypervelocity chance first -- the
    "hypervelocity-star" rate per star scaled by
    `(HVS_REFERENCE_RADIUS_PC / r)^2`, since the central black hole ejects
    them -- then the runaway chance (about 1.5%).

    Args:
        sector (SpaceSector): The populated sector, changed in place.
        galactic_center_dist_ly (float, optional): The sector's distance
            from the galactic center; `None` uses
            `physical_constants.GALACTIC_CENTER_DISTANCE_LY`.

    Returns:
        int: How many systems were flagged.
    """
    if galactic_center_dist_ly is None:
        galactic_center_dist_ly = physical_constants.GALACTIC_CENTER_DISTANCE_LY
    radius_pc = max(ly_to_pc(galactic_center_dist_ly), 1.0)
    hvs_chance = min(1.0, program_constants.phenomenon_rate_per_star("hypervelocity-star")
                     * (program_constants.HVS_REFERENCE_RADIUS_PC / radius_pc) ** 2)
    runaway_chance = program_constants.phenomenon_rate_per_star("runaway-star")
    flagged = 0
    for entry in sector.entries:
        system = entry.star_system
        if random.random() < hvs_chance:
            system.runaway_class = "hypervelocity"
            system.runaway_speed_kms = _log_uniform(*program_constants.HYPERVELOCITY_STAR_SPEED_RANGE_KMS)
        elif random.random() < runaway_chance:
            system.runaway_class = "runaway"
            system.runaway_speed_kms = _log_uniform(*program_constants.RUNAWAY_STAR_SPEED_RANGE_KMS)
        else:
            continue
        flagged += 1
    return flagged


NUCLEUS_ADDRESS = (0, 0, 0)
"""tuple: The `(ring, layer, slot)` of the one sector per galaxy that rolls
for an active nucleus (see `add_galactic_nucleus`)."""


def add_galactic_nucleus(sector, args, galactic_center_dist_ly):
    """
    Adds the galaxy's central supermassive black hole to `sector` at the
    galactic center: an active `Quasar` with
    `program_constants.QUASAR_ACTIVE_NUCLEUS_CHANCE`, otherwise a quiescent
    supermassive `BlackHole` (like Sagittarius A*), so every galaxy has one.

    There is only ever one nucleus, and only at the origin.
    `generate_and_save_sector_at` calls this for exactly one sector per
    galaxy, `NUCLEUS_ADDRESS` (ring 0, layer 0, slot 0 -- every ring-0,
    layer-0 cell has the galactic axis as its inner edge and the plane
    through its middle, so each touches the origin, and picking one keeps
    it to a single roll).

    The sector's local +X axis points radially out from the galactic axis
    to the sector's own center, and a layer-0 center sits on the plane
    (`galaxyGeometry.sector_orientation`), so the center sits at
    `(-galactic_center_dist_ly, 0, 0)` in the sector's own frame;
    `_db.insert_sector` converts that back to the galaxy origin.

    Args:
        sector (SpaceSector): The core sector, already populated.
        args (argparse.Namespace): Parsed arguments; only `args.markdown`
            is consulted.
        galactic_center_dist_ly (float): The sector center's distance from
            the galactic center, in light-years.

    Returns:
        SectorPhenomenonEntry: The quasar's or the black hole's entry.
    """
    nucleus_config = SystemConfig()
    nucleus_config.MARKDOWN = args.markdown
    position = (-galactic_center_dist_ly, 0.0, 0.0)
    if random.random() < program_constants.QUASAR_ACTIVE_NUCLEUS_CHANCE:
        return sector.add_phenomenon(Quasar(nucleus_config), "quasar", position=position)
    log.debug(f"Sector {sector.name!r}: galactic nucleus is quiescent (supermassive black hole, no quasar)")
    black_hole = BlackHole(nucleus_config, mass_class="supermassive")
    return sector.add_phenomenon(black_hole, "black-hole", position=position)


def _add_preplaced_systems(sector, args, fill, galactic_center_dist_ly):
    """Builds a full system around each of the sector's pre-placed bright
    stars (`fill.bright_rows`) and places it first, at its stored point.
    Each one's companion, planets and moons come from its own stored
    `seed`, without disturbing the run's random state."""
    for row in fill.bright_rows:
        cfg = build_system_config(args)
        cfg.POPULATION = row["population"]
        state = random.getstate()
        random.seed(row["seed"])
        try:
            system = StarSystem(system_config=cfg, galactic_center_dist_ly=galactic_center_dist_ly,
                                primary_star_params=brightStars.star_params(row))
        finally:
            random.setstate(state)
        entry = sector.add_preplaced_system(system, brightStars.local_position_ly(row, fill.center_pc),
                                            system_config=cfg)
        entry.bright_star_id = row["id"]


def generate_sector(args, galactic_center_dist_ly=None, cell=None, fill=None):
    """
    Builds a fully populated `SpaceSector` from parsed args, without
    rendering, printing, or saving anything -- the shared core `run_sector`
    and the `galaxy` subcommand both build on.

    Args:
        args (argparse.Namespace): Parsed arguments from the `sector`
            subcommand (or an equivalently-shaped namespace the `galaxy`
            subcommand builds itself -- see
            `add_shared_generation_options`/`validate_shared_generation_args`
            for the option surface this function actually reads, via
            `build_sector_configs`).
        galactic_center_dist_ly (float, optional): This sector's distance
            from the galactic center, in light-years, threaded down into
            every generated system's Hill-sphere calculation (see
            `Star.calculate_system_perimeter`'s docstring). `None` (the
            default) falls back to the fixed `GALACTIC_CENTER_DISTANCE_LY`
            constant -- the `sector` subcommand (no galaxy context) always
            calls this with the default, so it keeps producing "unplaced"
            sectors.
        cell (galaxyGeometry.SectorCell, optional): The sector's real
            cylindrical cell, in light-years, for a galaxy-placed sector --
            systems and phenomena are then placed inside it rather than in
            a cube (see `SpaceSector.cell`).
        fill (brightStars.FillContext, optional): For a galaxy-placed
            sector: its pre-placed bright stars are built and placed
            first, and every other system takes its age from the stellar
            population mix there and stays below the bright-star
            threshold; a `--density` count shrinks by the bright stars'
            share so the expected total is unchanged.

    Returns:
        tuple: `(sector_name, SpaceSector)` -- `sector_name` is
              `args.sector_name` or a freshly generated one (see
              `generate_sector_name`); the `SpaceSector` has every system
              that fit added (`SpaceSector.add_system`'s Hill-sphere-based
              random placement -- see the "Capacity" note below for what
              happens when not all of them do), plus a realistically
              sparse population of exotic phenomena (see
              `generate_sector_phenomena`). When `args.density` drove the
              system count (see below) and both that Poisson draw and
              `generate_sector_phenomena`'s own independent draws came
              back completely empty, one system is force-added anyway --
              see the "guaranteed non-empty" note below.

    Capacity
    --------
    A sector's cube can only physically hold so many systems before every
    point left is within some existing system's Hill sphere -- real space
    is finite, and `SpaceSector.add_system`'s random placement
    (`_random_position`) raises `ValueError` once it can't find room for
    one more within `program_constants.SECTOR_MAX_PLACEMENT_ATTEMPTS`
    tries. `--num-systems`/`--density` (and, compounding, `--min-habitable`
    forcing extra large stars with their own larger Hill spheres) can ask
    for more systems than a `program_constants.DEFAULT_SECTOR_EDGE_LY`
    cube can hold at realistic stellar spacing. Rather than let that
    `ValueError` propagate (previously discarding every system already
    generated for this sector, including the expensive planet/moon
    generation behind each one) or keep burning full `StarSystem`
    generations on placements that are no longer going to fit, this loop
    stops as soon as one system can't be placed and returns the sector as
    it stands -- with fewer systems than requested, not zero. See
    `sector_generation_summary_lines`'s "actual vs. expected" density
    line, which already accounts for this (it was previously only ever
    exercised by `--density`'s own Poisson undersampling).
    """
    sector_name = args.sector_name or generate_sector_name()
    sector = SpaceSector(name=sector_name, cell=cell)
    if fill is not None:
        _add_preplaced_systems(sector, args, fill, galactic_center_dist_ly)

    # `--density` resolves to a concrete system count per sector (this
    # sector's own volume, sampled fresh each call) rather than once at
    # parse time -- this is what lets each sector under `--num-sectors`/
    # the `galaxy` subcommand vary independently instead of sharing one
    # fixed count. A copy avoids mutating the caller's shared `args`
    # namespace, since this function runs once per sector.
    density_driven = args.density is not None
    if density_driven:
        args = copy.copy(args)
        # An astronomical --density (1e308) can overflow the mean to inf;
        # any count that large just fills the sector (see "Capacity").
        mean = min(sector.expected_system_count() * args.density, sys.float_info.max)
        if fill is not None:
            mean *= 1.0 - fill.bright_share()
        args.num_systems = _sample_poisson_count(mean)

    total = args.num_systems
    for i, cfg in enumerate(iter_sector_configs(args)):
        if fill is not None:
            fill.apply(cfg)
        with log.timed_phase(f"generate system {i + 1}/{total}"):
            system = StarSystem(system_config=cfg, galactic_center_dist_ly=galactic_center_dist_ly)

        try:
            with log.timed_phase(f"place system {i + 1}/{total}"):
                sector.add_system(system, system_config=cfg)
        except ValueError:
            # No room left for another system's Hill sphere in this
            # sector's cube -- see this function's own "Capacity"
            # docstring note. Every following config would almost
            # certainly fail the same way (this sector only gets fuller
            # from here), so stop generating and placing altogether
            # rather than pay for `total - i - 1` more full
            # StarSystem generations just to discard them too.
            log.normal(
                f"Sector '{sector_name}': ran out of room after placing {len(sector.entries)} of "
                f"{total} requested systems -- the {sector.edge_ly:.1f} ly cube has no space left "
                f"that clears every already-placed system's Hill sphere. Returning the sector as-is "
                f"rather than the full requested count."
            )
            break

    with log.timed_phase("generate_sector_phenomena"):
        generate_sector_phenomena(sector, args, galactic_center_dist_ly=galactic_center_dist_ly)
        flag_fast_stars(sector, galactic_center_dist_ly=galactic_center_dist_ly)

    if density_driven and not sector.entries and not sector.phenomena:
        # Guaranteed non-empty: a qualifying sector's own Poisson draws
        # (system count here, each phenomenon type in generate_sector_phenomena)
        # are independent, so all of them landing on zero simultaneously is
        # a real, expected outcome at low means (e.g. ~13% at mean 2) --
        # but a sector with literally nothing in it isn't useful to anyone
        # visiting it, so force exactly one system rather than leave it
        # empty. Only applies when the count came from --density (an
        # explicit --num-systems, including 0, is a deliberate request
        # this never overrides).
        fallback_config = build_system_config(args)
        if fill is not None:
            fill.apply(fallback_config)
        fallback_system = StarSystem(system_config=fallback_config, galactic_center_dist_ly=galactic_center_dist_ly)
        sector.add_system(fallback_system, system_config=fallback_config)

    return sector_name, sector


def sector_generation_summary_lines(sector, args_used):
    """
    Builds the per-sector status lines every sector-generating command
    (`sector`, `galaxy`'s three modes) prints right after a sector is
    saved: how many star systems of each spectral class and phenomena of
    each type were actually generated, plus how the sector's actual star
    count compares to what its own local density predicted -- so a run's
    output is enough on its own to sanity-check what came out of it,
    without a separate database query.

    Args:
        sector (SpaceSector): The freshly generated (and, by the time a
            caller has this, already saved -- so `sector.name` is final)
            sector.
        args_used (argparse.Namespace): The args actually passed to
            `generate_sector` to build this particular sector --
            `.density`/`.num_systems` reflect whichever one drove it
            (for `galaxy` mode, this is `_BatchDensity.resolve`'s own
            per-sector result, not necessarily the top-level parsed args).

    Returns:
        str: Two or three newline-joined, already-indented lines: systems
            by spectral class, phenomena by type (omitted when the sector
            has none), and actual vs. expected star density (`1.0` =
            real local stellar density, matching `--density`'s own
            convention).
    """
    system_types = Counter((entry.star_system.star.type or "?")[0] for entry in sector.entries)
    phenomenon_types = Counter(entry.phenomenon_type for entry in sector.phenomena)

    e_value = sector.expected_system_count()
    if args_used.density is not None:
        expected_density = args_used.density
    elif e_value > 0:
        expected_density = (args_used.num_systems or 0) / e_value
    else:
        expected_density = 0.0
    actual_density = (len(sector.entries) / e_value) if e_value > 0 else 0.0

    lines = []
    if system_types:
        types_str = ", ".join(f"{count} {cls}-type" for cls, count in sorted(system_types.items()))
    else:
        types_str = "none"
    lines.append(f"    Systems: {types_str}")

    if phenomenon_types:
        phenomena_str = ", ".join(
            f"{count} {TYPE_LABELS[ptype]}" for ptype, count in sorted(phenomenon_types.items())
        )
        lines.append(f"    Phenomena: {phenomena_str}")

    lines.append(
        f"    Star density: actual {actual_density:.2f}x local, expected {expected_density:.2f}x local"
    )

    return "\n".join(lines)


def run_sector(args):
    """
    Generates `args.num_sectors` sectors (1 by default), each one a
    fully populated `SpaceSector` (see `generate_sector` -- one
    `SystemConfig` per system via `build_sector_configs`, a full
    `StarSystem` from each, all placed in the sector via
    `SpaceSector.add_system`'s Hill-sphere-based placement) saved into
    the same database. No `galactic_center_dist_ly` is passed to
    `generate_sector` here -- the `sector` subcommand has no galaxy
    context (see the `galaxy` subcommand for that), so every sector it
    produces is "unplaced", regardless of `--num-sectors`.

    Each sector is saved to the database; only a short status line and
    summary per saved sector is printed, so a large `--num-sectors` run
    doesn't flood the console with rendered text nobody asked to see.

    Args:
        args (argparse.Namespace): Validated arguments (`command ==
            "sector"`).
    """
    mysql_config = _db.mysql_config_from_args(args)

    def saved(result, _seconds, _weight):
        if queue.parallel:
            RUN_COUNTS["sectors"] += 1
            RUN_COUNTS["systems"] += result["systems"]
            RUN_COUNTS["phenomena"] += result["phenomena"]
        phenomena_note = f", {result['phenomena']} phenomena" if result["phenomena"] else ""
        log.normal(
            f"Saved sector '{result['name']}' to the database (sector_id={result['sector_id']}, "
            f"{result['systems']} systems{phenomena_note}, "
            f"{mysql_config.database}@{mysql_config.host}:{mysql_config.port})."
        )
        log.normal(result["summary"])

    with _work_queue(args, f"Sectors ({args.num_sectors} unplaced)") as queue:
        for index in range(args.num_sectors):
            queue.submit("sector", f"unplaced-{index}", _unplaced_sector_task, args, on_done=saved)

    if args.num_sectors > 1:
        log.normal(f"Generated {args.num_sectors} sectors.")
    run_population_after(args)


def _unplaced_sector_task(args):
    """One `sector` subcommand sector, generated and saved -- a work queue
    task (see `_fill_sector_task`)."""
    _sector_name, sector = generate_sector(args)
    sector_id = _db.save_sector(sector, config=_db.mysql_config_from_args(args))
    _count_sector(sector)
    return {
        "sector_id": sector_id, "name": sector.name, "systems": len(sector.entries),
        "phenomena": len(sector.phenomena), "summary": sector_generation_summary_lines(sector, args),
    }


# ===========================================================================
# 3. Galaxy generation
# ===========================================================================

LARGE_RING_WARNING_THRESHOLD = 2000
"""int: `--ring I` requires `--limit` or `--yes` when ring `I` holds more
slots than this (see `ring_sector_count`) -- from ring 321, ~4,200 ly
out. Anything larger takes a real, unbounded amount of time and disk, so
it needs an explicit choice."""

def add_galaxy_arguments(parser):
    """
    Adds the `galaxy` subcommand's own mode/placement options -- on top
    of whatever `add_shared_generation_options` already added -- to
    `parser`: the `--ring`/`--block`/`--center-sector` mode-selection
    group, `--layer`, `--slot`, `--block-layer`, `--limit`, `--yes`,
    `--radius-pc`, `--max-ring`,
    `--min-start-density`, and the MySQL connection args.

    Args:
        parser (argparse.ArgumentParser): The parser to add options to.
    """
    mode_group = parser.add_mutually_exclusive_group(required=False)
    mode_group.add_argument('--ring', type=int, metavar='I',
                            help="Batch mode: generate every not-yet-generated sector in ring I (0-indexed "
                                 "cylindrical radius band, one sector edge wide) at --layer J.")
    mode_group.add_argument('--block', metavar='M.I.S.SLAB',
                            help="Block mode: generate every not-yet-generated sector the galaxy's outline "
                                 "allows inside one Galaxy Map drill-down block, given by its key "
                                 "(size.ring.wedge.slab, as in the map's stage links, e.g. 3.40.7.0). "
                                 "Needs --limit or --yes past LARGE_RING_WARNING_THRESHOLD sectors.")
    mode_group.add_argument('--center-sector', type=int, metavar='SECTOR_ID',
                            help="Local-neighborhood mode: generate every not-yet-generated sector "
                                 "within --radius-pc of the given, already galaxy-placed sector's own "
                                 "stored center.")

    parser.add_argument('--layer', type=int, metavar='J',
                        help="With --ring: the height band (signed; 0 is centered on the galactic "
                             "plane). Default: 0.")
    parser.add_argument('--slot', type=int, metavar='K',
                        help="Requires --ring I: single-address mode, generating exactly the one "
                             "not-yet-generated sector at (ring I, layer J, slot K) -- e.g. the address "
                             "behind a designation copied from the interactive Galaxy Map. Skips it "
                             "(exit status 1) if it doesn't qualify (would hold no real content at this "
                             "galaxy's own predicted density) or reports its existing sector_id if it "
                             "was already generated.")
    parser.add_argument('--column', action='store_true',
                        help="With --ring I --slot K: column mode, generating every not-yet-generated "
                             "sector at that ring and slot through every layer the galaxy's stored "
                             "outline reaches there (galaxy_column).")
    parser.add_argument('--shell', action='store_true',
                        help="With --ring I: shell mode, generating every not-yet-generated sector in "
                             "ring I through every layer the outline reaches (a cylindrical shell). Far "
                             "larger than one ring at one layer, so it needs --limit or --yes past "
                             "LARGE_RING_WARNING_THRESHOLD sectors.")
    parser.add_argument('--block-layer', type=int, metavar='J',
                        help="With --block: only the block's sectors on sector layer J (one layer of a "
                             "size-3 block is about 9 sectors).")
    parser.add_argument('--limit', type=int,
                        help="With --ring (or --ring --shell, or --block): generate only the first N "
                             "not-yet-generated slots.")
    parser.add_argument('--yes', action='store_true',
                        help="With --ring (or --ring --shell, or --block): skip the confirmation normally "
                             "required before generating more than LARGE_RING_WARNING_THRESHOLD sectors.")
    parser.add_argument('--radius-pc', type=finite_float,
                        help="With --center-sector: the neighborhood search radius, in parsecs. With "
                             "--ring --slot: after generating that address, also generate its "
                             "neighborhood within this radius. With "
                             "neither --ring nor --center-sector (random-start mode): overrides the "
                             "default 12 pc neighborhood radius around the randomly chosen starting "
                             "sector. Every generated sector also gets the bright stars "
                             "(100 L_sun and up) within 100 ly of it.")
    parser.add_argument('--max-ring', type=int,
                        help="With neither --ring nor --center-sector (random-start mode): the highest "
                             "ring the randomly chosen starting sector may land in. Default: anywhere "
                             "inside the galaxy's stored outline.")
    parser.add_argument('--min-start-density', type=finite_float,
                        help="With neither --ring nor --center-sector (random-start mode): require the "
                             "randomly chosen starting sector's own real relative_density (the same "
                             "'expected' figure printed alongside each saved sector) to be at least this "
                             "value before accepting it -- e.g. 1.0 for at least as dense as the galaxy's "
                             "own real local density. Retried the same way an already-occupied or "
                             "non-qualifying address is (see RANDOM_START_MAX_PLACEMENT_ATTEMPTS). Cannot "
                             "be combined with --density/--num-systems (those override every position's "
                             "density uniformly, leaving no per-position value to compare against).")
    _db.add_mysql_connection_args(parser)


def validate_galaxy_args(args, parser):
    """
    Validates `add_galaxy_arguments`'s own options, calling `parser.error`
    (which exits) on the first problem found, then shapes `args` into the
    namespace `generate_sector()`/`build_sector_configs()` expect.
    Doesn't include `add_shared_generation_options`'s own validation
    (`validate_shared_generation_args`) -- callers run both.

    Args:
        args (argparse.Namespace): Parsed arguments.
        parser (argparse.ArgumentParser): The parser to raise errors
                                          through (so the caller's own
                                          `--help`/usage text is shown).
    """
    block = getattr(args, "block", None)
    block_layer = getattr(args, "block_layer", None)
    if block is not None:
        try:
            args.block = parse_drill_key(block)
        except ValueError as exc:
            parser.error(f"--block: {exc}")
        if args.block.m == 1:
            parser.error("--block takes a block, not a single sector; use --ring --layer --slot.")
        others = [flag for flag, value in (
            ("--layer", args.layer), ("--slot", args.slot), ("--radius-pc", args.radius_pc),
            ("--max-ring", args.max_ring), ("--min-start-density", args.min_start_density),
        ) if value is not None] + [flag for flag, on in (("--column", args.column), ("--shell", args.shell)) if on]
        if others:
            parser.error(f"--block can't be combined with {', '.join(others)}.")
        if block_layer is not None and block_layer not in block_layers(args.block):
            layers = block_layers(args.block)
            parser.error(f"--block-layer must be one of the block's layers, {layers[0]} to {layers[-1]}.")
        if args.limit is not None and not 1 <= args.limit <= generationLimits.MAX_GENERATE_LIMIT:
            parser.error(f"--limit must be between 1 and {generationLimits.MAX_GENERATE_LIMIT}.")
        args.sector_name = args.system_file = args.num_orbits = args.name = None
        return
    if block_layer is not None:
        parser.error("--block-layer requires --block.")

    random_start = args.ring is None and args.center_sector is None

    if args.ring is not None and args.ring < 0:
        parser.error("--ring must be >= 0.")
    if args.ring is not None and args.ring > generationLimits.MAX_GENERATE_RING:
        parser.error(f"--ring must be at most {generationLimits.MAX_GENERATE_RING}.")
    if args.column and (args.ring is None or args.slot is None):
        parser.error("--column requires --ring and --slot.")
    if args.shell and args.ring is None:
        parser.error("--shell requires --ring.")
    if args.column and args.shell:
        parser.error("--column and --shell can't be combined.")
    if (args.column or args.shell) and args.layer is not None:
        parser.error("--column and --shell cover every layer, so they don't take --layer.")
    if args.shell and args.slot is not None:
        parser.error("--shell covers every slot of the ring; use --column for one slot.")
    if args.column and args.slot is not None and args.slot < 0:
        parser.error("--slot must be >= 0.")
    if args.column:
        # Validated as column mode below, not single-address mode.
        column_slot, args.slot = args.slot, None
    if args.layer is not None and args.ring is None:
        parser.error("--layer requires --ring.")
    if args.ring is not None and args.layer is None:
        args.layer = 0

    if args.slot is not None and args.ring is None:
        parser.error("--slot requires --ring.")
    if args.slot is not None and args.slot < 0:
        parser.error("--slot must be >= 0.")
    if args.slot is not None and args.limit is not None:
        parser.error("--limit only applies to --ring (batch mode), not --ring --slot (single-address mode).")
    if args.slot is not None and args.yes:
        parser.error("--yes only applies to --ring (batch mode), not --ring --slot (single-address mode).")
    if args.slot is not None and (args.density is not None or args.num_systems is not None):
        parser.error("--density/--num-systems can't be combined with --ring --slot (single-address "
                     "mode) -- ensure_sector_generated always uses this address's own real predicted "
                     "density, the same as when the galaxy map's own live view found it.")

    if args.center_sector is not None and args.radius_pc is None:
        parser.error("--center-sector requires --radius-pc.")
    if args.radius_pc is not None and args.ring is not None and args.slot is None:
        parser.error("--radius-pc only applies to --center-sector, --ring --slot, or random-start mode "
                     "(neither --ring nor --center-sector), not --ring alone.")
    if args.radius_pc is not None and args.radius_pc <= 0:
        parser.error("--radius-pc must be a positive number.")
    if args.radius_pc is not None and args.radius_pc > generationLimits.MAX_GENERATE_RADIUS_PC:
        parser.error(f"--radius-pc must be at most {generationLimits.MAX_GENERATE_RADIUS_PC:g}.")

    if args.limit is not None and args.ring is None:
        parser.error("--limit only applies to --ring.")
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be a positive integer.")
    if args.limit is not None and args.limit > generationLimits.MAX_GENERATE_LIMIT:
        parser.error(f"--limit must be at most {generationLimits.MAX_GENERATE_LIMIT}.")
    if args.yes and args.ring is None:
        parser.error("--yes only applies to --ring.")
    if args.column:
        if args.limit is not None or args.yes:
            parser.error("--limit and --yes don't apply to --column.")
        args.slot = column_slot

    if args.max_ring is not None and not random_start:
        parser.error("--max-ring only applies to random-start mode (neither --ring nor --center-sector).")
    if args.max_ring is not None and args.max_ring < 0:
        parser.error("--max-ring must be >= 0.")
    if args.max_ring is not None and args.max_ring > generationLimits.MAX_GENERATE_RING:
        parser.error(f"--max-ring must be at most {generationLimits.MAX_GENERATE_RING}.")

    if args.min_start_density is not None and not random_start:
        parser.error("--min-start-density only applies to random-start mode (neither --ring nor "
                     "--center-sector).")
    if args.min_start_density is not None and args.min_start_density <= 0:
        parser.error("--min-start-density must be a positive number.")
    if args.min_start_density is not None and (args.density is not None or args.num_systems is not None):
        parser.error("--min-start-density cannot be combined with --density/--num-systems -- those "
                     "override every position's density uniformly, leaving no per-position value for "
                     "--min-start-density to compare against.")

    # generate_sector()/build_sector_configs() expect a namespace shaped
    # like the `sector` subcommand's own output -- see this function's
    # docstring for why these are always None here.
    args.sector_name = None
    args.system_file = None
    args.num_orbits = None
    args.name = None


def _format_address(address):
    """`(ring, layer, slot)` as the words the CLI prints."""
    ring_index, layer_index, slot_index = address
    return f"ring {ring_index} layer {layer_index} slot {slot_index}"


class _BatchDensity:
    """
    Checks each address against the galaxy's stored outline and resolves
    each sector's own `--density` from the skeleton's real position-based
    `relative_density`, for `galaxy` mode's batch/local-neighborhood/
    random-start generation -- the same mechanism `ensure_sector_generated`
    uses for a single lazily-generated sector, applied across a whole run
    instead of every sector sharing one flat CLI value.

    The outline check (`galaxySkeleton.GalaxyBounds.contains`) always
    applies: nothing is ever generated outside the galaxy, even with an
    explicit `--density`/`--num-systems`. Past it, an explicit flag is a
    deliberate, uniform override for the whole run (`resolve` returns
    `args` unchanged); otherwise `resolve` also applies the exact
    `predicted_star_count >= 1` check, so batch runs never save
    all-but-certainly-empty sectors.

    The skeleton and outline are fetched once, on first use -- neither
    changes mid-run.
    """

    def __init__(self, config):
        self._config = config
        self._skeleton = None
        self._bounds = None

    def _load(self):
        if self._skeleton is None:
            conn = _db.get_connection(self._config)
            try:
                self._skeleton = _db.get_galaxy_shape(conn)
                self._bounds = _db.get_galaxy_bounds(conn)
            finally:
                conn.close()
            if self._skeleton is None:
                raise RuntimeError(
                    "The galaxy's skeleton has never been built (no galaxy_shape row) -- run "
                    "'generate.py plan' first, so every sector can be checked against the galaxy's "
                    "bounds before it is generated."
                )

    @property
    def skeleton(self):
        """The stored `_db.GalaxySkeletonInfo`."""
        self._load()
        return self._skeleton

    @property
    def bounds(self):
        """The stored outline, a `galaxySkeleton.GalaxyBounds`."""
        self._load()
        return self._bounds

    def resolve(self, args, address, position_pc):
        """
        Args:
            args (argparse.Namespace): The `galaxy` subcommand's own parsed
                (and validated) arguments.
            address (tuple): This sector's `(ring, layer, slot)`.
            position_pc (tuple): This sector's `(x, y, z)` center, parsecs.

        Returns:
            argparse.Namespace or None: `None` if this address is outside
                the galaxy's stored outline. Otherwise `args` itself when a
                density/count was given explicitly; else `None` if its own
                exact `predicted_star_count < 1.0`, else a fresh copy with
                `.density` set to this position's own `relative_density`
                and `.num_systems` cleared.
        """
        if not self.bounds.contains(address[0], address[1]):
            return None
        if args.density is not None or args.num_systems is not None:
            return args
        skeleton = self.skeleton

        star_count = predicted_star_count(position_pc, skeleton.shape, skeleton.expected_system_count_at_density_1)
        if star_count < 1.0:
            return None

        resolved = copy.copy(args)
        resolved.density = relative_density(position_pc, skeleton.shape)
        resolved.num_systems = None
        return resolved


def _edge_pc():
    """
    The sector edge length used for every grid computation, in parsecs --
    the one standard, `program_constants.DEFAULT_SECTOR_EDGE_PC` (not a
    CLI option for `galaxy` or `plan`, so the grid and the skeleton always
    agree).
    """
    return float(program_constants.DEFAULT_SECTOR_EDGE_PC)


def _fill_context(args, address, position_pc):
    """The `brightStars.FillContext` for one galaxy sector: its population
    mix, and its unfilled pre-placed bright stars down to its block's
    level (`_db.bright_star_fill_level`: the backfill's, else the galaxy
    scatter's). `None` without a stored skeleton (nothing to take the mix
    from)."""
    conn = _db.get_connection(_db.mysql_config_from_args(args))
    try:
        skeleton = _db.get_galaxy_shape(conn)
        if skeleton is None:
            return None
        level = _db.bright_star_fill_level(conn, *address)
        if level is None:
            return brightStars.FillContext(position_pc, skeleton.shape)
        rows = _db.bright_stars_for_sector(conn, *address)
    finally:
        conn.close()
    return brightStars.FillContext(position_pc, skeleton.shape, rows, min_luminosity_sol=level)


def _block_addresses(block):
    """Every sector `(ring, layer, slot)` of level-3 block `block`."""
    return [(sector.ring, sector.slab, sector.wedge)
            for layer in drill_slabs(block) for sector in drill_block_sectors(block, layer)]


def backfill_bright_stars(config, center_pc, radius_ly=None, min_luminosity_sol=None):
    """
    The bright-star backfill around a generated sector (GEN.23): every
    sector block (`galaxyDrill`'s level-3 blocks, 3x3x3 sectors) with a
    sector within `radius_ly` of `center_pc` gets every star from
    `min_luminosity_sol` up to the level it already holds (its own
    `bright_star_blocks` level, else the galaxy scatter's threshold, else
    no ceiling when no scatter ran), in its sectors that aren't filled
    yet. Each block is drawn whole, under a row lock, and then records
    its new level, so a block already that deep is skipped (one query
    for the whole sphere) and no star is ever drawn twice. A block's draw
    is seeded from the galaxy's scatter seed, the block and the band, so
    it is the same whichever sector reached it first.

    Args:
        config (MySQLConfig): Connection parameters.
        center_pc (tuple): The generated sector's center, galaxy-frame
            parsecs.
        radius_ly (float, optional): Defaults to
            `program_constants.BRIGHT_STAR_BACKFILL_RADIUS_LY` (100 ly).
        min_luminosity_sol (float, optional): Defaults to
            `program_constants.BRIGHT_STAR_BACKFILL_MIN_LUMINOSITY_SOL`
            (100 L_sun).

    Returns:
        dict: `blocks` (drawn now) and `stars` (placed now), both int;
            zeros without a stored skeleton or when the galaxy scatter
            already went that deep.
    """
    radius_pc = ly_to_pc(program_constants.BRIGHT_STAR_BACKFILL_RADIUS_LY if radius_ly is None else radius_ly)
    floor = float(program_constants.BRIGHT_STAR_BACKFILL_MIN_LUMINOSITY_SOL
                  if min_luminosity_sol is None else min_luminosity_sol)
    summary = {"blocks": 0, "stars": 0}
    conn = _db.get_connection(config)
    try:
        skeleton = _db.get_galaxy_shape(conn)
        if skeleton is None:
            return summary
        settings = _db.bright_star_scatter_settings(conn)
        galaxy_level, seed = settings if settings else (None, 0)
        if galaxy_level is not None and galaxy_level <= floor:
            return summary
        bounds = _db.get_galaxy_bounds(conn)
        blocks = {
            drill_parent(DrillBlock(1, ring, slot, layer))
            for ring, layer, slot, *_rest in enumerate_sectors_within_radius(center_pc, radius_pc, skeleton.edge_pc)
            if bounds.contains(ring, layer)
        }
        levels = _db.bright_star_block_levels(conn, blocks)
        todo = sorted(block for block in blocks
                      if levels.get((block.ring, block.wedge, block.slab), math.inf) > floor)
        for block in todo:
            level = _db.lock_bright_star_block(conn, block)
            if level is not None and level <= floor:
                conn.commit()
                continue
            ceiling = level if level is not None else galaxy_level
            addresses = [address for address in _block_addresses(block) if bounds.contains(address[0], address[1])]
            filled = _db.get_occupied_addresses(conn, {address[0] for address in addresses})
            rng = random.Random(f"{seed}:{block.ring}:{block.wedge}:{block.slab}:{floor:g}:{ceiling}")
            rows = list(brightStars.backfill_cells(
                skeleton.shape, [address for address in addresses if address not in filled], skeleton.edge_pc,
                skeleton.expected_system_count_at_density_1, floor, ceiling, rng,
            ))
            _db.insert_bright_stars(conn, rows)
            _db.set_bright_star_block_level(conn, block, floor)
            conn.commit()
            summary["blocks"] += 1
            summary["stars"] += len(rows)
    except BaseException:
        conn.rollback()
        raise
    finally:
        conn.close()
    if summary["blocks"]:
        log.debug(f"bright-star backfill: {summary['stars']} stars in {summary['blocks']} block(s) down to "
                  f"{floor:g} L_sun")
    return summary


def generate_and_save_sector_at(args, address, position_pc, edge_pc):
    """
    Generates one sector via `generate_sector` inside its real grid cell
    and saves it -- the per-sector unit of work every `galaxy` mode
    repeats. First backfills the bright stars within 100 ly of it
    (`backfill_bright_stars`, GEN.23), so this sector's own fill already
    builds its systems around them.

    Args:
        args (argparse.Namespace): Parsed arguments (see `add_galaxy_arguments`).
        address (tuple): This sector's `(ring, layer, slot)`.
        position_pc (tuple): Its `(x, y, z)` center, parsecs.
        edge_pc (float): The sector edge length, parsecs.

    Returns:
        tuple: `(sector_id, sector_name, sector)` of the newly saved sector
            -- `sector_name` is read back after saving, since name
            uniqueness (v22) may rename it on save.
    """
    ring_index, layer_index, slot_index = address
    x, y, z = position_pc
    radius_pc = galactic_radius_pc(position_pc)
    cell = SectorCell.for_ring(ring_index, pc_to_ly(edge_pc))

    backfill_bright_stars(_db.mysql_config_from_args(args), position_pc)
    fill = _fill_context(args, address, position_pc)
    _sector_name, sector = generate_sector(args, galactic_center_dist_ly=pc_to_ly(radius_pc), cell=cell, fill=fill)
    if address == NUCLEUS_ADDRESS:
        add_galactic_nucleus(sector, args, pc_to_ly(radius_pc))

    galaxy_position = {
        "center_x_pc": x, "center_y_pc": y, "center_z_pc": z,
        "galactic_radius_pc": radius_pc,
        "ring_index": ring_index, "layer_index": layer_index, "ring_slot_index": slot_index,
    }
    sector_id = _db.save_sector(sector, config=_db.mysql_config_from_args(args), galaxy_position=galaxy_position)
    _count_sector(sector)
    return sector_id, sector.name, sector


def _default_generation_args(config=None):
    """
    Builds a minimal `argparse.Namespace`, shaped exactly like the
    `sector` subcommand's own parsed output, for a caller with no actual
    command line to parse -- `ensure_sector_generated` in particular,
    which runs at sector-visit time.

    Built by feeding an empty argument list through
    `add_shared_generation_options`/`validate_shared_generation_args`
    (rather than hand-listing every field here) so this can never
    silently drift from whatever options those functions actually
    declare.

    Args:
        config (MySQLConfig, optional): Connection parameters, copied
            onto the returned namespace's `mysql_*` attributes. Defaults
            to `DEFAULT_MYSQL_CONFIG`.

    Returns:
        argparse.Namespace: Every shared generation option at its
            documented default (`num_systems` resolved to `10`) -- a caller
            driving generation by density should set `args.density` and
            clear `args.num_systems` back to `None` first.
    """
    parser = argparse.ArgumentParser(prefix_chars='-+')
    add_shared_generation_options(parser)
    args = parser.parse_args([])
    validate_shared_generation_args(args, parser)

    args.sector_name = None
    args.system_file = None
    args.num_orbits = None
    args.name = None

    config = config or _db.MySQLConfig()
    args.mysql_host = config.host
    args.mysql_port = config.port
    args.mysql_user = config.user
    args.mysql_password = config.password
    args.mysql_database = config.database
    return args


def ensure_sector_generated(ring_index, layer_index, ring_slot_index, config=None):
    """
    The galaxy map's "recalculate on visit" entry point: returns the
    sector already generated at this address if one exists; otherwise
    uses the stored skeleton (see section 4, `build_skeleton`) to decide,
    cheaply and exactly, whether this address is worth generating, and
    -- if so -- generates and saves it on the spot, using its own
    position's `relative_density` as the `--density` multiplier.

    A concurrent visit to the same never-generated address is handled by
    `sectors`'s `UNIQUE (ring_index, layer_index, ring_slot_index)`: the
    losing `INSERT` raises `pymysql.err.IntegrityError`, caught here and
    turned into "return what the other call just created".

    Returns:
        dict: `created` (bool), `qualifies` (bool), `sector_id` (int or
              `None`), `sector_name` (str or `None`, only when `created`).

    Raises:
        ValueError: If `ring_slot_index` is out of range for the ring.
        RuntimeError: If the galaxy's skeleton has never been built --
                     run `generate.py plan` first.
    """
    address = (ring_index, layer_index, ring_slot_index)
    conn = _db.get_connection(config)
    try:
        existing_id = _db.get_sector_id_at(conn, *address)
        if existing_id is not None:
            return {"created": False, "qualifies": True, "sector_id": existing_id, "sector_name": None}

        skeleton = _db.get_galaxy_shape(conn)
        if skeleton is None:
            raise RuntimeError(
                "The galaxy's skeleton has never been built (no galaxy_shape row) -- run "
                "'generate.py plan' first."
            )

        bounds = _db.get_galaxy_bounds(conn)
    finally:
        conn.close()

    position_pc = sector_position_pc(ring_index, layer_index, ring_slot_index, skeleton.edge_pc)
    if not bounds.contains(ring_index, layer_index):
        # Past the layer's stored outer ring -- galaxySkeleton's bound is
        # exact, so this is a certain "no".
        return {"created": False, "qualifies": False, "sector_id": None, "sector_name": None}

    density = relative_density(position_pc, skeleton.shape)
    star_count = predicted_star_count(position_pc, skeleton.shape, skeleton.expected_system_count_at_density_1)
    if star_count < 1.0:
        # Inside the layer (a safe superset) but this slot's own angle
        # didn't clear the exact threshold.
        return {"created": False, "qualifies": False, "sector_id": None, "sector_name": None}

    args = _default_generation_args(config=config)
    args.density = density
    args.num_systems = None

    try:
        sector_id, sector_name, _sector = generate_and_save_sector_at(
            args, address, position_pc, skeleton.edge_pc,
        )
    except pymysql.err.IntegrityError:
        conn = _db.get_connection(config)
        try:
            existing_id = _db.get_sector_id_at(conn, *address)
        finally:
            conn.close()
        if existing_id is None:
            raise
        return {"created": False, "qualifies": True, "sector_id": existing_id, "sector_name": None}

    return {"created": True, "qualifies": True, "sector_id": sector_id, "sector_name": sector_name}


def _log_saved(saved, address, suffix=""):
    designation = provisional_sector_designation(*address)
    log.normal(
        f"Saved sector '{saved['name']}' [{designation}] at {_format_address(address)}{suffix} "
        f"(sector_id={saved['sector_id']})."
    )
    log.normal(saved["summary"])


def _log_level(args):
    """The console severity `main` configured from `--quiet`/`--debug`,
    for the work queue's workers."""
    if getattr(args, "quiet", False):
        return log.SILENT
    if getattr(args, "debug", None) is not None:
        return log.DEBUG
    return log.NORMAL


def _work_queue(args, title):
    """
    The `workQueue.WorkQueue` a run hands its sectors to (PERF.8):
    `--workers` (or `PLANETGEN_WORKERS`) worker processes, by default 80%
    of the cores less one when MySQL runs on this machine, with the
    control database's lease so only one run's workers use the machine
    at a time. One worker generates every sector right here, in order.
    """
    mysql_config = _db.mysql_config_from_args(args)
    workers = workQueue.worker_count(getattr(args, "workers", None), mysql_config.host)
    log.debug(f"{title}: {workers} worker process(es) ({workQueue.cpu_count()} cores).")
    return workQueue.WorkQueue(
        title, workers=workers,
        control_config=_db.control_mysql_config(mysql_config) if workers > 1 else None,
        log_level=_log_level(args), debug_file=getattr(args, "debug", None) or None,
    )


def _fill_sector_task(payload):
    """
    One galaxy sector, start to finish -- a work queue task (PERF.8): runs
    in a worker process (or in this one, with one worker), generates the
    sector in its grid cell and saves it in one transaction, and returns
    what the run reports for it.

    Returns:
        dict: `sector_id`, `name` (as saved), `systems`, `phenomena` and
            `summary` (`sector_generation_summary_lines`).
    """
    sector_args = payload["args"]
    sector_id, sector_name, sector = generate_and_save_sector_at(
        sector_args, payload["address"], payload["position_pc"], payload["edge_pc"],
    )
    return {
        "sector_id": sector_id, "name": sector_name, "systems": len(sector.entries),
        "phenomena": len(sector.phenomena), "summary": sector_generation_summary_lines(sector, sector_args),
    }


def _submit_sector(queue, sector_args, address, position_pc, edge_pc, progress, task, suffix=""):
    """Queues one galaxy sector (`_fill_sector_task`); when it's saved,
    advances `task` and logs it."""
    def saved(result, _seconds, _weight):
        if queue.parallel:
            # A worker's own RUN_COUNTS die with it; the run's are here.
            RUN_COUNTS["sectors"] += 1
            RUN_COUNTS["systems"] += result["systems"]
            RUN_COUNTS["phenomena"] += result["phenomena"]
        progress.update(task, advance=1)
        _log_saved(result, address, suffix=suffix)

    payload = {"args": sector_args, "address": address, "position_pc": position_pc, "edge_pc": edge_pc}
    queue.submit("sector", ",".join(str(part) for part in address), _fill_sector_task, payload, on_done=saved)


def _require_inside(bounds, ring_index, layer_index, what):
    """Stops the run, before anything is generated, when `(ring_index,
    layer_index)` lies outside the galaxy's stored outline."""
    if not bounds.contains(ring_index, layer_index):
        log.error(f"Nothing generated: {what} -- {bounds.describe_miss(ring_index, layer_index)}.")
        raise SystemExit(1)


def run_ring_batch(args, edge_pc, progress):
    """
    Batch mode: generates every not-yet-generated, qualifying sector in
    ring `args.ring` at layer `args.layer` (up to `args.limit`, if given)
    -- when neither `--density` nor `--num-systems` was given, an address
    below the 1-star-per-sector threshold is skipped (see
    `_BatchDensity.resolve`).

    Args:
        args (argparse.Namespace): Parsed arguments; `args.ring` set.
        edge_pc (float): The sector edge length, in parsecs (`_edge_pc`).
        progress (rich.progress.Progress): `run_galaxy`'s shared progress
            display -- an outer "Sectors" task is added here and advanced
            once per slot visited.

    Raises:
        SystemExit: If the ring and layer lie outside the galaxy's stored
                   outline, or the ring's slot count exceeds
                   `LARGE_RING_WARNING_THRESHOLD` and neither `--limit`
                   nor `--yes` was given.
    """
    ring_index, layer_index = args.ring, args.layer
    total_slots = ring_sector_count(ring_index)

    if total_slots > LARGE_RING_WARNING_THRESHOLD and args.limit is None and not args.yes:
        log.error(
            f"Ring {ring_index} holds {total_slots} sector slots -- generating a whole ring this large "
            f"is likely impractical. Pass --limit N to generate only the first N not-yet-generated "
            f"slots, or --yes to confirm generating all {total_slots}."
        )
        raise SystemExit(1)

    mysql_config = _db.mysql_config_from_args(args)
    batch_density = _BatchDensity(mysql_config)
    _require_inside(batch_density.bounds, ring_index, layer_index, f"ring {ring_index} layer {layer_index}")

    conn = _db.get_connection(mysql_config)
    try:
        occupied = {a for a in _db.get_occupied_addresses(conn, [ring_index]) if a[1] == layer_index}
    finally:
        conn.close()

    to_generate = total_slots - len(occupied)
    if args.limit is not None and (args.density is not None or args.num_systems is not None):
        # Only a safe cap when every unoccupied slot is sure to generate
        # (an explicit density bypasses the qualification gate).
        to_generate = min(to_generate, args.limit)
    outer_task = progress.add_task(f"Sectors (ring {ring_index} layer {layer_index})", total=max(to_generate, 0))

    generated = 0
    skipped = 0
    with _work_queue(args, f"Sectors (ring {ring_index} layer {layer_index})") as queue:
        for slot_index in range(total_slots):
            if args.limit is not None and generated >= args.limit:
                break
            address = (ring_index, layer_index, slot_index)
            if address in occupied:
                continue

            position_pc = sector_position_pc(ring_index, layer_index, slot_index, edge_pc)
            sector_args = batch_density.resolve(args, address, position_pc)
            if sector_args is None:
                log.debug(f"{_format_address(address)}: skipped (below the 1-star-per-sector threshold, or "
                          f"outside its layer's stored extent)")
                skipped += 1
                progress.update(outer_task, advance=1)
                continue

            log.debug(f"{_format_address(address)}: generating (density={sector_args.density})")
            _submit_sector(queue, sector_args, address, position_pc, edge_pc, progress, outer_task)
            generated += 1

    skip_note = f", {skipped} skipped (below the star-count threshold)" if skipped else ""
    log.normal(
        f"Generated {generated} new sector(s) in ring {ring_index} layer {layer_index} "
        f"({total_slots} total slots, {len(occupied)} already existed{skip_note})."
    )


def _neighborhood_candidates(center, radius_pc, edge_pc, config, bounds):
    """
    Every address within `radius_pc` of `center` that lies inside the
    galaxy's outline, the set of those already occupied, and how many
    addresses in the sphere were left out for lying outside it -- the
    sphere is trimmed to the galaxy before anything is generated.
    """
    if bounds:
        # Nothing in the galaxy lies farther from `center` than this, so a
        # larger radius only enumerates (and discards) empty space --
        # `--radius-pc 1e300` used to walk ~1e300 rings before trimming.
        top_layer = max(abs(layer) for layer in bounds.outer_ring)
        reach_pc = (math.hypot(center[0], center[1]) + abs(center[2])
                    + (bounds.outer_ring_index + top_layer + 2) * edge_pc)
        radius_pc = min(radius_pc, reach_pc)
    candidates = []
    outside = 0
    for candidate in enumerate_sectors_within_radius(center, radius_pc, edge_pc):
        if bounds.contains(candidate[0], candidate[1]):
            candidates.append(candidate)
        else:
            outside += 1
    conn = _db.get_connection(config)
    try:
        occupied = _db.get_occupied_addresses(conn, {c[0] for c in candidates})
    finally:
        conn.close()
    return candidates, occupied, outside


def run_local_neighborhood(args, edge_pc, progress):
    """
    Local-neighborhood mode: generates every not-yet-generated, qualifying
    sector within `args.radius_pc` parsecs of `args.center_sector`'s own
    stored center -- see `run_ring_batch` for what "qualifying" means.

    Args:
        args (argparse.Namespace): Parsed arguments;
            `args.center_sector`/`args.radius_pc` must not be `None`.
        edge_pc (float): The sector edge length, in parsecs (`_edge_pc`).
        progress (rich.progress.Progress): `run_galaxy`'s shared progress
            display (see `run_ring_batch`). `run_random_start` also calls
            this directly after generating its own seed sector.

    Raises:
        SystemExit: If `args.center_sector` doesn't exist, or has never
                   been placed in a galaxy.
    """
    conn = _db.get_connection(_db.mysql_config_from_args(args))
    try:
        try:
            center_position = _db.get_sector_galaxy_position(conn, args.center_sector)
        except ValueError as exc:
            log.error(str(exc))
            raise SystemExit(1) from exc
    finally:
        conn.close()

    if center_position is None:
        log.error(
            f"sector_id={args.center_sector} has never been placed in a galaxy (its galaxy-position "
            f"columns are NULL) -- --center-sector requires an already galaxy-placed sector (one "
            f"generated via 'generate.py galaxy' itself, not 'generate.py sector'). Use --ring to "
            f"generate placed sectors from scratch instead."
        )
        raise SystemExit(1)

    center = (
        center_position["center_x_pc"], center_position["center_y_pc"], center_position["center_z_pc"],
    )
    mysql_config = _db.mysql_config_from_args(args)
    batch_density = _BatchDensity(mysql_config)
    center_ring, center_layer, _slot = sector_address_at(center, edge_pc)
    _require_inside(batch_density.bounds, center_ring, center_layer,
                    f"sector_id={args.center_sector} sits at {_format_address(sector_address_at(center, edge_pc))}")
    candidates, occupied, outside = _neighborhood_candidates(
        center, args.radius_pc, edge_pc, mysql_config, batch_density.bounds,
    )

    to_generate = sum(1 for c in candidates if c[:3] not in occupied)
    outer_task = progress.add_task("Sectors (local neighborhood)", total=to_generate)

    generated = 0
    skipped = 0
    already_existed = 0
    with _work_queue(args, f"Sectors (within {args.radius_pc:g} pc of sector {args.center_sector})") as queue:
        for ring_index, layer_index, slot_index, x, y, z, distance_pc in candidates:
            address = (ring_index, layer_index, slot_index)
            if address in occupied:
                already_existed += 1
                continue

            sector_args = batch_density.resolve(args, address, (x, y, z))
            if sector_args is None:
                log.debug(f"{_format_address(address)}: skipped (below the 1-star-per-sector threshold, or "
                          f"outside its layer's stored extent)")
                skipped += 1
                progress.update(outer_task, advance=1)
                continue

            log.debug(f"{_format_address(address)}: generating (density={sector_args.density})")
            _submit_sector(queue, sector_args, address, (x, y, z), edge_pc, progress, outer_task,
                           suffix=f", {distance_pc:.2f} pc from sector_id={args.center_sector}")
            generated += 1

    skip_note = f", {skipped} skipped (below the star-count threshold)" if skipped else ""
    outside_note = f", {outside} beyond the galaxy's edge left out" if outside else ""
    log.normal(
        f"Generated {generated} new sector(s) within {args.radius_pc} pc of sector_id={args.center_sector} "
        f"({len(candidates)} candidate slot(s) inside the galaxy{outside_note}, {already_existed} already "
        f"existed{skip_note})."
    )


def generate_sector_neighborhood(center_sector_id, radius_ly=None, config=None):
    """
    Non-CLI counterpart to `run_local_neighborhood` -- for the admin web
    UI's "generate more sectors around this one" action
    (`html/api/routes.py`'s `generate_sector_neighborhood_route`). Same
    work, a plain result dict instead of prints, and a catchable
    `ValueError` instead of `SystemExit` for an invalid/unplaced sector.

    The default radius is `program_constants.DEFAULT_GENERATE_RADIUS_PC`
    (12 pc, about 100 candidate addresses; GEN.23). The 100 ly around each
    generated sector gets only its bright stars (`backfill_bright_stars`).
    A larger `radius_ly` (100 ly holds roughly 1,500-2,000 addresses) can
    take minutes to hours.

    Args:
        center_sector_id (int): The already galaxy-placed sector to
            generate a neighborhood around.
        radius_ly (float, optional): Defaults to
            `program_constants.DEFAULT_GENERATE_RADIUS_PC` (12 pc).
        config (MySQLConfig, optional): Connection parameters.

    Returns:
        dict: `generated`, `already_existed`, `skipped`, `candidates`
            (addresses inside the galaxy), `outside_galaxy` (addresses in
            the sphere past the galaxy's edge, left out) -- all int.

    Raises:
        ValueError: If `center_sector_id` doesn't exist, has never been
                   placed in a galaxy, or lies outside its outline.
        RuntimeError: If the galaxy's skeleton has never been built.
    """
    edge_pc = _edge_pc()
    radius_pc = (
        ly_to_pc(radius_ly) if radius_ly is not None
        else program_constants.DEFAULT_GENERATE_RADIUS_PC
    )
    config = config or _db.DEFAULT_MYSQL_CONFIG

    conn = _db.get_connection(config)
    try:
        center_position = _db.get_sector_galaxy_position(conn, center_sector_id)
    finally:
        conn.close()

    if center_position is None:
        raise ValueError(
            f"sector_id={center_sector_id} has never been placed in a galaxy (its galaxy-position "
            f"columns are NULL) -- generating a neighborhood requires an already galaxy-placed sector "
            f"(one generated via 'generate.py galaxy', not 'generate.py sector')."
        )

    center = (center_position["center_x_pc"], center_position["center_y_pc"], center_position["center_z_pc"])
    batch_density = _BatchDensity(config)
    center_ring, center_layer, _slot = sector_address_at(center, edge_pc)
    if not batch_density.bounds.contains(center_ring, center_layer):
        raise ValueError(
            f"sector_id={center_sector_id} lies outside the galaxy: "
            f"{batch_density.bounds.describe_miss(center_ring, center_layer)}."
        )
    candidates, occupied, outside = _neighborhood_candidates(
        center, radius_pc, edge_pc, config, batch_density.bounds,
    )

    args = _default_generation_args(config=config)
    # Density-driven from the skeleton (_BatchDensity), not the flat
    # num_systems=10 default.
    args.density = None
    args.num_systems = None

    generated = 0
    skipped = 0
    already_existed = 0
    for ring_index, layer_index, slot_index, x, y, z, _distance_pc in candidates:
        address = (ring_index, layer_index, slot_index)
        if address in occupied:
            already_existed += 1
            continue
        sector_args = batch_density.resolve(args, address, (x, y, z))
        if sector_args is None:
            skipped += 1
            continue
        generate_and_save_sector_at(sector_args, address, (x, y, z), edge_pc)
        generated += 1

    return {
        "generated": generated,
        "already_existed": already_existed,
        "skipped": skipped,
        "candidates": len(candidates),
        "outside_galaxy": outside,
    }


def run_random_start(args, edge_pc, progress):
    """
    Random-start mode (no `--ring`/`--center-sector` given): picks a
    random, not-yet-occupied, qualifying sector address -- drawn only from
    inside the galaxy's stored outline (`GalaxyBounds.random_address`,
    every sector equally likely, optionally only out to `--max-ring`), so
    the start and the neighborhood around it are always in the galaxy --
    retried up to `program_constants.RANDOM_START_MAX_PLACEMENT_ATTEMPTS`
    times, generates it, then generates every not-yet-generated sector
    within `args.radius_pc` of it (default
    `program_constants.DEFAULT_GENERATE_RADIUS_PC`, 12 pc) via
    `run_local_neighborhood`. `--min-start-density` tightens the retry:
    a qualifying address below it is retried too.

    Raises:
        SystemExit: If no suitable address was found within the attempt
                   budget.
    """
    radius_pc = (
        args.radius_pc if args.radius_pc is not None
        else program_constants.DEFAULT_GENERATE_RADIUS_PC
    )

    mysql_config = _db.mysql_config_from_args(args)
    batch_density = _BatchDensity(mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        sector_args = None
        for _ in range(program_constants.RANDOM_START_MAX_PLACEMENT_ATTEMPTS):
            address = batch_density.bounds.random_address(random, max_ring=args.max_ring)
            if _db.get_sector_id_at(conn, *address) is not None:
                continue
            position_pc = sector_position_pc(*address, edge_pc)
            sector_args = batch_density.resolve(args, address, position_pc)
            if sector_args is None:
                continue
            if args.min_start_density is not None and sector_args.density < args.min_start_density:
                continue
            break
        else:
            density_note = (
                f", meeting --min-start-density {args.min_start_density} (try lowering it or --max-ring "
                f"-- most volume-weighted draws land in the galaxy's own sparser outskirts)"
                if args.min_start_density is not None else ""
            )
            within = f"within {args.max_ring} rings" if args.max_ring is not None else "in the galaxy"
            log.error(
                f"Could not find an unoccupied, qualifying sector address {within} after {program_constants.RANDOM_START_MAX_PLACEMENT_ATTEMPTS} attempts{density_note} "
                f"-- this galaxy may already be almost entirely generated within that range, or that range "
                f"may hold too little real stellar density; try a larger --max-ring."
            )
            raise SystemExit(1)
    finally:
        conn.close()

    sector_id, sector_name, sector = generate_and_save_sector_at(sector_args, address, position_pc, edge_pc)
    designation = provisional_sector_designation(*address)
    log.normal(
        f"Saved random starting sector '{sector_name}' [{designation}] at {_format_address(address)} "
        f"(sector_id={sector_id})."
    )
    log.normal(sector_generation_summary_lines(sector, sector_args))

    args.center_sector = sector_id
    args.radius_pc = radius_pc
    run_local_neighborhood(args, edge_pc, progress)


def run_single_slot(args, edge_pc, progress):
    """
    Single-address mode: generates exactly the one sector at `(--ring I,
    --layer J, --slot K)` via the same `ensure_sector_generated` the
    galaxy map's "recalculate on visit" path uses -- the direct route from
    an address copied out of the interactive 3D Galaxy Map. With
    `--radius-pc`, then generates that sector's neighborhood
    (`run_local_neighborhood`), the Sector Map's "Generate neighborhood".

    Raises:
        SystemExit: If `--slot` is out of range for the ring, or the
                   address doesn't qualify.
    """
    address = (args.ring, args.layer, args.slot)
    total_slots = ring_sector_count(args.ring)
    if not (0 <= args.slot < total_slots):
        log.error(
            f"--slot {args.slot} is out of range for ring {args.ring} (holds {total_slots} slots, "
            f"0..{total_slots - 1})."
        )
        raise SystemExit(1)

    mysql_config = _db.mysql_config_from_args(args)
    _require_inside(_BatchDensity(mysql_config).bounds, args.ring, args.layer, _format_address(address))
    task = progress.add_task(f"Sector ({_format_address(address)})", total=1)
    result = ensure_sector_generated(*address, config=mysql_config)
    progress.update(task, advance=1)

    if not result["qualifies"]:
        log.error(
            f"{_format_address(address)} doesn't qualify -- it would hold no real content at this "
            f"galaxy's own predicted density (below the 1-star-per-sector threshold, or outside the "
            f"layer's stored extent)."
        )
        raise SystemExit(1)

    designation = provisional_sector_designation(*address)
    if result["created"]:
        log.normal(
            f"Saved sector '{result['sector_name']}' [{designation}] at {_format_address(address)} "
            f"(sector_id={result['sector_id']})."
        )
    else:
        log.normal(
            f"Sector [{designation}] at {_format_address(address)} already existed "
            f"(sector_id={result['sector_id']})."
        )

    if args.radius_pc is not None:
        # Then its neighborhood, the same way random-start mode follows
        # its seed sector.
        args.center_sector = result["sector_id"]
        run_local_neighborhood(args, edge_pc, progress)


def _layers_reaching(bounds, ring_index):
    """Every layer the galaxy's stored outline reaches at `ring_index`,
    lowest first (the column `galaxy_column` stores for that ring)."""
    return sorted(layer for layer in bounds.outer_ring if bounds.contains(ring_index, layer))


def _generate_addresses(args, addresses, what, edge_pc, progress, batch_density):
    """
    Generates every not-yet-generated, qualifying address of `addresses`
    in order (up to `args.limit`, when set) -- the loop behind column and
    shell modes; see `run_ring_batch` for what "qualifying" means.
    """
    conn = _db.get_connection(_db.mysql_config_from_args(args))
    try:
        occupied = set(_db.get_occupied_addresses(conn, {a[0] for a in addresses}))
    finally:
        conn.close()

    pending = [a for a in addresses if a not in occupied]
    limit = getattr(args, "limit", None)
    task = progress.add_task(f"Sectors ({what})", total=len(pending))
    generated = skipped = 0
    with _work_queue(args, f"Sectors ({what})") as queue:
        for address in pending:
            if limit is not None and generated >= limit:
                break
            position_pc = sector_position_pc(*address, edge_pc)
            sector_args = batch_density.resolve(args, address, position_pc)
            if sector_args is None:
                skipped += 1
                progress.update(task, advance=1)
                continue
            _submit_sector(queue, sector_args, address, position_pc, edge_pc, progress, task)
            generated += 1

    skip_note = f", {skipped} skipped (below the star-count threshold)" if skipped else ""
    log.normal(
        f"Generated {generated} new sector(s) in {what} ({len(addresses)} total, "
        f"{len(addresses) - len(pending)} already existed{skip_note})."
    )


def block_layers(block):
    """
    Every sector layer a drill-down block spans, lowest first: a size-3
    block's three (`drill_slabs`), and for a bigger one the layers of
    each of its child slabs in turn.
    """
    m, _ring, _wedge, slab = block
    if m == 3:
        return drill_slabs(block)
    child_m = DRILL_LEVELS[DRILL_LEVELS.index(m) + 1]
    layers = []
    for child_slab in drill_slabs(block):
        layers.extend(block_layers(type(block)(child_m, 0, 0, child_slab)))
    return layers


def block_addresses(block, layer=None):
    """
    Every sector address `(ring, layer, slot)` inside drill-down block
    `block` (only on sector layer `layer` when given), ring by ring --
    a size-3 block's sectors (`drill_block_sectors`, design doc section
    3.5), and a bigger block's children's, recursively.
    """
    if block.m == 3:
        layers = drill_slabs(block) if layer is None else [layer]
        for sector_layer in layers:
            for sector in drill_block_sectors(block, sector_layer):
                yield (sector.ring, sector_layer, sector.wedge)
        return
    for _slab, children in drill_children(block):
        for child in children:
            if layer is not None and layer not in block_layers(child):
                continue
            yield from block_addresses(child, layer)


def run_block(args, edge_pc, progress):
    """
    Block mode: every sector the galaxy's outline allows inside one
    drill-down block (`--block`, optionally one `--block-layer`). Needs
    `--limit` or `--yes` when that is more than
    `LARGE_RING_WARNING_THRESHOLD` sectors.

    Raises:
        SystemExit: If the outline allows nothing in the block, or it is
                   too large and neither `--limit` nor `--yes` was given.
    """
    batch_density = _BatchDensity(_db.mysql_config_from_args(args))
    bounds = batch_density.bounds
    key = format_drill_key(args.block)
    where = f"block {key}" + (f" layer {args.block_layer}" if args.block_layer is not None else "")
    confirmed = args.limit is not None or args.yes
    addresses = []
    for address in block_addresses(args.block, args.block_layer):
        if not bounds.contains(address[0], address[1]):
            continue
        addresses.append(address)
        if not confirmed and len(addresses) > LARGE_RING_WARNING_THRESHOLD:
            # A big block can hold millions of sectors: stop counting here.
            log.error(
                f"{where.capitalize()} holds more than {LARGE_RING_WARNING_THRESHOLD} sector slots -- pass "
                f"--limit N to generate only the first N, or --yes to confirm generating all of them."
            )
            raise SystemExit(1)
    if not addresses:
        log.error(f"The galaxy's outline allows no sector in {where}.")
        raise SystemExit(1)
    _generate_addresses(args, addresses, where, edge_pc, progress, batch_density)


def run_column(args, edge_pc, progress):
    """
    Column mode: every sector at `(--ring I, --slot K)` through every
    layer the galaxy's outline reaches at ring I.

    Raises:
        SystemExit: If `--slot` is out of range for the ring, or the
                   outline doesn't reach the ring at all.
    """
    total_slots = ring_sector_count(args.ring)
    if not (0 <= args.slot < total_slots):
        log.error(
            f"--slot {args.slot} is out of range for ring {args.ring} (holds {total_slots} slots, "
            f"0..{total_slots - 1})."
        )
        raise SystemExit(1)
    batch_density = _BatchDensity(_db.mysql_config_from_args(args))
    layers = _layers_reaching(batch_density.bounds, args.ring)
    if not layers:
        _require_inside(batch_density.bounds, args.ring, 0, f"ring {args.ring}")
    addresses = [(args.ring, layer, args.slot) for layer in layers]
    _generate_addresses(args, addresses, f"the column at ring {args.ring} slot {args.slot}",
                        edge_pc, progress, batch_density)


def run_shell(args, edge_pc, progress):
    """
    Shell mode: every sector of ring `--ring I` through every layer the
    galaxy's outline reaches there. Needs `--limit` or `--yes` when that
    is more than `LARGE_RING_WARNING_THRESHOLD` sectors.

    Raises:
        SystemExit: If the outline doesn't reach the ring, or the shell is
                   too large and neither `--limit` nor `--yes` was given.
    """
    batch_density = _BatchDensity(_db.mysql_config_from_args(args))
    layers = _layers_reaching(batch_density.bounds, args.ring)
    if not layers:
        _require_inside(batch_density.bounds, args.ring, 0, f"ring {args.ring}")
    total_slots = ring_sector_count(args.ring)
    total = total_slots * len(layers)
    if total > LARGE_RING_WARNING_THRESHOLD and args.limit is None and not args.yes:
        log.error(
            f"The shell at ring {args.ring} holds {total} sector slots ({total_slots} slots across "
            f"{len(layers)} layers) -- pass --limit N to generate only the first N, or --yes to "
            f"confirm generating all of them."
        )
        raise SystemExit(1)
    addresses = [(args.ring, layer, slot) for layer in layers for slot in range(total_slots)]
    _generate_addresses(args, addresses, f"the shell at ring {args.ring}", edge_pc, progress, batch_density)


def run_galaxy(args):
    """
    Dispatches to block, column, shell, single-address, ring-batch,
    local-neighborhood, or random-start mode, owning the one `rich.progress.Progress` display
    they share. First checks the galaxy has been planned at the standard
    sector edge, since every mode validates its addresses against that
    plan's stored outline before generating anything.

    Args:
        args (argparse.Namespace): Validated arguments (`command ==
            "galaxy"`).
    """
    edge_pc = _edge_pc()
    conn = _db.get_connection(_db.mysql_config_from_args(args))
    try:
        bounds = _db.get_galaxy_bounds(conn)
    finally:
        conn.close()
    if bounds is None:
        log.error("The galaxy has never been planned -- run 'generate.py plan' first, so every sector "
                  "can be checked against the galaxy's bounds before it is generated.")
        raise SystemExit(1)
    if not math.isclose(bounds.edge_pc, edge_pc):
        log.error(f"The stored plan was built with {bounds.edge_pc:g} pc sectors, not the standard "
                  f"{edge_pc:g} pc -- re-run 'generate.py plan'.")
        raise SystemExit(1)
    if not bounds:
        log.error("The stored plan has no layers at all (nothing in this galaxy shape expects a star "
                  "per sector) -- re-run 'generate.py plan' with a different shape.")
        raise SystemExit(1)

    with _generation_progress() as progress:
        log.set_console(progress.console)
        try:
            if getattr(args, "block", None) is not None:
                run_block(args, edge_pc, progress)
            elif args.column:
                run_column(args, edge_pc, progress)
            elif args.shell:
                run_shell(args, edge_pc, progress)
            elif args.slot is not None:
                run_single_slot(args, edge_pc, progress)
            elif args.ring is not None:
                run_ring_batch(args, edge_pc, progress)
            elif args.center_sector is not None:
                run_local_neighborhood(args, edge_pc, progress)
            else:
                run_random_start(args, edge_pc, progress)
        finally:
            log.reset_console()
    run_population_after(args)


# ===========================================================================
# 4. Galaxy density skeleton
# ===========================================================================

def add_plan_arguments(parser):
    """
    Adds every option the `plan` subcommand accepts (besides `--version`)
    to `parser` -- the galaxy shape parameter group plus the scan options
    and the MySQL connection args.

    Args:
        parser (argparse.ArgumentParser): The parser to add options to.
    """
    shape_group = parser.add_argument_group("galaxy shape (galaxyDensity.GalaxyShape)")
    shape_group.add_argument('--disk-scale-length-pc', type=finite_float, default=2800.0,
                             help="Exponential disk radial scale length, parsecs. Default: 2800 "
                                  "(real Milky Way scale).")
    shape_group.add_argument('--disk-scale-height-pc', type=finite_float, default=350.0,
                             help="Disk vertical scale height, parsecs. Default: 350.")
    shape_group.add_argument('--bulge-scale-radius-pc', type=finite_float, default=200.0,
                             help="Bulge exponential scale radius, parsecs. Default: 200.")
    shape_group.add_argument('--bulge-amplitude', type=finite_float, default=1.0,
                             help="Bulge amplitude, relative to the disk term. Default: 1.0.")
    shape_group.add_argument('--arm-count', type=int, default=2,
                             help="Number of spiral arms. Default: 2 (grand-design).")
    shape_group.add_argument('--pitch-angle-deg', type=finite_float, default=15.0,
                             help="Spiral arm pitch angle, degrees. Default: 15.")
    shape_group.add_argument('--arm-amplitude', type=finite_float, default=0.4,
                             help="Arm/inter-arm density contrast amplitude, in [0, 1). Default: 0.4.")
    shape_group.add_argument('--calibration-radius-pc', type=finite_float, default=None,
                             help="In-plane radius the relative_density=1.0 calibration point sits at. "
                                  "Defaults to build_galaxy_shape's own default (2.82x disk scale length).")

    parser.add_argument('--max-ring', type=int, default=DEFAULT_MAX_RING,
                        help=f"Hard cap on how far out a layer is scanned. Default: {DEFAULT_MAX_RING}.")

    bright_group = parser.add_argument_group("bright-star pre-placement (stellarObjects.brightStars)")
    bright_group.add_argument('--bright-star-min-luminosity', type=finite_float,
                              default=program_constants.BRIGHT_STAR_MIN_LUMINOSITY_SOL,
                              help="Every star at least this bright (solar luminosities) is generated "
                                   "and placed galaxy-wide after the plan. Default: "
                                   f"{program_constants.BRIGHT_STAR_MIN_LUMINOSITY_SOL:g}.")
    bright_group.add_argument('--no-bright-stars', action='store_true',
                              help="Build the plan without scattering bright stars.")
    bright_group.add_argument('--bright-stars-only', action='store_true',
                              help="Re-scatter the bright stars on the stored plan without rebuilding it.")
    bright_group.add_argument('--force', action='store_true',
                              help="Scatter even when sectors are already filled, leaving those sectors out.")
    _db.add_mysql_connection_args(parser)
    add_logging_arguments(parser)


def validate_plan_args(args, parser):
    """
    Validates `add_plan_arguments`'s own options, calling `parser.error`
    (which exits) on the first problem found.

    Args:
        args (argparse.Namespace): Parsed arguments.
        parser (argparse.ArgumentParser): The parser to raise errors through.
    """
    if args.arm_amplitude < 0 or args.arm_amplitude >= 1:
        parser.error("--arm-amplitude must be in [0, 1).")
    if args.max_ring < 1:
        parser.error("--max-ring must be a positive integer.")
    for option in ("disk_scale_length_pc", "disk_scale_height_pc", "bulge_scale_radius_pc"):
        if getattr(args, option) <= 0:
            parser.error(f"--{option.replace('_', '-')} must be a positive number.")
    if args.arm_count < 1:
        parser.error("--arm-count must be a positive integer.")
    if not 0 < abs(args.pitch_angle_deg) <= 90:
        parser.error("--pitch-angle-deg must be non-zero and at most 90 degrees in magnitude.")
    if args.calibration_radius_pc is not None and args.calibration_radius_pc <= 0:
        parser.error("--calibration-radius-pc must be a positive number.")
    # Anything else that still can't be normalized (e.g. a calibration
    # radius so far out the density there underflows to zero).
    try:
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
    except (ArithmeticError, ValueError) as exc:
        parser.error(f"these galaxy shape parameters can't be normalized ({exc}).")
    if not (math.isfinite(shape.k_norm) and shape.k_norm > 0):
        parser.error(f"these galaxy shape parameters can't be normalized (k_norm={shape.k_norm!r}).")
    if args.no_bright_stars and args.bright_stars_only:
        parser.error("--no-bright-stars and --bright-stars-only can't be combined.")
    if not args.no_bright_stars:
        try:
            bright_star_fraction(args.bright_star_min_luminosity)
        except ValueError as exc:
            parser.error(f"--bright-star-min-luminosity: {exc}")


def build_skeleton(args):
    """
    Builds and persists the galaxy's skeleton: `galaxy_shape` (the shape
    parameters, calibration constant, edge length and outer ring -- one
    singleton row) and `galaxy_layer` (the galaxy's outline: one row per
    layer, highest to lowest, holding the last ring that layer reaches --
    see `galaxySkeleton.build_layer_extents`), plus `galaxy_column` (the
    same outline per ring: the highest and lowest layer each ring
    reaches). Together they are the bounds every generation path checks
    first. The edge is always the standard
    `program_constants.DEFAULT_SECTOR_EDGE_PC`.

    No sector content or individual addresses are stored -- that stays
    lazy (`ensure_sector_generated`). Re-running replaces the whole
    skeleton; a full Milky-Way-scale build takes a few milliseconds.

    Args:
        args (argparse.Namespace): Parsed arguments.

    Returns:
        dict: `outer_ring_index`, `layer_count`, `top_layer_index`,
              `total_candidate_sectors`, `elapsed_s`, `edge_confirmed`.
    """
    edge_pc = _edge_pc()
    e_value = expected_system_count_at_density_1(program_constants.DEFAULT_SECTOR_EDGE_LY)
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

    log.normal(
        f"Building skeleton: edge_pc={edge_pc:g} disk_scale_length_pc={shape.disk_scale_length_pc} "
        f"disk_scale_height_pc={shape.disk_scale_height_pc} "
        f"bulge_scale_radius_pc={shape.bulge_scale_radius_pc} "
        f"bulge_amplitude={shape.bulge_amplitude} arm_count={shape.arm_count} "
        f"k_norm={shape.k_norm:.4f} threshold_rho={threshold_rho:.6f}"
    )

    t0 = time.perf_counter()
    extents, outer_ring_index, edge_confirmed = build_layer_extents(
        shape, edge_pc, threshold_rho, max_ring=args.max_ring,
    )
    elapsed = time.perf_counter() - t0

    mysql_config = _db.mysql_config_from_args(args)
    # A new outline invalidates every pre-placed bright star (schema v43).
    conn = _db.get_connection(mysql_config)
    try:
        _db.clear_bright_stars(conn)
    finally:
        conn.close()
    _db.replace_galaxy_layers(extents, config=mysql_config)
    _db.save_galaxy_shape(
        shape, edge_pc=edge_pc, outer_ring_index=outer_ring_index,
        expected_system_count_at_density_1=e_value, config=mysql_config,
    )

    return {
        "outer_ring_index": outer_ring_index,
        "layer_count": len(extents),
        "top_layer_index": extents[0][0] if extents else None,
        "total_candidate_sectors": candidate_sector_count(extents),
        "elapsed_s": elapsed,
        "edge_confirmed": edge_confirmed,
    }


def scatter_bright_stars(args):
    """
    Pre-places every bright star on the stored plan (`brightStars.scatter`)
    into `bright_stars`, replacing any earlier scatter, one layer per
    commit, with a progress bar.

    Refuses when sectors are already filled (they'd never get their bright
    stars) unless `--force`, which leaves those sectors out.

    Returns:
        dict or None: `counts` (per population), `total` and `elapsed_s`;
            `None` when it refused.
    """
    mysql_config = _db.mysql_config_from_args(args)
    conn = _db.get_connection(mysql_config)
    try:
        skeleton = _db.get_galaxy_shape(conn)
        if skeleton is None:
            raise RuntimeError("The galaxy's skeleton has never been built -- run 'generate.py plan' first.")
        extents = _db.get_galaxy_layers(conn)
        filled = _db.filled_sector_addresses(conn)
        if filled and not args.force:
            log.normal(
                f"Not scattering bright stars: {len(filled):,} sectors are already filled and would never "
                f"get theirs. Reset the galaxy first, or re-run with --force to leave those sectors out."
            )
            return None
        min_luminosity_sol = float(args.bright_star_min_luminosity)
        seed = random.SystemRandom().getrandbits(63)
        _db.clear_bright_stars(conn)

        counts = {population: 0 for population in brightStars.POPULATIONS}
        t0 = time.perf_counter()
        with _generation_progress() as progress:
            log.set_console(progress.console)
            try:
                task = progress.add_task("Bright stars (layers)", total=len(extents))
                batch = []
                for row in brightStars.scatter(
                    skeleton.shape, extents, skeleton.edge_pc, skeleton.expected_system_count_at_density_1,
                    min_luminosity_sol, seed, skip_addresses=filled,
                    on_layer=lambda done, _total: progress.update(task, completed=done),
                ):
                    counts[row[6]] += 1
                    batch.append(row)
                    if len(batch) >= 10000:
                        _db.insert_bright_stars(conn, batch)
                        conn.commit()
                        batch = []
                if batch:
                    _db.insert_bright_stars(conn, batch)
                _db.record_bright_star_scatter(conn, min_luminosity_sol, seed)
                conn.commit()
            finally:
                log.reset_console()
    finally:
        conn.close()
    elapsed = time.perf_counter() - t0
    total = sum(counts.values())
    log.normal(
        f"Placed {total:,} bright stars (at least {min_luminosity_sol:g} L_sun) in {elapsed:.1f}s: "
        + ", ".join(f"{count:,} {population}" for population, count in counts.items()) + "."
    )
    log.debug(f"bright-star scatter: {total} stars, seed {seed}, threshold {min_luminosity_sol:g} L_sun")
    return {"counts": counts, "total": total, "elapsed_s": elapsed}


def run_plan(args):
    """
    Builds and persists the galaxy's density skeleton, then scatters its
    bright stars (unless `--no-bright-stars`; `--bright-stars-only` skips
    the rebuild).

    Args:
        args (argparse.Namespace): Validated arguments (`command ==
            "plan"`).
    """
    if getattr(args, "bright_stars_only", False):
        scatter_bright_stars(args)
        return
    summary = build_skeleton(args)
    if summary["layer_count"]:
        layers = f"{summary['layer_count']} layers ({summary['top_layer_index']} to {-summary['top_layer_index']})"
    else:
        layers = "no layers (nothing clears the one-star-per-sector threshold)"
    log.normal(
        f"Skeleton built in {summary['elapsed_s']:.2f}s: {layers}, outer edge = ring "
        f"{summary['outer_ring_index']}, ~{summary['total_candidate_sectors']:,} candidate sectors."
    )
    if not summary["edge_confirmed"]:
        log.normal(
            f"WARNING: the galactic plane still qualified at --max-ring ({args.max_ring}) -- the "
            f"galaxy's true edge was not reached. Re-run with a larger --max-ring if these shape "
            f"parameters really do produce a galaxy this large."
        )
    if summary["layer_count"] and not getattr(args, "no_bright_stars", False):
        scatter_bright_stars(args)


# ===========================================================================
# 5. Exotic phenomena
# ===========================================================================

TYPE_LABELS = {
    "black-hole": "black hole",
    "neutron-star": "neutron star",
    "nebula": "nebula",
    "supernova-remnant": "supernova remnant",
    "rogue-planet": "rogue planet",
    "comet": "interstellar comet",
    "asteroid-field": "asteroid field",
    "quasar": "quasar",
}
"""dict: `--type` value -> human-readable label, used in the `phenomenon`
subcommand's own status output (not part of the generated phenomenon's
own text)."""


def add_phenomenon_arguments(parser):
    """
    Adds every option the `phenomenon` subcommand accepts (besides
    `--version`) to `parser` -- `--type`, `--anchor-system`,
    `--num-orbits`, `--sector-id`, the MySQL connection args, `--markdown`,
    `--name`, and the logging options (see `add_logging_arguments`).

    Args:
        parser (argparse.ArgumentParser): The parser to add options to.
    """
    parser.add_argument('--type', type=str, choices=list(program_constants.PHENOMENON_TYPE_CHOICES),
                         help="The kind of phenomenon to generate. Omit to pick uniformly at random "
                              "(never a quasar, which only exists at a galaxy's center).")

    parser.add_argument('--anchor-system', action='store_true',
                         help="Only valid with --type black-hole or --type neutron-star: build a full star "
                              "system (with orbiting planets, if any) anchored by the compact remnant, "
                              "instead of describing it standalone.")

    parser.add_argument('--num-orbits', type=int,
                         help="With --anchor-system, force an exact number of orbital slots around the "
                              "compact remnant.")

    parser.add_argument('--sector-id', type=int,
                         help="Valid with --type nebula, asteroid-field, black-hole, or neutron-star (not "
                              "combined with --anchor-system): place the generated phenomenon in the galaxy "
                              "near this already-generated, already galaxy-placed sector (see "
                              "stellarObjects._db.compute_phenomenon_placement). Also links a "
                              "supernova-remnant/rogue-planet/comet to that sector, without a computed "
                              "galaxy position (those types have no placement columns of their own). A "
                              "quasar needs a ring-0, layer-0 (galactic core) sector and is placed at the galactic "
                              "center; only one per galaxy. Omit to generate it unplaced/unlinked, as before.")

    # Database persistence
    _db.add_mysql_connection_args(parser)

    # Output in Markdown format
    parser.add_argument('--markdown', '-m', action='store_true', help="Output in Markdown format.")

    # Name
    parser.add_argument('--name', type=str, help="Force the name of the generated phenomenon.")

    # Logging (--debug, --quiet/--silent)
    add_logging_arguments(parser)


def validate_phenomenon_args(args, parser):
    """
    Validates `add_phenomenon_arguments`'s own options, calling
    `parser.error` (which exits) on the first problem found.

    Args:
        args (argparse.Namespace): Parsed arguments.
        parser (argparse.ArgumentParser): The parser to raise errors
                                          through (so the caller's own
                                          `--help`/usage text is shown).
    """
    if args.num_orbits is not None and not args.anchor_system:
        parser.error("--num-orbits requires --anchor-system.")
    if args.num_orbits is not None and args.num_orbits < 0:
        parser.error("--num-orbits must be zero or a positive integer.")
    if args.num_orbits is not None and args.num_orbits > generationLimits.MAX_NUM_ORBITS:
        parser.error(f"--num-orbits must be at most {generationLimits.MAX_NUM_ORBITS}.")
    if args.anchor_system and args.type not in (None, "black-hole", "neutron-star"):
        parser.error("--anchor-system is only valid with --type black-hole or --type neutron-star.")
    if args.sector_id is not None and args.anchor_system:
        parser.error("--sector-id cannot be combined with --anchor-system -- an anchored compact remnant "
                     "belongs to its own StarSystem, not directly to a sector.")


def generate_phenomenon(phenomenon_type, system_config, anchor_system, name=None):
    """
    Builds one generated phenomenon object for `phenomenon_type`.

    Args:
        phenomenon_type (str): One of `program_constants.PHENOMENON_TYPE_CHOICES`.
        system_config (SystemConfig): The shared config to generate with.
        anchor_system (bool): Only consulted for `"black-hole"`/
            `"neutron-star"` -- if True, wraps the generated compact
            remnant in a full `StarSystem` (see
            `StarSystem.__init__`'s `compact_remnant` parameter) instead
            of returning it standalone.
        name (str, optional): An explicit name for the phenomenon (or, for
            an anchored system, for the compact remnant itself).

    Returns:
        The generated phenomenon: a `BlackHole`, `NeutronStar`, `StarSystem`
        (only when `anchor_system` is True), `Nebula`, `SupernovaRemnant`,
        `RoguePlanet`, `InterstellarComet`, `AsteroidField`, or `Quasar`.

    Raises:
        ValueError: If `phenomenon_type` isn't one of the recognized choices.
    """
    if phenomenon_type == "black-hole":
        remnant = BlackHole(system_config, name=name)
        return StarSystem(system_config=system_config, compact_remnant=remnant) if anchor_system else remnant
    if phenomenon_type == "neutron-star":
        remnant = NeutronStar(system_config, name=name)
        return StarSystem(system_config=system_config, compact_remnant=remnant) if anchor_system else remnant
    if phenomenon_type == "nebula":
        return Nebula(system_config, name=name)
    if phenomenon_type == "supernova-remnant":
        return SupernovaRemnant(system_config, name=name)
    if phenomenon_type == "rogue-planet":
        return RoguePlanet(system_config, name=name)
    if phenomenon_type == "comet":
        return InterstellarComet(system_config, name=name)
    if phenomenon_type == "asteroid-field":
        return AsteroidField(system_config, name=name)
    if phenomenon_type == "quasar":
        return Quasar(system_config, name=name)

    raise ValueError(f"Unknown phenomenon type: {phenomenon_type!r}")


def run_phenomenon(args):
    """
    Generates one exotic phenomenon and saves it to the database.

    Args:
        args (argparse.Namespace): Validated arguments (`command ==
            "phenomenon"`).
    """
    phenomenon_type = args.type or random.choice(program_constants.RANDOM_PHENOMENON_TYPE_CHOICES)

    if args.sector_id is not None:
        # Checked before generating anything, as a clean one-line error.
        conn = _db.get_connection(_db.mysql_config_from_args(args))
        try:
            _db.get_sector_galaxy_position(conn, args.sector_id)
        except ValueError as exc:
            log.error(f"Error: --sector-id {args.sector_id}: {exc}")
            raise SystemExit(1) from exc
        finally:
            conn.close()

    system_config = SystemConfig()
    system_config.MARKDOWN = args.markdown
    if args.num_orbits is not None:
        system_config.NUM_ORBITS = args.num_orbits

    phenomenon = generate_phenomenon(phenomenon_type, system_config, args.anchor_system, name=args.name)

    mysql_config = _db.mysql_config_from_args(args)
    phenomenon_id = _db.save_phenomenon(phenomenon, system_config, phenomenon_type, config=mysql_config,
                                         sector_id=args.sector_id)
    RUN_COUNTS["phenomena"] += 1
    log.normal(f"Saved {TYPE_LABELS[phenomenon_type]} to the database (id={phenomenon_id}, "
               f"{mysql_config.database}@{mysql_config.host}:{mysql_config.port}).")


# ===========================================================================
# 6. Population and politics (POP.1 to POP.4)
# ===========================================================================

def add_population_arguments(parser):
    """
    Adds the `population` subcommand's options: `--rescan`,
    `--territories-only`, and the MySQL connection and logging args.

    Args:
        parser (argparse.ArgumentParser): The parser to add options to.
    """
    parser.add_argument('--rescan', action='store_true',
                        help="Forget every species, polity and territory and scan every planet again "
                             "(new names, ages and borders).")
    parser.add_argument('--territories-only', action='store_true',
                        help="Only recompute which polity owns which system.")
    _db.add_mysql_connection_args(parser)
    add_logging_arguments(parser)


def validate_population_args(args, parser):
    """`--rescan` and `--territories-only` contradict each other."""
    if args.rescan and args.territories_only:
        parser.error("--rescan and --territories-only cannot be combined.")


def _population_summary(counts):
    return (f"{counts['new_species']:,} new species; {counts['species']:,} species in all, "
            f"{counts['spacefaring']:,} spacefaring; {counts['polities']:,} polities holding "
            f"{counts['owned_systems']:,} systems.")


def run_population(args):
    """
    Runs the population pass (`population.run_pass`): names the dominant
    species of every new life world, dates civilizations, founds polities
    and recomputes territories. See docs/design/population-and-politics.md.

    Args:
        args (argparse.Namespace): Validated arguments (`command ==
            "population"`).
    """
    conn = _db.get_connection(_db.mysql_config_from_args(args))
    try:
        counts = population.run_pass(conn, rescan=args.rescan, territories_only=args.territories_only)
    finally:
        conn.close()
    log.normal(f"Population: {_population_summary(counts)}")


def run_population_after(args):
    """The population pass after a `sector` or `galaxy` run, only with
    `--population` (off by default, Boss 2026-10-01)."""
    if not getattr(args, "population", False):
        return
    conn = _db.get_connection(_db.mysql_config_from_args(args))
    try:
        counts = population.run_pass(conn)
    finally:
        conn.close()
    log.normal(f"Population: {_population_summary(counts)}")


# ===========================================================================
# 7. Unified CLI
# ===========================================================================

def build_parser():
    """
    Builds the top-level parser and every subcommand's own subparser.

    Returns:
        tuple: `(parser, subparsers)` -- `parser` is the top-level
              `argparse.ArgumentParser`; `subparsers` is a `dict` mapping
              each command name (`"system"`, `"sector"`, `"galaxy"`,
              `"plan"`, `"phenomenon"`) to its own subparser, so
              `process_args` can route `parser.error`-raised messages
              (and `--help`/usage text) through the right one.
    """
    parser = argparse.ArgumentParser(
        prog='generate.py',
        description="Unified Generation CLI",
        epilog="Generates a star system, sector, galaxy, galaxy density skeleton, or exotic phenomenon, "
               "and saves it to the database -- see this script's module docstring for each subcommand's "
               "own section. Run 'generate.py <command> --help' for that command's own full option list.")
    parser.add_argument('--version', action=VersionAction, banner=version_banner('generate.py'))

    subparsers = parser.add_subparsers(dest='command', required=True)

    system_parser = subparsers.add_parser(
        'system', prefix_chars='-+',
        description="System Generation Options",
        help="Generate a single star system.")
    add_system_arguments(system_parser)

    sector_parser = subparsers.add_parser(
        'sector', prefix_chars='-+',
        description="Sector Generation Options",
        help="Generate one or more independent sectors.")
    add_shared_generation_options(sector_parser)
    add_sector_arguments(sector_parser)

    galaxy_parser = subparsers.add_parser(
        'galaxy', prefix_chars='-+',
        description="Galaxy Generation Options",
        help="Generate many sectors as one galaxy.")
    add_shared_generation_options(galaxy_parser)
    add_galaxy_arguments(galaxy_parser)

    plan_parser = subparsers.add_parser(
        'plan',
        description="Galaxy Density Skeleton Builder",
        help="Build/replace the galaxy's density skeleton.")
    add_plan_arguments(plan_parser)

    phenomenon_parser = subparsers.add_parser(
        'phenomenon',
        description="Exotic Stellar Phenomenon Generation Options",
        help="Generate a single exotic stellar phenomenon.")
    add_phenomenon_arguments(phenomenon_parser)

    population_parser = subparsers.add_parser(
        'population',
        description="Population and Politics Pass",
        help="Name species, date civilizations and draw territories from what is stored.")
    add_population_arguments(population_parser)

    return parser, {
        'system': system_parser,
        'sector': sector_parser,
        'galaxy': galaxy_parser,
        'plan': plan_parser,
        'phenomenon': phenomenon_parser,
        'population': population_parser,
    }


def process_args():
    """
    Parses command-line arguments for `generate.py`, then validates
    whichever subcommand was chosen through its own `validate_*_args`
    function(s).

    Returns:
        argparse.Namespace: Parsed (and validated) arguments, with
            `args.command` set to the chosen subcommand name.
    """
    parser, command_parsers = build_parser()
    args = parser.parse_args()
    command_parser = command_parsers[args.command]

    validate_logging_args(args, command_parser)

    port = getattr(args, "mysql_port", None)
    if port is not None and not 1 <= port <= 65535:
        command_parser.error("--mysql-port must be between 1 and 65535.")

    if args.command == 'system':
        validate_system_args(args, command_parser)
    elif args.command == 'sector':
        validate_shared_generation_args(args, command_parser)
        validate_sector_args(args, command_parser)
    elif args.command == 'galaxy':
        validate_shared_generation_args(args, command_parser)
        validate_galaxy_args(args, command_parser)
    elif args.command == 'plan':
        validate_plan_args(args, command_parser)
    elif args.command == 'phenomenon':
        validate_phenomenon_args(args, command_parser)
    elif args.command == 'population':
        validate_population_args(args, command_parser)

    return args


_COMMAND_HANDLERS = {
    'system': run_system,
    'sector': run_sector,
    'galaxy': run_galaxy,
    'plan': run_plan,
    'phenomenon': run_phenomenon,
    'population': run_population,
}


def main():
    """
    The main entry point for the unified generation CLI. Parses and
    validates command-line arguments, configures logging severity from
    `--debug`/`--quiet`/`--silent`, seeds the random number generator
    cryptographically, then dispatches to the chosen subcommand's own
    `run_*` function.
    """
    args = process_args()

    if args.quiet:
        level = log.SILENT
    elif args.debug is not None:
        level = log.DEBUG
    else:
        level = log.NORMAL
    try:
        log.configure(level, debug_file=(args.debug or None))
    except OSError as exc:
        _fatal(f"cannot open --debug file {args.debug!r}: {exc.strerror or exc}", logger_ready=False)
    log.debug("Command: %s, options: %s", args.command,
              {key: ("<withheld>" if "password" in key else value) for key, value in sorted(vars(args).items())})

    seed = secrets.randbits(128)
    random.seed(seed)
    log.debug(f"Seeded the random number generator with {seed} (cryptographically random; no --seed "
              f"option exists to reproduce this run).")

    # One start and one finish line per run in the activity log (SEC.28);
    # a `system --output` run writes no database, so it isn't logged.
    logged = not getattr(args, "output", None)
    try:
        database = _db.mysql_config_from_args(args).database
    except AttributeError:  # a subcommand without the --mysql-* options
        database = _db.DEFAULT_MYSQL_CONFIG.database
    started = time.monotonic()
    if logged:
        activitylog.event("GEN", "generate.start", user=_run_user(), command=args.command, db=database)
    status = "failed"
    try:
        _COMMAND_HANDLERS[args.command](args)
        status = "ok"
    except pymysql.err.MySQLError as exc:
        _fatal(f"database error: {exc}")
    except OSError as exc:
        # e.g. an unwritable --output path.
        where = f" ({exc.filename})" if getattr(exc, "filename", None) else ""
        _fatal(f"{exc.strerror or exc}{where}")
    except KeyboardInterrupt:
        status = "interrupted"
        raise
    finally:
        if logged:
            activitylog.event("GEN", "generate.finish", user=_run_user(), command=args.command, db=database,
                              status=status, seconds=round(time.monotonic() - started, 1),
                              **{key: RUN_COUNTS[key] for key in ("sectors", "systems", "phenomena")})


def _run_user():
    """Who ran this: the login name (`getpass.getuser`), or `None`."""
    try:
        return getpass.getuser()
    except Exception:  # noqa: BLE001 -- no user name in this environment
        return None


def _fatal(message, logger_ready=True):
    """A one-line error and exit status 1, instead of a traceback. Goes
    through `log.error` (shown even under --quiet/--silent) once the
    logger is configured, else straight to stderr."""
    if logger_ready:
        log.error(f"Error: {message}")
    else:
        print(f"generate.py: error: {message}", file=sys.stderr)
    raise SystemExit(1)


if __name__ == "__main__":
    main()
