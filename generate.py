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
   `galaxy_ring_band` summary `ensure_sector_generated` consults to
   decide, cheaply and exactly, whether a given address is worth
   generating at all, without ever enumerating the galaxy's ~10 billion
   candidate sector slots.
5. Exotic phenomena (`phenomenon`) -- black holes, neutron stars,
   nebulae, supernova remnants, rogue planets, interstellar comets, and
   standalone asteroid fields, generated on demand; section 2 also
   reuses this section's `generate_phenomenon` to seed every sector with
   its own sparse, science-based population of the same seven types.
6. The unified CLI itself (argument parsing/validation, dispatch,
   `main`).
"""

import argparse
import copy
import logging
import math
import os
import random
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

from stellarObjects import _db, log, program_constants
from stellarObjects._version import VersionAction, version_banner
from stellarObjects.asteroidFieldData import AsteroidField
from stellarObjects.compactRemnant import BlackHole, NeutronStar
from stellarObjects.config import SystemConfig
from stellarObjects.galaxyDensity import build_galaxy_shape, predicted_star_count, relative_density
from stellarObjects.galaxyGeometry import (
    SectorCell, enumerate_sectors_within_radius, galactic_radius_pc,
    provisional_sector_designation, ring_bounds_pc, ring_sector_count, sector_position_pc,
)
from stellarObjects.galaxySkeleton import (
    DEFAULT_EMPTY_STREAK_TO_STOP, DEFAULT_MAX_RING, build_ring_bands, expected_system_count_at_density_1,
)
from stellarObjects.nebulaData import Nebula
from stellarObjects.roguePlanetData import InterstellarComet, RoguePlanet
from stellarObjects.spaceSector import SpaceSector, _sample_poisson_count
from stellarObjects.supernovaRemnantData import SupernovaRemnant
from stellarObjects.systemData import StarSystem
from stellarObjects.utils import generate_sector_name, ly_to_pc, pc_to_ly

# Suppress transformers warnings
logging.getLogger("transformers").setLevel(logging.ERROR)


def _generation_progress():
    """
    Builds the shared `rich.progress.Progress` used by `run_galaxy`'s three
    modes -- one "Sectors" task per run tracking how many sectors have been
    generated so far. Deliberately not used by `run_sector`'s own
    `--num-sectors` loop: that command's own run is normally short enough
    (and its sector count small enough) that a bar added more visual noise
    than it was worth, whereas a `galaxy` run (a whole ring, a
    neighborhood, or a random start's default 100 ly one) can mean
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

    Returns:
        Progress: Not yet started.
    """
    return Progress(
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
    parser.add_argument('--flavor-chance-system', type=float,
                        help="Override the default FLAVOR_CHANCE_SYSTEM constant.")

    # Override Flavor Chance Planet
    parser.add_argument('--flavor-chance-planet', type=float,
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

    if args.intelligent_life is not None and args.habitable_world is False:
        parser.error("+intelligent_life/-intelligent_life cannot be combined with -habitable_world.")

    if args.num_orbits is not None and args.num_orbits < 0:
        parser.error("--num-orbits must be zero or a positive integer.")

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
    Generates one star system and saves it to the database.

    Args:
        args (argparse.Namespace): Validated arguments (`command ==
            "system"`).
    """
    system_config = build_system_config(args)
    system = StarSystem(system_config=system_config)

    mysql_config = _db.mysql_config_from_args(args)
    star_system_id = _db.save_system(system, system_config, config=mysql_config)
    log.normal(f"Saved system '{system.star.name}' to the database (star_system_id={star_system_id}, "
               f"{mysql_config.database}@{mysql_config.host}:{mysql_config.port}).")


# ===========================================================================
# 2. Sector generation
# ===========================================================================

_PHENOMENON_FACTORIES = {
    "black-hole": lambda config, galactic_center_dist_ly: BlackHole(
        config, galactic_center_dist_ly=galactic_center_dist_ly),
    "neutron-star": lambda config, galactic_center_dist_ly: NeutronStar(
        config, galactic_center_dist_ly=galactic_center_dist_ly),
    "nebula": lambda config, galactic_center_dist_ly: Nebula(config),
    "supernova-remnant": lambda config, galactic_center_dist_ly: SupernovaRemnant(config),
    "rogue-planet": lambda config, galactic_center_dist_ly: RoguePlanet(config),
    "comet": lambda config, galactic_center_dist_ly: InterstellarComet(config),
    "asteroid-field": lambda config, galactic_center_dist_ly: AsteroidField(config),
}
"""dict: `program_constants.PHENOMENON_TYPE_CHOICES` entry -> a
`(config, galactic_center_dist_ly)` factory building one fresh instance of
that phenomenon type -- `generate_sector_phenomena`'s per-type dispatch.
Only `"black-hole"`/`"neutron-star"` actually consult
`galactic_center_dist_ly` (threaded into their own Hill-sphere/galactic-
orbit calculations, exactly like every star in the sector already gets via
`generate_sector`'s own `galactic_center_dist_ly` parameter); every other
phenomenon type has no galaxy-frame-distance-dependent physics of its own,
so it's accepted and ignored, keeping every factory the same shape."""


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
    parser.add_argument('--density', type=float, default=None,
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
    parser.add_argument('--flavor-chance-system', type=float,
                        help="Override the default FLAVOR_CHANCE_SYSTEM constant.")
    parser.add_argument('--flavor-chance-planet', type=float,
                        help="Override the default FLAVOR_CHANCE_PLANET constant.")
    parser.add_argument('--max-planet-flavor', action='store_true',
                        help="Sets the maximum flavor text total for planets to 99.")

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
    configs = [build_system_config(args) for _ in range(args.num_systems)]

    if args.min_habitable > len(configs):
        log.error(
            f"Error: --min-habitable ({args.min_habitable}) exceeds this sector's generated system "
            f"count ({len(configs)}); with --density, the count is randomly sampled per sector and can "
            f"land below --min-habitable. Try a smaller --min-habitable, a higher --density, or "
            f"--num-systems for an exact count instead."
        )
        raise SystemExit(1)

    if args.min_habitable > 0:
        already_habitable = [i for i, cfg in enumerate(configs) if cfg.HABITABLE_WORLD is True]
        still_needed = args.min_habitable - len(already_habitable)
        if still_needed > 0:
            candidates = [i for i in range(len(configs)) if i not in already_habitable]
            for i in random.sample(candidates, k=still_needed):
                configs[i].HABITABLE_WORLD = True
                # Mirrors build_system_config's own habitable-world +
                # asteroid-belt normalization, reapplied here since it ran before
                # this override existed.
                if configs[i].ASTEROID_BELT is True:
                    if configs[i].LARGE_STAR is False:
                        log.error(
                            "Error: --min-habitable requires forcing a habitable world onto a system that also "
                            "has +asteroid_belt forced sector-wide; that combination needs a large star, but "
                            "-large_star was also forced sector-wide."
                        )
                        raise SystemExit(1)
                    configs[i].LARGE_STAR = True

    return configs


def generate_sector_phenomena(sector, args, galactic_center_dist_ly=None):
    """
    Populates an already-built `sector` with a realistically sparse
    population of exotic stellar phenomena (see section 5's own seven
    generated types), sampled independently per type from
    `program_constants.PHENOMENON_RATE_PER_STAR_SYSTEM` -- each type's own
    expected count per star system, scaled by however many star systems
    this sector actually ended up with (see that constant's own docstring
    for where each rate comes from). The Poisson draw
    (`spaceSector._sample_poisson_count`) is the same mechanism
    `--density`'s own system count already uses, so a denser sector gets
    proportionally more phenomena too, and most sectors -- realistically --
    get none at all.

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
    system_count = len(sector.entries)
    new_entries = []

    for phenomenon_type, rate_per_system in program_constants.PHENOMENON_RATE_PER_STAR_SYSTEM.items():
        count = _sample_poisson_count(rate_per_system * system_count)
        for _ in range(count):
            phenomenon_config = SystemConfig()
            phenomenon_config.MARKDOWN = args.markdown
            phenomenon = _PHENOMENON_FACTORIES[phenomenon_type](phenomenon_config, galactic_center_dist_ly)
            try:
                new_entries.append(sector.add_phenomenon(phenomenon, phenomenon_type))
            except ValueError:
                # Only a massive type (black-hole/neutron-star) can raise
                # here (see SpaceSector._random_position) -- no room left
                # to place it without violating another massive object's
                # Hill sphere. Skip just this one draw rather than this
                # phenomenon type, or the whole sector.
                continue

    return new_entries


def generate_sector(args, galactic_center_dist_ly=None, cell=None):
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

    # `--density` resolves to a concrete system count per sector (this
    # sector's own volume, sampled fresh each call) rather than once at
    # parse time -- this is what lets each sector under `--num-sectors`/
    # the `galaxy` subcommand vary independently instead of sharing one
    # fixed count. A copy avoids mutating the caller's shared `args`
    # namespace, since this function runs once per sector.
    density_driven = args.density is not None
    if density_driven:
        args = copy.copy(args)
        args.num_systems = _sample_poisson_count(sector.expected_system_count() * args.density)

    with log.timed_phase("build_sector_configs"):
        configs = build_sector_configs(args)

    for i, cfg in enumerate(configs):
        with log.timed_phase(f"generate system {i + 1}/{len(configs)}"):
            system = StarSystem(system_config=cfg, galactic_center_dist_ly=galactic_center_dist_ly)

        try:
            with log.timed_phase(f"place system {i + 1}/{len(configs)}"):
                sector.add_system(system, system_config=cfg)
        except ValueError:
            # No room left for another system's Hill sphere in this
            # sector's cube -- see this function's own "Capacity"
            # docstring note. Every following config would almost
            # certainly fail the same way (this sector only gets fuller
            # from here), so stop generating and placing altogether
            # rather than pay for `len(configs) - i - 1` more full
            # StarSystem generations just to discard them too.
            log.normal(
                f"Sector '{sector_name}': ran out of room after placing {len(sector.entries)} of "
                f"{len(configs)} requested systems -- the {sector.edge_ly:.1f} ly cube has no space left "
                f"that clears every already-placed system's Hill sphere. Returning the sector as-is "
                f"rather than the full requested count."
            )
            break

    with log.timed_phase("generate_sector_phenomena"):
        generate_sector_phenomena(sector, args, galactic_center_dist_ly=galactic_center_dist_ly)

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
    for _i in range(args.num_sectors):
        _sector_name, sector = generate_sector(args)
        systems = [entry.star_system for entry in sector.entries]

        mysql_config = _db.mysql_config_from_args(args)
        sector_id = _db.save_sector(sector, config=mysql_config)

        phenomena_note = f", {len(sector.phenomena)} phenomena" if sector.phenomena else ""
        log.normal(
            f"Saved sector '{sector.name}' to the database (sector_id={sector_id}, "
            f"{len(systems)} systems{phenomena_note}, "
            f"{mysql_config.database}@{mysql_config.host}:{mysql_config.port})."
        )
        log.normal(sector_generation_summary_lines(sector, args))

    if args.num_sectors > 1:
        log.normal(f"Generated {args.num_sectors} sectors.")


# ===========================================================================
# 3. Galaxy generation
# ===========================================================================

LARGE_RING_WARNING_THRESHOLD = 2000
"""int: `--ring I` requires `--limit` or `--yes` when ring `I` holds more
slots than this (see `ring_sector_count`) -- about ring 318, ~3,700 ly
out. Anything larger takes a real, unbounded amount of time and disk, so
it needs an explicit choice."""

RANDOM_START_MAX_HEIGHT_PC = 1000.0
"""float: How far above or below the plane random-start mode draws its
seed sector, parsecs -- a generous thick-disk half-height. Draws outside
the stored skeleton's band are simply retried."""


def add_galaxy_arguments(parser):
    """
    Adds the `galaxy` subcommand's own mode/placement options -- on top
    of whatever `add_shared_generation_options` already added -- to
    `parser`: the `--ring`/`--center-sector` mode-selection group,
    `--layer`, `--slot`, `--limit`, `--yes`, `--radius-pc`, `--max-ring`,
    `--min-start-density`, and the MySQL connection args.

    Args:
        parser (argparse.ArgumentParser): The parser to add options to.
    """
    mode_group = parser.add_mutually_exclusive_group(required=False)
    mode_group.add_argument('--ring', type=int, metavar='I',
                            help="Batch mode: generate every not-yet-generated sector in ring I (0-indexed "
                                 "cylindrical radius band, one sector edge wide) at --layer J.")
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
    parser.add_argument('--limit', type=int,
                        help="With --ring: generate only the first N not-yet-generated slots.")
    parser.add_argument('--yes', action='store_true',
                        help="With --ring: skip the confirmation normally required before generating a "
                             "ring whose slot count exceeds LARGE_RING_WARNING_THRESHOLD.")
    parser.add_argument('--radius-pc', type=float,
                        help="With --center-sector: the neighborhood search radius, in parsecs. With "
                             "neither --ring nor --center-sector (random-start mode): overrides the "
                             "default 100 ly neighborhood radius around the randomly chosen starting "
                             "sector.")
    parser.add_argument('--max-ring', type=int,
                        help="With neither --ring nor --center-sector (random-start mode): the highest "
                             f"ring the randomly chosen starting sector may land in. Defaults to the ring "
                             f"at a real Milky-Way-scale galaxy radius "
                             f"({program_constants.GALAXY_RADIUS_PC:,.0f} pc).")
    parser.add_argument('--min-start-density', type=float,
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
    random_start = args.ring is None and args.center_sector is None

    if args.ring is not None and args.ring < 0:
        parser.error("--ring must be >= 0.")
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
    if args.slot is not None and args.radius_pc is not None:
        parser.error("--radius-pc doesn't apply to --ring --slot (single-address mode).")
    if args.slot is not None and (args.density is not None or args.num_systems is not None):
        parser.error("--density/--num-systems can't be combined with --ring --slot (single-address "
                     "mode) -- ensure_sector_generated always uses this address's own real predicted "
                     "density, the same as when the galaxy map's own live view found it.")

    if args.center_sector is not None and args.radius_pc is None:
        parser.error("--center-sector requires --radius-pc.")
    if args.radius_pc is not None and args.ring is not None:
        parser.error("--radius-pc only applies to --center-sector or random-start mode (neither "
                     "--ring nor --center-sector), not --ring.")
    if args.radius_pc is not None and args.radius_pc <= 0:
        parser.error("--radius-pc must be a positive number.")

    if args.limit is not None and args.ring is None:
        parser.error("--limit only applies to --ring.")
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be a positive integer.")
    if args.yes and args.ring is None:
        parser.error("--yes only applies to --ring.")

    if args.max_ring is not None and not random_start:
        parser.error("--max-ring only applies to random-start mode (neither --ring nor --center-sector).")
    if args.max_ring is not None and args.max_ring < 0:
        parser.error("--max-ring must be >= 0.")

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
    Resolves each sector's own `--density` from the galaxy skeleton's real
    position-based `relative_density`, for `galaxy` mode's batch/local-
    neighborhood/random-start generation -- the same mechanism
    `ensure_sector_generated` uses for a single lazily-generated sector,
    applied across a whole run instead of every sector sharing one flat
    CLI value.

    A no-op (`resolve` returns `args` unchanged) whenever the operator
    explicitly passed `--density`/`--num-systems` -- an explicit flag is a
    deliberate, uniform override for the whole run.

    Otherwise `resolve` also applies `ensure_sector_generated`'s
    qualification gate (the stored ring band, then the exact
    `predicted_star_count`), returning `None` for an address that doesn't
    qualify, so batch runs never save all-but-certainly-empty sectors.

    The skeleton is fetched once, on first use, and each ring's band is
    cached -- neither changes mid-run.
    """

    def __init__(self, config):
        self._config = config
        self._skeleton = None
        self._bands_cache = {}

    def _get_skeleton(self):
        if self._skeleton is None:
            conn = _db.get_connection(self._config)
            try:
                self._skeleton = _db.get_galaxy_shape(conn)
            finally:
                conn.close()
            if self._skeleton is None:
                raise RuntimeError(
                    "Neither --density nor --num-systems was given, and the galaxy's skeleton has "
                    "never been built (no galaxy_shape row) -- run 'generate.py plan' first, or pass "
                    "--density/--num-systems explicitly to skip per-sector skeleton density."
                )
        return self._skeleton

    def _get_band(self, ring_index):
        if ring_index not in self._bands_cache:
            conn = _db.get_connection(self._config)
            try:
                self._bands_cache[ring_index] = _db.get_galaxy_ring_band(conn, ring_index)
            finally:
                conn.close()
        return self._bands_cache[ring_index]

    def resolve(self, args, address, position_pc):
        """
        Args:
            args (argparse.Namespace): The `galaxy` subcommand's own parsed
                (and validated) arguments.
            address (tuple): This sector's `(ring, layer, slot)`.
            position_pc (tuple): This sector's `(x, y, z)` center, parsecs.

        Returns:
            argparse.Namespace or None: `args` itself when a density/count
                was given explicitly. Otherwise `None` if this address
                doesn't qualify (outside its ring's stored band, or its own
                exact `predicted_star_count < 1.0`), else a fresh copy with
                `.density` set to this position's own `relative_density`
                and `.num_systems` cleared.
        """
        if args.density is not None or args.num_systems is not None:
            return args
        skeleton = self._get_skeleton()

        band = self._get_band(address[0])
        if band is None or not (band[0] <= address[1] <= band[1]):
            return None

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
    `program_constants.DEFAULT_SECTOR_EDGE_LY` converted (sector edge
    length is not an exposed CLI option for `galaxy`).
    """
    return ly_to_pc(program_constants.DEFAULT_SECTOR_EDGE_LY)


def generate_and_save_sector_at(args, address, position_pc, edge_pc):
    """
    Generates one sector via `generate_sector` inside its real grid cell
    and saves it -- the per-sector unit of work every `galaxy` mode
    repeats.

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

    _sector_name, sector = generate_sector(args, galactic_center_dist_ly=pc_to_ly(radius_pc), cell=cell)

    galaxy_position = {
        "center_x_pc": x, "center_y_pc": y, "center_z_pc": z,
        "galactic_radius_pc": radius_pc,
        "ring_index": ring_index, "layer_index": layer_index, "ring_slot_index": slot_index,
    }
    sector_id = _db.save_sector(sector, config=_db.mysql_config_from_args(args), galaxy_position=galaxy_position)
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

        band = _db.get_galaxy_ring_band(conn, ring_index)
    finally:
        conn.close()

    position_pc = sector_position_pc(ring_index, layer_index, ring_slot_index, skeleton.edge_pc)
    if band is None or not (band[0] <= layer_index <= band[1]):
        # Outside the ring's stored band -- galaxySkeleton.find_ring_band's
        # bound is exact, so this is a certain "no".
        return {"created": False, "qualifies": False, "sector_id": None, "sector_name": None}

    density = relative_density(position_pc, skeleton.shape)
    star_count = predicted_star_count(position_pc, skeleton.shape, skeleton.expected_system_count_at_density_1)
    if star_count < 1.0:
        # Inside the band (a safe superset) but this slot's own angle
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


def _log_saved(sector_id, sector_name, sector, sector_args, address, suffix=""):
    designation = provisional_sector_designation(*address)
    log.normal(
        f"Saved sector '{sector_name}' [{designation}] at {_format_address(address)}{suffix} "
        f"(sector_id={sector_id})."
    )
    log.normal(sector_generation_summary_lines(sector, sector_args))


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
        SystemExit: If the ring's slot count exceeds
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
    conn = _db.get_connection(mysql_config)
    try:
        occupied = {a for a in _db.get_occupied_addresses(conn, [ring_index]) if a[1] == layer_index}
    finally:
        conn.close()

    batch_density = _BatchDensity(mysql_config)

    to_generate = total_slots - len(occupied)
    if args.limit is not None and (args.density is not None or args.num_systems is not None):
        # Only a safe cap when every unoccupied slot is sure to generate
        # (an explicit density bypasses the qualification gate).
        to_generate = min(to_generate, args.limit)
    outer_task = progress.add_task(f"Sectors (ring {ring_index} layer {layer_index})", total=max(to_generate, 0))

    generated = 0
    skipped = 0
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
                      f"outside the ring's stored band)")
            skipped += 1
            progress.update(outer_task, advance=1)
            continue

        log.debug(f"{_format_address(address)}: generating (density={sector_args.density})")
        sector_id, sector_name, sector = generate_and_save_sector_at(sector_args, address, position_pc, edge_pc)
        generated += 1
        progress.update(outer_task, advance=1)
        _log_saved(sector_id, sector_name, sector, sector_args, address)

    skip_note = f", {skipped} skipped (below the star-count threshold)" if skipped else ""
    log.normal(
        f"Generated {generated} new sector(s) in ring {ring_index} layer {layer_index} "
        f"({total_slots} total slots, {len(occupied)} already existed{skip_note})."
    )


def _neighborhood_candidates(center, radius_pc, edge_pc, config):
    """Every address within `radius_pc` of `center`, plus the set of
    those already occupied."""
    candidates = list(enumerate_sectors_within_radius(center, radius_pc, edge_pc))
    conn = _db.get_connection(config)
    try:
        occupied = _db.get_occupied_addresses(conn, {c[0] for c in candidates})
    finally:
        conn.close()
    return candidates, occupied


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
    candidates, occupied = _neighborhood_candidates(center, args.radius_pc, edge_pc, mysql_config)
    batch_density = _BatchDensity(mysql_config)

    to_generate = sum(1 for c in candidates if c[:3] not in occupied)
    outer_task = progress.add_task("Sectors (local neighborhood)", total=to_generate)

    generated = 0
    skipped = 0
    already_existed = 0
    for ring_index, layer_index, slot_index, x, y, z, distance_pc in candidates:
        address = (ring_index, layer_index, slot_index)
        if address in occupied:
            already_existed += 1
            continue

        sector_args = batch_density.resolve(args, address, (x, y, z))
        if sector_args is None:
            log.debug(f"{_format_address(address)}: skipped (below the 1-star-per-sector threshold, or "
                      f"outside the ring's stored band)")
            skipped += 1
            progress.update(outer_task, advance=1)
            continue

        log.debug(f"{_format_address(address)}: generating (density={sector_args.density})")
        sector_id, sector_name, sector = generate_and_save_sector_at(sector_args, address, (x, y, z), edge_pc)
        generated += 1
        progress.update(outer_task, advance=1)
        _log_saved(sector_id, sector_name, sector, sector_args, address,
                   suffix=f", {distance_pc:.2f} pc from sector_id={args.center_sector}")

    skip_note = f", {skipped} skipped (below the star-count threshold)" if skipped else ""
    log.normal(
        f"Generated {generated} new sector(s) within {args.radius_pc} pc of sector_id={args.center_sector} "
        f"({len(candidates)} candidate slot(s) found, {already_existed} already existed{skip_note})."
    )


def generate_sector_neighborhood(center_sector_id, radius_ly=None, config=None):
    """
    Non-CLI counterpart to `run_local_neighborhood` -- for the admin web
    UI's "generate more sectors around this one" action
    (`html/api/routes.py`'s `generate_sector_neighborhood_route`). Same
    work, a plain result dict instead of prints, and a catchable
    `ValueError` instead of `SystemExit` for an invalid/unplaced sector.

    The default 100 ly radius is large relative to one 11.5 ly sector: its
    sphere holds roughly 2,000-3,000 candidate addresses, so this can take
    minutes to hours depending on the server and how many already exist.

    Args:
        center_sector_id (int): The already galaxy-placed sector to
            generate a neighborhood around.
        radius_ly (float, optional): Defaults to
            `program_constants.RANDOM_START_NEIGHBORHOOD_RADIUS_LY` (100 ly).
        config (MySQLConfig, optional): Connection parameters.

    Returns:
        dict: `generated`, `already_existed`, `skipped`, `candidates` (all int).

    Raises:
        ValueError: If `center_sector_id` doesn't exist, or has never been
                   placed in a galaxy.
        RuntimeError: If the galaxy's skeleton has never been built.
    """
    edge_pc = _edge_pc()
    radius_pc = (
        ly_to_pc(radius_ly) if radius_ly is not None
        else ly_to_pc(program_constants.RANDOM_START_NEIGHBORHOOD_RADIUS_LY)
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
    candidates, occupied = _neighborhood_candidates(center, radius_pc, edge_pc, config)

    args = _default_generation_args(config=config)
    # Density-driven from the skeleton (_BatchDensity), not the flat
    # num_systems=10 default.
    args.density = None
    args.num_systems = None
    batch_density = _BatchDensity(config)

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
    }


def _pick_random_address(max_ring_index, edge_pc):
    """
    A random grid address uniformly by volume within a disk of
    `max_ring_index + 1` rings and `+/- RANDOM_START_MAX_HEIGHT_PC`:
    `R = R_max * sqrt(u)` is uniform by area over the disk, then the
    height and angle are uniform.

    Returns:
        tuple: `(ring, layer, slot)`.
    """
    r_max = ring_bounds_pc(max_ring_index, edge_pc)[1]
    r = r_max * math.sqrt(random.random())
    theta = random.uniform(0.0, 2 * math.pi)
    z = random.uniform(-RANDOM_START_MAX_HEIGHT_PC, RANDOM_START_MAX_HEIGHT_PC)
    ring_index = min(max_ring_index, int(r / edge_pc))
    layer_index = int(math.floor(z / edge_pc + 0.5))
    n = ring_sector_count(ring_index)
    slot_index = min(n - 1, int(theta * n / (2 * math.pi)))
    return ring_index, layer_index, slot_index


def run_random_start(args, edge_pc, progress):
    """
    Random-start mode (no `--ring`/`--center-sector` given): picks a
    random, not-yet-occupied, qualifying sector address (`_pick_random_address`,
    up to `--max-ring` or the ring at `program_constants.GALAXY_RADIUS_PC`),
    retried up to `program_constants.RANDOM_START_MAX_PLACEMENT_ATTEMPTS`
    times, generates it, then generates every not-yet-generated sector
    within `args.radius_pc` of it (default
    `program_constants.RANDOM_START_NEIGHBORHOOD_RADIUS_LY`, 100 ly) via
    `run_local_neighborhood`. `--min-start-density` tightens the retry:
    a qualifying address below it is retried too.

    Raises:
        SystemExit: If no suitable address was found within the attempt
                   budget.
    """
    max_ring_index = (
        args.max_ring if args.max_ring is not None
        else int(program_constants.GALAXY_RADIUS_PC / edge_pc)
    )
    radius_pc = (
        args.radius_pc if args.radius_pc is not None
        else ly_to_pc(program_constants.RANDOM_START_NEIGHBORHOOD_RADIUS_LY)
    )

    mysql_config = _db.mysql_config_from_args(args)
    batch_density = _BatchDensity(mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        sector_args = None
        for _ in range(program_constants.RANDOM_START_MAX_PLACEMENT_ATTEMPTS):
            address = _pick_random_address(max_ring_index, edge_pc)
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
            log.error(
                f"Could not find an unoccupied, qualifying sector address within {max_ring_index} "
                f"rings after {program_constants.RANDOM_START_MAX_PLACEMENT_ATTEMPTS} attempts{density_note} "
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
    an address copied out of the interactive 3D Galaxy Map.

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

    task = progress.add_task(f"Sector ({_format_address(address)})", total=1)
    mysql_config = _db.mysql_config_from_args(args)
    result = ensure_sector_generated(*address, config=mysql_config)
    progress.update(task, advance=1)

    if not result["qualifies"]:
        log.error(
            f"{_format_address(address)} doesn't qualify -- it would hold no real content at this "
            f"galaxy's own predicted density (below the 1-star-per-sector threshold, or outside the "
            f"ring's stored band)."
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


def run_galaxy(args):
    """
    Dispatches to single-address, ring-batch, local-neighborhood, or
    random-start mode, owning the one `rich.progress.Progress` display
    they share.

    Args:
        args (argparse.Namespace): Validated arguments (`command ==
            "galaxy"`).
    """
    edge_pc = _edge_pc()

    with _generation_progress() as progress:
        log.set_console(progress.console)
        try:
            if args.slot is not None:
                run_single_slot(args, edge_pc, progress)
            elif args.ring is not None:
                run_ring_batch(args, edge_pc, progress)
            elif args.center_sector is not None:
                run_local_neighborhood(args, edge_pc, progress)
            else:
                run_random_start(args, edge_pc, progress)
        finally:
            log.reset_console()


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
                        help=f"Sector edge length (ring width and layer height), light-years. Default: "
                             f"{program_constants.DEFAULT_SECTOR_EDGE_LY}.")
    parser.add_argument('--empty-streak-to-stop', type=int, default=DEFAULT_EMPTY_STREAK_TO_STOP,
                        help=f"Consecutive empty rings before concluding the galaxy's edge has been "
                             f"reached. Default: {DEFAULT_EMPTY_STREAK_TO_STOP}.")
    parser.add_argument('--max-ring', type=int, default=DEFAULT_MAX_RING,
                        help=f"Hard cap on rings scanned, regardless of --empty-streak-to-stop. "
                             f"Default: {DEFAULT_MAX_RING}.")
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
    if args.empty_streak_to_stop < 1:
        parser.error("--empty-streak-to-stop must be a positive integer.")
    if args.max_ring < 1:
        parser.error("--max-ring must be a positive integer.")


def build_skeleton(args):
    """
    Builds and persists the galaxy's skeleton: `galaxy_shape` (the shape
    parameters, calibration constant, edge length and outer ring -- one
    singleton row) and `galaxy_ring_band` (one row per ring that can hold
    content: its layer range, see `galaxySkeleton.find_ring_band`).

    No sector content or individual addresses are stored -- that stays
    lazy (`ensure_sector_generated`). Re-running replaces the whole
    skeleton; a full Milky-Way-scale build takes about half a second.

    Args:
        args (argparse.Namespace): Parsed arguments.

    Returns:
        dict: `rings_scanned`, `outer_ring_index`, `total_bands`,
              `total_candidate_sectors`, `elapsed_s`, `edge_confirmed`.
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

    log.normal(
        f"Building skeleton: disk_scale_length_pc={shape.disk_scale_length_pc} "
        f"disk_scale_height_pc={shape.disk_scale_height_pc} "
        f"bulge_scale_radius_pc={shape.bulge_scale_radius_pc} "
        f"bulge_amplitude={shape.bulge_amplitude} arm_count={shape.arm_count} "
        f"k_norm={shape.k_norm:.4f} threshold_rho={threshold_rho:.6f}"
    )

    t0 = time.perf_counter()
    bands, outer_ring_index, edge_confirmed = build_ring_bands(
        shape, edge_pc, threshold_rho,
        empty_streak_to_stop=args.empty_streak_to_stop, max_ring=args.max_ring,
    )
    elapsed = time.perf_counter() - t0

    mysql_config = _db.mysql_config_from_args(args)
    _db.replace_galaxy_ring_bands(bands, config=mysql_config)
    _db.save_galaxy_shape(
        shape, edge_pc=edge_pc, outer_ring_index=outer_ring_index,
        expected_system_count_at_density_1=e_value, config=mysql_config,
    )

    rings_scanned = (outer_ring_index + 1 + args.empty_streak_to_stop) if edge_confirmed else args.max_ring + 1
    total_candidate_sectors = sum((hi - lo + 1) * ring_sector_count(ring) for ring, lo, hi in bands)
    return {
        "rings_scanned": rings_scanned,
        "outer_ring_index": outer_ring_index,
        "total_bands": len(bands),
        "total_candidate_sectors": total_candidate_sectors,
        "elapsed_s": elapsed,
        "edge_confirmed": edge_confirmed,
    }


def run_plan(args):
    """
    Builds and persists the galaxy's density skeleton.

    Args:
        args (argparse.Namespace): Validated arguments (`command ==
            "plan"`).
    """
    summary = build_skeleton(args)
    log.normal(
        f"Skeleton built in {summary['elapsed_s']:.2f}s: scanned {summary['rings_scanned']} rings, "
        f"outer edge = ring {summary['outer_ring_index']}, {summary['total_bands']} band(s) stored, "
        f"~{summary['total_candidate_sectors']:,} candidate sectors."
    )
    if not summary["edge_confirmed"]:
        log.normal(
            f"WARNING: reached --max-ring ({args.max_ring}) without a run of "
            f"{args.empty_streak_to_stop} consecutive empty rings -- the galaxy's true edge was not "
            f"confirmed. Re-run with a larger --max-ring if these shape parameters really do produce "
            f"a galaxy this large."
        )


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
                         help="The kind of phenomenon to generate. Omit to pick uniformly at random.")

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
                              "galaxy position (those types have no placement columns of their own). Omit "
                              "to generate it unplaced/unlinked, as before.")

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
        `RoguePlanet`, `InterstellarComet`, or `AsteroidField`.

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

    raise ValueError(f"Unknown phenomenon type: {phenomenon_type!r}")


def run_phenomenon(args):
    """
    Generates one exotic phenomenon and saves it to the database.

    Args:
        args (argparse.Namespace): Validated arguments (`command ==
            "phenomenon"`).
    """
    phenomenon_type = args.type or random.choice(program_constants.PHENOMENON_TYPE_CHOICES)

    system_config = SystemConfig()
    system_config.MARKDOWN = args.markdown
    if args.num_orbits is not None:
        system_config.NUM_ORBITS = args.num_orbits

    phenomenon = generate_phenomenon(phenomenon_type, system_config, args.anchor_system, name=args.name)

    mysql_config = _db.mysql_config_from_args(args)
    phenomenon_id = _db.save_phenomenon(phenomenon, system_config, phenomenon_type, config=mysql_config,
                                         sector_id=args.sector_id)
    log.normal(f"Saved {TYPE_LABELS[phenomenon_type]} to the database (id={phenomenon_id}, "
               f"{mysql_config.database}@{mysql_config.host}:{mysql_config.port}).")


# ===========================================================================
# 6. Unified CLI
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

    return parser, {
        'system': system_parser,
        'sector': sector_parser,
        'galaxy': galaxy_parser,
        'plan': plan_parser,
        'phenomenon': phenomenon_parser,
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

    return args


_COMMAND_HANDLERS = {
    'system': run_system,
    'sector': run_sector,
    'galaxy': run_galaxy,
    'plan': run_plan,
    'phenomenon': run_phenomenon,
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
    log.configure(level, debug_file=(args.debug or None))
    log.debug("Command: %s, options: %s", args.command,
              {key: ("<withheld>" if "password" in key else value) for key, value in sorted(vars(args).items())})

    seed = secrets.randbits(128)
    random.seed(seed)
    log.debug(f"Seeded the random number generator with {seed} (cryptographically random; no --seed "
              f"option exists to reproduce this run).")

    _COMMAND_HANDLERS[args.command](args)


if __name__ == "__main__":
    main()
