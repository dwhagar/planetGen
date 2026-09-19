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
   (shell/slot address, `galactic_center_dist_ly`) and persisting that
   position and this sector's own Voronoi prism vertices. Three modes:
   `--shell` (batch), `--center-sector`/`--radius-pc` (local
   neighborhood), or neither (random start). `ensure_sector_generated`
   is this same per-sector logic exposed as a non-CLI, visit-triggered
   entry point against the galaxy skeleton section 4 builds.
4. Galaxy density skeleton (`plan`) -- the compact `galaxy_shape`/
   `galaxy_shell_band` summary `ensure_sector_generated` consults to
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
import multiprocessing
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
    enumerate_sectors_within_radius, galactic_radius_pc, provisional_sector_designation,
    sector_position_pc, shell_radius_pc, shell_sector_count,
)
from stellarObjects.galaxySkeleton import expected_system_count_at_density_1, find_shell_bands
from stellarObjects.nebulaData import Nebula
from stellarObjects.roguePlanetData import InterstellarComet, RoguePlanet
from stellarObjects.sectorGeometry import prism_vertices
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
    than it was worth, whereas a `galaxy` run (a whole shell, a
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
    generated sector itself, once per shell slot / neighborhood address).

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
        # count here -- run_shell_batch/run_local_neighborhood compute each
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


def generate_sector(args, galactic_center_dist_ly=None):
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

    Returns:
        tuple: `(sector_name, SpaceSector)` -- `sector_name` is
              `args.sector_name` or a freshly generated one (see
              `generate_sector_name`); the `SpaceSector` already has every
              system added (`SpaceSector.add_system`'s Hill-sphere-based
              random placement), plus a realistically sparse population of
              exotic phenomena (see `generate_sector_phenomena`). When
              `args.density` drove the system count (see below) and both
              that Poisson draw and `generate_sector_phenomena`'s own
              independent draws came back completely empty, one system is
              force-added anyway -- see the "guaranteed non-empty" note
              below.
    """
    sector_name = args.sector_name or generate_sector_name()
    sector = SpaceSector(name=sector_name)

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

    configs = build_sector_configs(args)

    systems = [
        StarSystem(system_config=cfg, galactic_center_dist_ly=galactic_center_dist_ly) for cfg in configs
    ]

    for system, cfg in zip(systems, configs):
        sector.add_system(system, system_config=cfg)

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

LARGE_SHELL_WARNING_THRESHOLD = 2000
"""int: `--shell K` requires `--limit` or `--yes` when shell `K` holds
more sector slots than this -- see `shell_sector_count`. Chosen so a
handful of inner shells (through roughly shell 20, ~5,400 slots -- close
to but slightly above this threshold at the boundary, so the check is
deliberately conservative) can be generated without an extra flag, while
anything large enough to take a real, unbounded amount of time and disk
space requires an explicit, deliberate choice."""


def add_galaxy_arguments(parser):
    """
    Adds the `galaxy` subcommand's own mode/placement options -- on top
    of whatever `add_shared_generation_options` already added -- to
    `parser`: the `--shell`/`--center-sector` mode-selection group,
    `--limit`, `--yes`, `--radius-pc`, `--max-shell`, `--min-start-density`,
    and the MySQL connection args.

    Args:
        parser (argparse.ArgumentParser): The parser to add options to.
    """
    mode_group = parser.add_mutually_exclusive_group(required=False)
    mode_group.add_argument('--shell', type=int, metavar='K',
                            help="Batch mode: generate every not-yet-generated sector slot in radial shell K "
                                 "(0-indexed; see docs/design/galaxy-coordinate-system.md section 3).")
    mode_group.add_argument('--center-sector', type=int, metavar='SECTOR_ID',
                            help="Local-neighborhood mode: generate every not-yet-generated sector slot "
                                 "within --radius-pc of the given, already galaxy-placed sector's own "
                                 "stored center.")

    parser.add_argument('--limit', type=int,
                        help="With --shell: generate only the first N not-yet-generated slots of the "
                             "shell, rather than the whole shell.")
    parser.add_argument('--yes', action='store_true',
                        help="With --shell: skip the confirmation normally required before generating a "
                             "shell whose total slot count exceeds LARGE_SHELL_WARNING_THRESHOLD.")
    parser.add_argument('--radius-pc', type=float,
                        help="With --center-sector: the neighborhood search radius, in parsecs. With "
                             "neither --shell nor --center-sector (random-start mode): overrides the "
                             "default 100 ly neighborhood radius around the randomly chosen starting "
                             "sector.")
    parser.add_argument('--max-shell', type=int,
                        help="With neither --shell nor --center-sector (random-start mode): the highest "
                             f"shell index the randomly chosen starting sector may land in. Defaults to "
                             f"the shell nearest a real Milky-Way-scale galaxy radius "
                             f"({program_constants.GALAXY_RADIUS_PC:,.0f} pc).")
    parser.add_argument('--min-start-density', type=float,
                        help="With neither --shell nor --center-sector (random-start mode): require the "
                             "randomly chosen starting sector's own real relative_density (the same "
                             "'expected' figure printed alongside each saved sector) to be at least this "
                             "value before accepting it -- e.g. 1.0 for at least as dense as the galaxy's "
                             "own real local density, 2.0 for twice that. Retried the same way an "
                             "already-occupied or otherwise-non-qualifying address is (see "
                             "RANDOM_START_MAX_PLACEMENT_ATTEMPTS); a high threshold combined with a large "
                             "--max-shell can take many more attempts to satisfy, since a volume-weighted "
                             "random draw favors the galaxy's own sparser outskirts. Cannot be combined "
                             "with --density/--num-systems (those override every position's density "
                             "uniformly, leaving no per-position value to compare against).")
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
    random_start = args.shell is None and args.center_sector is None

    if args.shell is not None and args.shell < 0:
        parser.error("--shell must be >= 0.")

    if args.center_sector is not None and args.radius_pc is None:
        parser.error("--center-sector requires --radius-pc.")
    if args.radius_pc is not None and args.shell is not None:
        parser.error("--radius-pc only applies to --center-sector or random-start mode (neither "
                     "--shell nor --center-sector), not --shell.")
    if args.radius_pc is not None and args.radius_pc <= 0:
        parser.error("--radius-pc must be a positive number.")

    if args.limit is not None and args.shell is None:
        parser.error("--limit only applies to --shell.")
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be a positive integer.")
    if args.yes and args.shell is None:
        parser.error("--yes only applies to --shell.")

    if args.max_shell is not None and not random_start:
        parser.error("--max-shell only applies to random-start mode (neither --shell nor "
                     "--center-sector).")
    if args.max_shell is not None and args.max_shell < 0:
        parser.error("--max-shell must be >= 0.")

    if args.min_start_density is not None and not random_start:
        parser.error("--min-start-density only applies to random-start mode (neither --shell nor "
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


class _BatchDensity:
    """
    Resolves each sector's own `--density` from the galaxy skeleton's real
    position-based `relative_density`, for `galaxy` mode's batch/local-
    neighborhood/random-start generation -- the same mechanism
    `ensure_sector_generated` already uses to pick a single lazily-
    generated sector's density (see its own docstring), applied here
    across a whole run instead of every sector in a batch quietly sharing
    one flat CLI value/default regardless of where it actually sits in the
    spiral structure `generate.py plan` computed.

    A no-op (`resolve` returns `args` unchanged) whenever the operator
    explicitly passed `--density`/`--num-systems` -- an explicit flag is
    still an intentional, uniform override for the whole run, not
    something this should second-guess. Only takes over for the "neither
    given" case `validate_shared_generation_args` now leaves both `None`
    for in `galaxy` mode specifically (see that function's own comment) --
    `sector` mode (no galaxy position to compute a density from at all)
    still gets its flat default of 10 there, unaffected.

    In that "neither given" case, `resolve` also applies the same
    qualification gate `ensure_sector_generated` already applies to a
    single lazily-generated sector (`galaxy_shell_band`'s stored candidate
    bands, then the exact `predicted_star_count`) -- returning `None`
    rather than a resolved namespace for a slot that doesn't qualify. Prior
    to this, `run_shell_batch`/`run_local_neighborhood`/`run_random_start`
    generated and saved a real (all-but-certainly-empty, 0-system,
    0-phenomena) sector row for *every* not-yet-occupied slot regardless of
    how far below the 1-star-per-sector threshold its own position's
    density fell -- most of a realistic galaxy's volume, off the spiral
    arms/disk plane -- so a batch or neighborhood run could come back
    "full of empty sectors" even though the single-sector lazy path
    (`ensure_sector_generated`) never has that problem. Skipping here
    instead keeps that same guarantee for every generation path.

    Fetches the stored skeleton (`generate.py plan`'s own output) once, on
    first use, and reuses it for the rest of the run -- a shell/local-
    neighborhood/random-start batch can mean thousands of sectors, and the
    skeleton itself never changes mid-run. Stored candidate bands are
    likewise fetched once per shell and cached (`_bands_cache`), since a
    shell's own band never changes mid-run either.
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

    def _get_bands(self, shell_index):
        if shell_index not in self._bands_cache:
            conn = _db.get_connection(self._config)
            try:
                self._bands_cache[shell_index] = _db.get_galaxy_shell_bands(conn, shell_index)
            finally:
                conn.close()
        return self._bands_cache[shell_index]

    def resolve(self, args, shell_index, slot_index, position_pc):
        """
        Args:
            args (argparse.Namespace): The `galaxy` subcommand's own parsed
                (and validated) arguments.
            shell_index (int): This sector's shell index.
            slot_index (int): This sector's slot index within that shell --
                needed (alongside `position_pc`) to check the stored
                candidate bands.
            position_pc (tuple): This one sector's `(x, y, z)` galaxy-frame
                position, in parsecs.

        Returns:
            argparse.Namespace or None: `args` itself when a density/count
                was given explicitly (no gating applied -- see this class's
                own docstring). Otherwise, `None` if this position doesn't
                qualify (outside every stored candidate band, or its own
                exact `predicted_star_count < 1.0`) -- the caller should
                skip this slot entirely rather than generate and save a
                sector for it. Otherwise a fresh copy (never mutates the
                shared `args`, the same reasoning `generate_sector`'s own
                `--density` handling already follows) with `.density` set
                to this sector's own `relative_density` and `.num_systems`
                cleared back to `None`, matching `ensure_sector_generated`'s
                own convention exactly.
        """
        if args.density is not None or args.num_systems is not None:
            return args
        skeleton = self._get_skeleton()

        bands = self._get_bands(shell_index)
        if not any(lo <= slot_index <= hi for lo, hi in bands):
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
    The (uniform, per this design's scope) sector edge length used for
    every shell/enumeration computation, in parsecs -- derived from
    `program_constants.DEFAULT_SECTOR_EDGE_LY`, the same fixed edge length
    every sector `SpaceSector` produces today (sector edge length is not
    currently an exposed CLI option anywhere in this package).

    Returns:
        float: The sector edge length, in parsecs.
    """
    return ly_to_pc(program_constants.DEFAULT_SECTOR_EDGE_LY)


def generate_and_save_sector_at(args, shell_index, shell_slot_index, position_pc, edge_pc):
    """
    Generates one sector via `generate_sector` at the given galaxy-frame
    position and shell address, and saves it -- the one per-sector unit
    of work both `run_shell_batch` and `run_local_neighborhood` repeat.

    Args:
        args (argparse.Namespace): Parsed arguments (see `add_galaxy_arguments`).
        shell_index (int): This sector's shell index.
        shell_slot_index (int): This sector's slot index within that shell.
        position_pc (tuple): `(x, y, z)` in parsecs, this sector's
                             galaxy-frame center (from `sector_position_pc`
                             or an `enumerate_sectors_within_radius`
                             result).
        edge_pc (float): The sector edge length, in parsecs (`_edge_pc`) --
                         threaded in rather than recomputed, since both
                         callers already have it.

    Returns:
        tuple: `(sector_id, sector_name, sector)` of the newly saved sector
            -- `sector_name` is read back from `sector.name` *after* saving
            (see the `Returns` note below), not the pre-save name
            `generate_sector` returned, since `stellarObjects._db`'s
            name-uniqueness machinery (v22) may rename it on save if it
            collides with something already in the database. `sector` (the
            full `SpaceSector`) is returned too so a caller can report its
            own per-type system/phenomenon breakdown without a separate
            query.
    """
    x, y, z = position_pc
    radius_pc = galactic_radius_pc(position_pc)
    galactic_center_dist_ly = pc_to_ly(radius_pc)

    _sector_name, sector = generate_sector(args, galactic_center_dist_ly=galactic_center_dist_ly)

    vertices_pc = prism_vertices(shell_index, shell_slot_index, edge_pc)
    galaxy_position = {
        "center_x_pc": x, "center_y_pc": y, "center_z_pc": z,
        "galactic_radius_pc": radius_pc,
        "shell_index": shell_index, "shell_slot_index": shell_slot_index,
        "vertices_pc": vertices_pc,
    }
    sector_id = _db.save_sector(sector, config=_db.mysql_config_from_args(args), galaxy_position=galaxy_position)
    # sector.name, not the discarded _sector_name above -- save_sector may
    # have just renamed it (a collision with an already-saved sector).
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
            onto the returned namespace's `mysql_*` attributes so
            `generate_and_save_sector_at`'s own
            `_db.mysql_config_from_args(args)` call works on this
            synthetic namespace exactly like it does on a real parsed
            result. Defaults to `DEFAULT_MYSQL_CONFIG`.

    Returns:
        argparse.Namespace: Every shared generation option at its
            documented default (`num_systems` resolved to `10`, since
            neither `--density` nor `--num-systems` was "given") -- a
            caller overriding density-driven generation (as
            `ensure_sector_generated` does) should set `args.density`
            and clear `args.num_systems` back to `None` before calling
            `generate_sector`/`generate_and_save_sector_at`.
    """
    parser = argparse.ArgumentParser(prefix_chars='-+')
    add_shared_generation_options(parser)
    args = parser.parse_args([])
    validate_shared_generation_args(args, parser)

    # Shaped like the `sector` subcommand's own output -- see that
    # function's docstring for why these fields must exist even though
    # their value is never used here.
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


def ensure_sector_generated(shell_index, shell_slot_index, config=None):
    """
    The galaxy map's "recalculate on visit" entry point: returns the
    sector already generated at this address if one exists; otherwise
    uses the stored skeleton (see section 4, `build_skeleton`) to decide,
    cheaply and exactly, whether this address is even worth generating,
    and -- if so -- generates and saves it on the spot. This is the
    backend half of "recompute whatever's needed as soon as a sector is
    visited"; wiring an actual UI/game-loop trigger to call this is
    separate, future work.

    A sector's own `relative_density` at its position (from the stored
    skeleton) becomes the `--density` multiplier `generate_sector` uses,
    rather than a uniform default -- so a bulge sector and a sparse
    outer-disk sector generate proportionally different system counts,
    matching what the skeleton itself predicted when deciding this address
    was worth visiting at all.

    A concurrent visit to the same never-before-generated address is
    possible (this isn't a single-process batch script) -- handled via
    `sectors`'s own `UNIQUE (shell_index, shell_slot_index)` constraint
    (schema.sql's "v8" note): whichever caller's `INSERT` loses the race
    gets `pymysql.err.IntegrityError` back from `generate_and_save_sector_at`,
    caught here and turned into "return what the other call just created"
    rather than a crash or a duplicate row.

    Args:
        shell_index (int): This sector's shell index.
        shell_slot_index (int): This sector's slot index within that
                                shell.
        config (MySQLConfig, optional): Connection parameters. Defaults
                                        to `DEFAULT_MYSQL_CONFIG`.

    Returns:
        dict: `created` (bool -- whether this call generated a new sector,
              as opposed to finding an existing one or confirming this
              address holds nothing), `qualifies` (bool -- whether this
              address holds, or would hold, any content at all),
              `sector_id` (int or `None`), `sector_name` (str or `None`,
              only set when `created` is `True`).

    Raises:
        RuntimeError: If the galaxy's skeleton has never been built --
                     run `generate.py plan` first.
    """
    conn = _db.get_connection(config)
    try:
        existing_id = _db.get_sector_id_at(conn, shell_index, shell_slot_index)
        if existing_id is not None:
            return {"created": False, "qualifies": True, "sector_id": existing_id, "sector_name": None}

        skeleton = _db.get_galaxy_shape(conn)
        if skeleton is None:
            raise RuntimeError(
                "The galaxy's skeleton has never been built (no galaxy_shape row) -- run "
                "'generate.py plan' first."
            )

        bands = _db.get_galaxy_shell_bands(conn, shell_index)
    finally:
        conn.close()

    if not any(lo <= shell_slot_index <= hi for lo, hi in bands):
        # Outside every stored candidate band -- galaxySkeleton.find_shell_bands's
        # own bound is exact, not sampled, so this is a certain "no",
        # not a heuristic one; no need to fall through to the exact check.
        return {"created": False, "qualifies": False, "sector_id": None, "sector_name": None}

    position_pc = sector_position_pc(shell_index, shell_slot_index, skeleton.edge_pc)
    density = relative_density(position_pc, skeleton.shape)
    star_count = predicted_star_count(position_pc, skeleton.shape, skeleton.expected_system_count_at_density_1)
    if star_count < 1.0:
        # Inside the shell's candidate band (a safe superset) but this
        # particular slot's own theta didn't clear the exact threshold --
        # see galaxySkeleton's own module docstring for why that's expected.
        return {"created": False, "qualifies": False, "sector_id": None, "sector_name": None}

    args = _default_generation_args(config=config)
    args.density = density
    args.num_systems = None

    try:
        sector_id, sector_name, _sector = generate_and_save_sector_at(
            args, shell_index, shell_slot_index, position_pc, skeleton.edge_pc,
        )
    except pymysql.err.IntegrityError:
        conn = _db.get_connection(config)
        try:
            existing_id = _db.get_sector_id_at(conn, shell_index, shell_slot_index)
        finally:
            conn.close()
        if existing_id is None:
            raise
        return {"created": False, "qualifies": True, "sector_id": existing_id, "sector_name": None}

    return {"created": True, "qualifies": True, "sector_id": sector_id, "sector_name": sector_name}


def run_shell_batch(args, edge_pc, progress):
    """
    Batch mode: generates every not-yet-generated sector slot in shell
    `args.shell` (up to `args.limit`, if given) that qualifies -- when
    neither `--density` nor `--num-systems` was given, a slot whose own
    position falls below the 1-star-per-sector threshold is skipped
    entirely (see `_BatchDensity.resolve`) rather than saved as an
    all-but-certainly-empty sector.

    Args:
        args (argparse.Namespace): Parsed arguments; `args.shell` must
            not be `None`.
        edge_pc (float): The sector edge length, in parsecs (`_edge_pc`).
        progress (rich.progress.Progress): `run_galaxy`'s shared progress
            display -- an outer "Sectors" task is added to it here (total
            = however many not-yet-generated slots this batch will
            actually generate) and advanced once per sector. Every status
            line below is logged via `log.normal`/`log.debug`, which routes
            through `progress.console` rather than a raw stdout write --
            see `_generation_progress`'s own docstring for why that matters
            while a `Progress` is live.

    Raises:
        SystemExit: If the shell's total slot count exceeds
                   `LARGE_SHELL_WARNING_THRESHOLD` and neither `--limit`
                   nor `--yes` was given.
    """
    shell_index = args.shell
    total_slots = shell_sector_count(shell_index)

    if total_slots > LARGE_SHELL_WARNING_THRESHOLD and args.limit is None and not args.yes:
        log.error(
            f"Shell {shell_index} holds {total_slots} sector slots -- generating a whole shell this "
            f"large is likely impractical. Pass --limit N to generate only the first N not-yet-generated "
            f"slots, or --yes to confirm generating all {total_slots}."
        )
        raise SystemExit(1)

    mysql_config = _db.mysql_config_from_args(args)
    conn = _db.get_connection(mysql_config)
    try:
        occupied = _db.get_occupied_shell_slots(conn, [shell_index])
    finally:
        conn.close()

    batch_density = _BatchDensity(mysql_config)

    to_generate = total_slots - len(occupied)
    if args.limit is not None and (args.density is not None or args.num_systems is not None):
        # Only a safe cap when every not-yet-occupied slot is guaranteed to
        # actually generate (an explicit --density/--num-systems bypasses
        # _BatchDensity's own qualification gate entirely) -- otherwise an
        # unknown number of slots this scan reaches may be skipped for
        # falling below the star-count threshold before --limit many
        # sectors are actually generated, so the total can't be capped to
        # --limit in advance without the bar overrunning it.
        to_generate = min(to_generate, args.limit)
    outer_task = progress.add_task(f"Sectors (shell {shell_index})", total=max(to_generate, 0))

    generated = 0
    skipped = 0
    for slot_index in range(total_slots):
        if args.limit is not None and generated >= args.limit:
            break
        if (shell_index, slot_index) in occupied:
            continue

        position_pc = sector_position_pc(shell_index, slot_index, edge_pc)
        sector_args = batch_density.resolve(args, shell_index, slot_index, position_pc)
        if sector_args is None:
            # Below the 1-star-per-sector threshold (or outside every
            # stored candidate band) -- skip it entirely rather than save
            # an all-but-certainly-empty sector, matching
            # ensure_sector_generated's own gating (see _BatchDensity).
            log.debug(f"Shell {shell_index} slot {slot_index}: skipped (below the 1-star-per-sector "
                      f"threshold, or outside every stored candidate band)")
            skipped += 1
            progress.update(outer_task, advance=1)
            continue

        log.debug(f"Shell {shell_index} slot {slot_index}: generating (density={sector_args.density})")
        sector_id, sector_name, sector = generate_and_save_sector_at(
            sector_args, shell_index, slot_index, position_pc, edge_pc,
        )
        generated += 1
        progress.update(outer_task, advance=1)
        designation = provisional_sector_designation(
            shell_index, slot_index, edge_pc, program_constants.DEFAULT_SECTOR_EDGE_LY,
        )
        log.normal(
            f"Saved sector '{sector_name}' [{designation}] at shell {shell_index} slot {slot_index} "
            f"(sector_id={sector_id})."
        )
        log.normal(sector_generation_summary_lines(sector, sector_args))

    already_existed = len(occupied)
    skip_note = f", {skipped} skipped (below the star-count threshold)" if skipped else ""
    log.normal(
        f"Generated {generated} new sector(s) in shell {shell_index} "
        f"({total_slots} total slots, {already_existed} already existed{skip_note})."
    )


def run_local_neighborhood(args, edge_pc, progress):
    """
    Local-neighborhood mode: generates every not-yet-generated, qualifying
    sector slot within `args.radius_pc` parsecs of `args.center_sector`'s
    own stored galaxy-frame center -- see `run_shell_batch`'s own
    docstring for what "qualifying" means and when it applies.

    Args:
        args (argparse.Namespace): Parsed arguments;
            `args.center_sector`/`args.radius_pc` must not be `None`.
        edge_pc (float): The sector edge length, in parsecs (`_edge_pc`).
        progress (rich.progress.Progress): `run_galaxy`'s shared progress
            display -- see `run_shell_batch`'s own `progress` docstring;
            same role here, an outer "Sectors" task over this mode's own
            not-yet-generated candidate count. `run_random_start` also
            calls this directly (after generating its own seed sector),
            adding a second "Sectors" task to the same `Progress` rather
            than a `run_shell_batch`-shaped one.

    Raises:
        SystemExit: If `args.center_sector` doesn't exist, or exists but
                   has never been placed in a galaxy (its galaxy-position
                   columns are NULL -- e.g. a sector generated via the
                   `sector` subcommand rather than `galaxy`).
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
            f"generated via 'generate.py galaxy' itself, not 'generate.py sector'). Use --shell to "
            f"generate placed sectors from scratch instead."
        )
        raise SystemExit(1)

    center = (
        center_position["center_x_pc"], center_position["center_y_pc"], center_position["center_z_pc"],
    )
    candidates = list(enumerate_sectors_within_radius(center, args.radius_pc, edge_pc))

    candidate_shells = sorted({shell_index for shell_index, _slot, _x, _y, _z, _dist in candidates})
    mysql_config = _db.mysql_config_from_args(args)
    conn = _db.get_connection(mysql_config)
    try:
        occupied = _db.get_occupied_shell_slots(conn, candidate_shells)
    finally:
        conn.close()

    batch_density = _BatchDensity(mysql_config)

    to_generate = sum(
        1 for shell_index, slot_index, _x, _y, _z, _dist in candidates
        if (shell_index, slot_index) not in occupied
    )
    outer_task = progress.add_task("Sectors (local neighborhood)", total=to_generate)

    generated = 0
    skipped = 0
    already_existed = 0
    for shell_index, slot_index, x, y, z, distance_pc in candidates:
        if (shell_index, slot_index) in occupied:
            already_existed += 1
            continue

        sector_args = batch_density.resolve(args, shell_index, slot_index, (x, y, z))
        if sector_args is None:
            # Below the 1-star-per-sector threshold (or outside every
            # stored candidate band) -- skip it entirely rather than save
            # an all-but-certainly-empty sector, matching
            # ensure_sector_generated's own gating (see _BatchDensity).
            log.debug(f"Shell {shell_index} slot {slot_index}: skipped (below the 1-star-per-sector "
                      f"threshold, or outside every stored candidate band)")
            skipped += 1
            progress.update(outer_task, advance=1)
            continue

        log.debug(f"Shell {shell_index} slot {slot_index}: generating (density={sector_args.density})")
        sector_id, sector_name, sector = generate_and_save_sector_at(
            sector_args, shell_index, slot_index, (x, y, z), edge_pc,
        )
        generated += 1
        progress.update(outer_task, advance=1)
        designation = provisional_sector_designation(
            shell_index, slot_index, edge_pc, program_constants.DEFAULT_SECTOR_EDGE_LY,
        )
        log.normal(
            f"Saved sector '{sector_name}' [{designation}] at shell {shell_index} slot {slot_index}, "
            f"{distance_pc:.2f} pc from sector_id={args.center_sector} (sector_id={sector_id})."
        )
        log.normal(sector_generation_summary_lines(sector, sector_args))

    skip_note = f", {skipped} skipped (below the star-count threshold)" if skipped else ""
    log.normal(
        f"Generated {generated} new sector(s) within {args.radius_pc} pc of sector_id={args.center_sector} "
        f"({len(candidates)} candidate slot(s) found, {already_existed} already existed{skip_note})."
    )


def generate_sector_neighborhood(center_sector_id, radius_ly=None, config=None):
    """
    Non-CLI counterpart to `run_local_neighborhood`'s core logic -- for a
    caller with no `argparse.Namespace` of its own (the admin web UI's
    "generate more sectors around this one"
    action, `html/api/routes.py`'s `generate_sector_neighborhood_route`),
    rather than the `galaxy` subcommand's `--center-sector` mode. Same
    underlying work (`_db.get_sector_galaxy_position`,
    `enumerate_sectors_within_radius`, `_db.get_occupied_shell_slots`,
    `generate_and_save_sector_at`), a plain result dict instead of prints,
    and a catchable `ValueError` instead of `SystemExit` for an invalid/
    unplaced sector -- there's no CLI here for `SystemExit` to exit out of.

    The default 100 ly radius is genuinely large relative to one sector's
    edge (`program_constants.DEFAULT_SECTOR_EDGE_LY`, 11.5 ly) -- its
    sphere holds on the order of **2,000-3,000 candidate sector slots**
    (confirmed by measurement, not just geometry: `(4/3)*pi*100**3 /
    11.5**3 ≈ 2750`), same as `galaxy`'s own random-start mode already
    generates today from the CLI. Called with no `radius_ly` override,
    this can take minutes to hours depending on the server and how many
    of those slots are already occupied -- every caller (the CLI, and
    especially `generate_sector_neighborhood_route`'s web-triggered,
    synchronous-request version of this) needs to account for that, not
    assume "generate a neighborhood" is a quick call.

    Args:
        center_sector_id (int): The already galaxy-placed sector to
            generate a neighborhood around.
        radius_ly (float, optional): Defaults to
            `program_constants.RANDOM_START_NEIGHBORHOOD_RADIUS_LY`
            (100 ly) -- the same default radius `galaxy`'s own
            random-start mode uses.
        config (MySQLConfig, optional): Connection parameters. Defaults
                                        to `DEFAULT_MYSQL_CONFIG`.

    Returns:
        dict: `generated` (int -- newly created sectors), `already_existed`
            (int -- candidate slots that already had a sector), `skipped`
            (int -- candidate slots below the 1-star-per-sector threshold,
            per `_BatchDensity`'s own gating), `candidates` (int -- total
            slots within the radius).

    Raises:
        ValueError: If `center_sector_id` doesn't exist, or exists but has
                   never been placed in a galaxy (its galaxy-position
                   columns are NULL -- e.g. a sector generated via the
                   `sector` subcommand rather than `galaxy`).
        RuntimeError: If the galaxy's skeleton has never been built --
                     see `_BatchDensity`/run `generate.py plan` first.
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
    candidates = list(enumerate_sectors_within_radius(center, radius_pc, edge_pc))

    candidate_shells = sorted({shell_index for shell_index, _slot, _x, _y, _z, _dist in candidates})
    conn = _db.get_connection(config)
    try:
        occupied = _db.get_occupied_shell_slots(conn, candidate_shells)
    finally:
        conn.close()

    args = _default_generation_args(config=config)
    # Drives each sector's own system count from the galaxy skeleton's real
    # position-based relative_density (_BatchDensity), the same way
    # run_local_neighborhood/ensure_sector_generated already do -- without
    # this, `_default_generation_args`'s own flat num_systems=10 default
    # would apply uniformly regardless of local density.
    args.density = None
    args.num_systems = None

    batch_density = _BatchDensity(config)

    generated = 0
    skipped = 0
    already_existed = 0
    for shell_index, slot_index, x, y, z, _distance_pc in candidates:
        if (shell_index, slot_index) in occupied:
            already_existed += 1
            continue
        sector_args = batch_density.resolve(args, shell_index, slot_index, (x, y, z))
        if sector_args is None:
            skipped += 1
            continue
        generate_and_save_sector_at(sector_args, shell_index, slot_index, (x, y, z), edge_pc)
        generated += 1

    return {
        "generated": generated,
        "already_existed": already_existed,
        "skipped": skipped,
        "candidates": len(candidates),
    }


def _pick_random_shell_index(max_shell_index, edge_pc):
    """
    Picks a shell index in `[0, max_shell_index]`, weighted by that shell's
    own volume (~proportional to `(shell_index + 0.5) ** 2`, since every
    shell shares the same radial thickness `edge_pc`) rather than uniformly
    across shell indices -- the latter would hugely overrepresent the
    sparse galactic core, where a shell holds far fewer sector slots than
    one at the same index gap farther out (see `shell_sector_count`).
    Sampling a random radius `r = r_max * u^(1/3)` (`u` uniform in `[0,
    1)`) is the standard closed-form way to draw a point uniformly *by
    volume* within a sphere -- the same cube-root correction
    `spaceSector._random_point_in_annulus` uses for its own inner-shell-vs-
    outer-shell volume bias -- then mapping that radius to the nearest
    shell index inverts `shell_radius_pc`.

    Used by `run_random_start` to choose where its "random sector" lands:
    a uniformly random point within the galaxy's own volume, not an
    arbitrary or core-biased one.

    Args:
        max_shell_index (int): The highest shell index that may be chosen.
        edge_pc (float): The sector edge length, in parsecs.

    Returns:
        int: The chosen shell index, in `[0, max_shell_index]`.
    """
    r_max = shell_radius_pc(max_shell_index, edge_pc)
    r = r_max * random.random() ** (1 / 3)
    shell_index = round(r / edge_pc - 0.5)
    return max(0, min(max_shell_index, shell_index))


def run_random_start(args, edge_pc, progress):
    """
    Random-start mode (no `--shell`/`--center-sector` given): picks a
    random, not-yet-occupied sector address somewhere within a real
    Milky-Way-scale galaxy, generates it, then falls straight through to
    `run_local_neighborhood`'s own logic to generate every not-yet-
    generated sector within `args.radius_pc` of it too -- by default,
    `program_constants.RANDOM_START_NEIGHBORHOOD_RADIUS_LY` (100 ly),
    converted to parsecs, in every direction, per this feature's own
    request.

    The address is chosen by `_pick_random_shell_index` (a uniformly
    random point by volume within the galaxy's sphere, up to `--max-shell`
    or `program_constants.GALAXY_RADIUS_PC`'s own nearest shell) plus a
    uniformly random slot within that shell -- retried, per
    `program_constants.RANDOM_START_MAX_PLACEMENT_ATTEMPTS`, whenever the
    drawn address already has a sector (overwhelmingly unlikely for any
    galaxy that isn't already almost entirely generated; see that
    constant's own docstring), or (when neither `--density` nor
    `--num-systems` was given) doesn't itself qualify per `_BatchDensity`'s
    own gating -- likely for any given draw, since most of a real galaxy's
    volume sits off the spiral arms/disk plane, but rare enough in
    aggregate across the whole sphere that a retry almost always lands on
    a qualifying address well within the attempt budget. `--min-start-density`
    tightens that same retry loop further: a qualifying address whose own
    `relative_density` still falls short of it is retried exactly like a
    non-qualifying one.

    Args:
        args (argparse.Namespace): Parsed arguments;
            `args.shell`/`args.center_sector` must both be `None`.
        edge_pc (float): The sector edge length, in parsecs (`_edge_pc`).
        progress (rich.progress.Progress): `run_galaxy`'s shared progress
            display, threaded straight through to `run_local_neighborhood`
            for the surrounding neighborhood's own "Sectors (local
            neighborhood)" task -- the seed sector generated below gets no
            task of its own (a single-sector 0-to-1 bar is done before it
            can even render a meaningful rate/ETA, so it only ever added
            noise, not a genuine progress display).

    Raises:
        SystemExit: If no unoccupied, qualifying address meeting
                   `args.min_start_density` (if given) could be found
                   within `program_constants.RANDOM_START_MAX_PLACEMENT_ATTEMPTS`
                   attempts.
    """
    max_shell_index = (
        args.max_shell if args.max_shell is not None
        else int(program_constants.GALAXY_RADIUS_PC / edge_pc - 0.5)
    )
    radius_pc = (
        args.radius_pc if args.radius_pc is not None
        else ly_to_pc(program_constants.RANDOM_START_NEIGHBORHOOD_RADIUS_LY)
    )

    mysql_config = _db.mysql_config_from_args(args)
    # `run_local_neighborhood` below builds its own separate `_BatchDensity`
    # for the surrounding neighborhood; this one is just for picking (and
    # generating) the seed sector itself.
    batch_density = _BatchDensity(mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        sector_args = None
        for _ in range(program_constants.RANDOM_START_MAX_PLACEMENT_ATTEMPTS):
            shell_index = _pick_random_shell_index(max_shell_index, edge_pc)
            slot_index = random.randrange(shell_sector_count(shell_index))
            if _db.get_sector_id_at(conn, shell_index, slot_index) is not None:
                continue
            position_pc = sector_position_pc(shell_index, slot_index, edge_pc)
            # Same qualification gate `run_shell_batch`/`run_local_neighborhood`
            # apply -- an unoccupied address whose own position doesn't
            # clear the 1-star-per-sector threshold is retried exactly like
            # an already-occupied one, rather than generated as an
            # all-but-certainly-empty seed sector.
            sector_args = batch_density.resolve(args, shell_index, slot_index, position_pc)
            if sector_args is None:
                continue
            if args.min_start_density is not None and sector_args.density < args.min_start_density:
                # Qualifies (>= 1 predicted star), but not dense enough to
                # satisfy the operator's own --min-start-density -- retried
                # the same way as any other rejected draw.
                continue
            break
        else:
            density_note = (
                f", meeting --min-start-density {args.min_start_density} (try lowering it or --max-shell "
                f"-- a high density threshold combined with a large --max-shell means most volume-weighted "
                f"draws land in the galaxy's own sparser outskirts, far short of it)"
                if args.min_start_density is not None else ""
            )
            log.error(
                f"Could not find an unoccupied, qualifying sector address within {max_shell_index} "
                f"shells after {program_constants.RANDOM_START_MAX_PLACEMENT_ATTEMPTS} attempts{density_note} "
                f"-- this galaxy may already be almost entirely generated within that range, or that range "
                f"may hold too little real stellar density; try a larger --max-shell."
            )
            raise SystemExit(1)
    finally:
        conn.close()

    sector_id, sector_name, sector = generate_and_save_sector_at(
        sector_args, shell_index, slot_index, position_pc, edge_pc,
    )
    designation = provisional_sector_designation(
        shell_index, slot_index, edge_pc, program_constants.DEFAULT_SECTOR_EDGE_LY,
    )
    log.normal(
        f"Saved random starting sector '{sector_name}' [{designation}] at shell {shell_index} slot "
        f"{slot_index} (sector_id={sector_id})."
    )
    log.normal(sector_generation_summary_lines(sector, sector_args))

    args.center_sector = sector_id
    args.radius_pc = radius_pc
    run_local_neighborhood(args, edge_pc, progress)


def run_galaxy(args):
    """
    Dispatches to shell-batch, local-neighborhood, or random-start mode.

    Owns the one `rich.progress.Progress` display shared across whichever
    mode runs -- each mode adds its own "Sectors" task to it (see
    `run_shell_batch`/`run_local_neighborhood`/`run_random_start`'s own
    `progress` docstrings).

    Args:
        args (argparse.Namespace): Validated arguments (`command ==
            "galaxy"`).
    """
    edge_pc = _edge_pc()

    with _generation_progress() as progress:
        log.set_console(progress.console)
        try:
            if args.shell is not None:
                run_shell_batch(args, edge_pc, progress)
            elif args.center_sector is not None:
                run_local_neighborhood(args, edge_pc, progress)
            else:
                run_random_start(args, edge_pc, progress)
        finally:
            log.reset_console()


# ===========================================================================
# 4. Galaxy density skeleton
# ===========================================================================

DEFAULT_EMPTY_STREAK_TO_STOP = 50
"""int: How many consecutive empty shells confirms the galaxy's true edge
has been reached -- see `galaxySkeleton`'s own docstring on why the
qualifying region is expected to shrink monotonically outward, plus a
safety margin against a razor-thin band a coarse scan could miss for one
single shell."""

DEFAULT_MAX_SHELL = 100000
"""int: Hard cap on how many shells `plan` will ever scan, regardless of
`--empty-streak-to-stop` -- protects against a pathological shape
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


def add_plan_arguments(parser):
    """
    Adds every option the `plan` subcommand accepts (besides `--version`)
    to `parser` -- the galaxy shape parameter group plus the scan/
    parallelism options and the MySQL connection args.

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
    add_logging_arguments(parser)


def validate_plan_args(args, parser):
    """
    Validates `add_plan_arguments`'s own options, calling `parser.error`
    (which exits) on the first problem found, then defaults `--workers`/
    `--chunk-size` from `os.cpu_count()`.

    Args:
        args (argparse.Namespace): Parsed arguments.
        parser (argparse.ArgumentParser): The parser to raise errors
                                          through (so the caller's own
                                          `--help`/usage text is shown).
    """
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


def build_skeleton(args):
    """
    Runs the full parallel skeleton build and persists it: `galaxy_shape`
    (the galaxy's shape parameters, calibration constant, sector edge
    length, and outer edge -- one singleton row) and `galaxy_shell_band`
    (one row per contiguous *candidate* slot-index band per shell -- a
    safe, cheap-to-compute superset of where a shell's qualifying sectors
    could be, not an exact per-sector list; see
    `galaxySkeleton.find_shell_bands`'s own docstring).

    No sector content, vertices, or even individual sector addresses are
    generated or stored here -- that stays lazy, happening only when a
    sector is actually visited (`ensure_sector_generated`). Re-running
    this always replaces the entire skeleton wholesale -- there is no
    incremental update, since a full build is already well under a
    minute even at real Milky-Way scale (per-shell work here is a
    closed-form calculation, not a per-sector scan).

    Shells are scanned outward in order (needed to detect the galaxy's
    true edge -- a run of consecutive empty shells), but each shell's own
    band-finding is fully independent of every other shell's, so shells
    are dispatched to a `multiprocessing.Pool` in fixed-size chunks
    (`--chunk-size`, default scaled to `--workers`).

    Args:
        args (argparse.Namespace): Parsed arguments.

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

    log.normal(
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


def run_plan(args):
    """
    Builds and persists the galaxy's density skeleton.

    Args:
        args (argparse.Namespace): Validated arguments (`command ==
            "plan"`).
    """
    summary = build_skeleton(args)
    log.normal(
        f"Skeleton built in {summary['elapsed_s']:.2f}s: scanned {summary['shells_scanned']} shells, "
        f"outer edge = shell {summary['outer_shell_index']}, {summary['total_bands']} band(s) stored, "
        f"~{summary['total_candidate_slots']:,} candidate sector slots."
    )
    if not summary["edge_confirmed"]:
        log.normal(
            f"WARNING: reached --max-shell ({args.max_shell}) without a run of "
            f"{args.empty_streak_to_stop} consecutive empty shells -- the galaxy's true edge was not "
            f"confirmed. Re-run with a larger --max-shell if these shape parameters really do produce "
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

    seed = secrets.randbits(128)
    random.seed(seed)
    log.debug(f"Seeded the random number generator with {seed} (cryptographically random; no --seed "
              f"option exists to reproduce this run).")

    _COMMAND_HANDLERS[args.command](args)


if __name__ == "__main__":
    main()
