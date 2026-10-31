# planetgen/cli/generate.py

"""
Generate
========

The single entry point for every generator, dispatched by subcommand
(run from the checkout's `src/` as `python3 -m planetgen.cli.generate`,
or as the installed `planetgen` command):

    planetgen system [options]      -- one star system
    planetgen sector [options]      -- one or more independent sectors
    planetgen sectors [options]     -- sectors of the planned galaxy, in one of the modes
    planetgen galaxy [options]      -- many sectors placed as one galaxy
    planetgen plan [options]        -- the galaxy's density skeleton
    planetgen link --since TIME     -- link the sectors saved since TIME to their neighbours
    planetgen backfill --since TIME -- the bright-star backfill out from the sectors saved since TIME
    planetgen settle --since TIME   -- save the paths of the sectors saved since TIME and of those around them
    planetgen scatter --KIND        -- one scatter (mass, luminosity, phenomena, nebula-field, band)
    planetgen fill-map-luminosity   -- fill the map luminosity of the scattered phenomena
    planetgen rescore-habitability  -- re-derive the stored habitability numbers
    planetgen fill-nebula-blocks    -- make the Galaxy Map's nebula block totals again
    planetgen fill-block-totals     -- make the Galaxy Map's block light totals again
    planetgen clear-scatter         -- empty one kind of scatter, so a scatter can draw it again
    planetgen clean-nebula-field    -- count (or with --apply delete) the stored field clouds the keep rate rejects
    planetgen phenomenon [options]  -- one exotic stellar phenomenon
    planetgen population [options]  -- species, civilizations and polities
    planetgen territories [options] -- which polity owns which system
    planetgen check-math            -- the math check (update.sh runs it; so does the admin page)
    planetgen versions              -- the version keys each update recorded
    planetgen fingerprint [options] -- a digest of the generated content (GEN.58)

Run `planetgen <command> --help` for that command's own full option
list. This module is the command line: every command's options and
their validation, the parser, the run history, and `main`. Each
command's work is in `planetgen.generation.run_<command>`.
"""

import argparse
import getpass
import json
import math
import os
import secrets
import signal
import sys
import threading
import time

import pymysql
from planetgen.cli import benchmark
from planetgen.cli import checkdb as checkdb_cli, queue as queue_cli, warm as warm_cli

from planetgen.queue import commands as queueCommands
from planetgen.db import check as dbCheck, fingerprint, parity, repair as dbRepair, store, table_rebuild
from planetgen.db import fingerprint_tree as fingerprintTree
from planetgen.admin import activity_log, updateguard
from planetgen.galaxy.span import Span, SpanError, parse_range
from planetgen.generation import directives, limits, luminosity_floor, prevalence
from planetgen.galaxy import bring_to_epoch as bringToEpoch, drill, seed as galaxySeed, spherical, version_check, version_key
from planetgen import tuning as program_constants
from planetgen.util import log
from planetgen._version import VersionAction, version_banner
from planetgen.galaxy import region_run, scatter_limits
from planetgen.galaxy import dials
from planetgen.galaxy import shape_draw as shapeDraw
from planetgen.galaxy.skeleton import DEFAULT_REACH_SHELL
from planetgen.generation.star import STAR_TYPE_PATTERN
from planetgen.generation import star_population
from planetgen.generation.star_population import bright_star_fraction
from planetgen.generation import stages
from planetgen.generation import finish_run as finishRun
from planetgen.generation import job as generationJob
from planetgen.generation import run_common
from planetgen.generation import run_galaxy
from planetgen.generation import run_phenomenon
from planetgen.generation import run_plan
from planetgen.generation import run_population
from planetgen.datasets import fetch
from planetgen.db import dataversion
from planetgen.generation import run_sector
from planetgen.generation import run_system
from planetgen.util import draw


class TristateAction(argparse.Action):
    """
    Sets `namespace.dest` to True when invoked as `+name`, or False when
    invoked as `-name`. Leaving the option off the command line leaves the
    `default` (None) in place, meaning "let the generator decide".
    """

    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, option_string.startswith('+'))


class SingleSystemOnlyAction(argparse.Action):
    """
    A forcing option (`+name`/`-name`) given to `sector` or `galaxy`: it
    applies only to a single system now (GEN.51), so the run stops with a
    message naming the option rather than argparse's bare "unrecognized
    arguments" -- which is what a queued or saved command line from before
    the change would otherwise get.
    """

    def __call__(self, parser, namespace, values, option_string=None):
        parser.error(f"{option_string} forces one system, so it only works with 'planetgen system' "
                     f"(and the one-off system page); sector and galaxy runs no longer take forcing options.")


def prevalence_setting(text):
    """
    `argparse` type for `--prevalence`: `FEATURE=PERCENT` (GEN.52), a
    feature from `prevalence.FEATURES` and a percentage of at least -100.

    Returns:
        tuple: `(feature, percent)`.
    """
    feature, sep, value = text.partition("=")
    feature = feature.strip().lower().replace("-", "_")
    if not sep or feature not in prevalence.FEATURES:
        raise argparse.ArgumentTypeError(
            f"expected FEATURE=PERCENT with FEATURE one of {', '.join(prevalence.FEATURES)}; got {text!r}")
    percent = finite_float(value.strip().rstrip("%"))
    if percent < prevalence.MIN_PERCENT:
        raise argparse.ArgumentTypeError(f"{feature}={value}: a prevalence can't be below -100%")
    return feature, percent


def compact_mass_limit(text):
    """`argparse` type of `--compact-min-mass` (GEN.195): "star" (use the stellar mass limit) or a number."""
    if text.strip().lower() == program_constants.COMPACT_MIN_MASS_STAR:
        return program_constants.COMPACT_MIN_MASS_STAR
    return finite_float(text)


def check_compact_mass_limit(args, parser):
    """Errors unless `--compact-min-mass` is one of `tuning.COMPACT_MIN_MASS_PRESETS` or "star"."""
    value = getattr(args, "compact_min_mass", None)
    if value is None or value == program_constants.COMPACT_MIN_MASS_STAR:
        return
    if not any(math.isclose(value, preset) for preset in program_constants.COMPACT_MIN_MASS_PRESETS):
        parser.error("--compact-min-mass must be one of "
                     + ", ".join(f"{preset:g}" for preset in program_constants.COMPACT_MIN_MASS_PRESETS)
                     + f", or {program_constants.COMPACT_MIN_MASS_STAR}.")


COMPACT_MIN_MASS_HELP = (
    "The lowest mass (solar masses) of a neutron star or black hole placed galaxy-wide: one of "
    + ", ".join(f"{preset:g}" for preset in program_constants.COMPACT_MIN_MASS_PRESETS)
    + ", or 'star' to use the stellar mass limit (--phenomenon-min-mass), which is also the default "
      "when this is left out. Lighter ones are drawn when their sector is filled. A lower limit places "
      "far more rows: 1 solar mass is about 1.2 billion rows (161 GB), 2 about 390 million (54 GB), "
      "4 about 220 million (30 GB), 6 about 200 million (28 GB). The central supermassive black hole "
      "or quasar is always created.")


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


def _validate_name(args, parser):
    """--name must fit the database (with room for the planets and moons
    named after it), checked here rather than failing at the save."""
    name = getattr(args, "name", None)
    if name is not None and len(name) > store.SYSTEM_NAME_MAX_LENGTH:
        parser.error(f"--name is {len(name)} characters; a system name can be at most "
                     f"{store.SYSTEM_NAME_MAX_LENGTH}.")


def add_logging_arguments(parser):
    """
    Adds the logging options every subcommand shares to `parser`:
    `--debug` and `--quiet`/`--silent`, and `--disregard-missing-data`
    (GEN.599: not logging, but a run-wide switch every subcommand takes, so a
    command that reads no downloaded data accepts it and ignores it).

    `--debug` alone enables debug-severity output to the console (stdout), the
    lines of the run's queue workers included; giving it a filename
    (`--debug FILE`) additionally mirrors that same output to `FILE`.
    Debug mode in `config.json` (or `PLANETGEN_DEBUG`) does the same for
    every command, and also writes the debug log. `--debug-areas` limits
    what is written (`planetgen.util.log.AREAS`); by default it is
    everything. `--quiet`/`--silent` (two spellings of the same flag) restrict
    output to errors only. `validate_logging_args` enforces the one
    combination that doesn't make sense: `--quiet`/`--silent` together with
    a filename-less `--debug`, which would have nowhere to send debug
    output once the console is silenced.

    Args:
        parser (argparse.ArgumentParser): The parser to add options to.
    """
    parser.add_argument('--debug', nargs='?', const='', default=None, metavar='FILE',
                        help="Log every decision the run makes, and why (every random draw, branch choice, stage, "
                             "unit, SQL statement and cache hit, workers included), with ISO timestamps, to stdout. "
                             "If FILE is given, also mirror that output to FILE. Debug mode in config.json does "
                             "this for every command.")
    parser.add_argument('--debug-areas', metavar='AREAS', default=None,
                        help="Limit --debug to some kinds of output: " + ", ".join(log.AREAS) + " (comma "
                             "separated; a leading - leaves one out, as in -draw). The default is all of them.")
    parser.add_argument('--record-file', metavar='FILE', default=None,
                        help="Also write every log line to FILE as a record: one per line, the time in UTC "
                             "(PLANETGEN_RECORD_FILE does the same).")
    parser.add_argument('--record-json', action='store_true',
                        help="Write the --record-file records as JSON lines.")
    parser.add_argument('--quiet', '--silent', dest='quiet', action='store_true',
                        help="Suppress all output except errors.")
    parser.add_argument('--disregard-missing-data', action='store_true',
                        help="Generate even when a downloaded data file is missing (refused by default), taking "
                             "each missing dataset's labelled fallback; the run says which. Setup downloads the files.")


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
    if args.debug_areas is not None:
        try:
            log.parse_areas(args.debug_areas)
        except ValueError as exc:
            parser.error(f"--debug-areas: {exc}")


_RING_LIMIT = f"{run_galaxy.LARGE_RUN_WARNING_THRESHOLD:,}"
"""str: The size past which a shell or block run needs `--limit` or `--yes`, as the help says it."""

STRICT_HELP = ("Stop instead of warning: a run the console would otherwise go ahead with after a "
               "warning (the database disk too small, a shell or block past "
               f"{_RING_LIMIT} sectors without --limit or --yes, an address outside the "
               "galaxy's outline, --min-habitable above a sector's drawn count, a forced body no "
               "system had room for) ends with an error and nothing generated.")


MAX_STARS_HELP = (f"The most stars a system may hold: the cap on the star-count draw (default {program_constants.MAX_STARS_DEFAULT}, "
                  f"at most {program_constants.MAX_STARS_CEILING}). Raise it to let a rare system hold more; "
                  f"--prevalence multi_star=PERCENT changes how often five or more stars come up.")


def add_system_arguments(parser):
    """
    Adds every option the `system` subcommand accepts (besides
    `--version`, which only the top-level parser offers) to `parser` --
    the tri-state `+name`/`-name` flags (see `TRISTATE_OPTIONS`),
    `--system-file`, the MySQL connection args, `--format`,
    `--star-type`, `--num-orbits`, `--name`, `--age`, the flavor-text
    overrides, and the logging options (see `add_logging_arguments`).

    Args:
        parser (argparse.ArgumentParser): The parser to add options to.
                                          Must accept `+`/`-` prefix chars
                                          (`prefix_chars='-+'`) for the
                                          tri-state options.
    """
    for name, _attr, description in run_system.TRISTATE_OPTIONS:
        parser.add_argument(f'-{name}', f'+{name}', dest=name, action=TristateAction,
                            nargs=0, default=None,
                            help=f"+{name} forces the system to have {description}; "
                                 f"-{name} forces the system to not have {description}.")

    # Load system options from a JSON file
    parser.add_argument('--system-file', '-f', type=str,
                        help="Load system generation options from a JSON file. Command-line options "
                             "override the values it sets.")

    # Database persistence
    store.add_mysql_connection_args(parser)

    parser.add_argument('--strict', action='store_true', help=STRICT_HELP)

    # Write the description instead of saving (a one-off system)
    parser.add_argument('--format', choices=run_system.OUTPUT_FORMATS, default=run_system.OUTPUT_FORMATS[0],
                        help="With --output, how to write the system: plain text (the default), Markdown "
                             "or the Document as JSON.")
    parser.add_argument('--output', '-o', type=str, metavar='FILE',
                        help="Write the system's description (in --format) to FILE "
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

    # Number of stars (GEN.267)
    parser.add_argument('--stars', type=int,
                        help="Force the number of stars in the system, from 1 up to the cap (--max-stars). "
                             "Refused when it is over the cap or contradicts +binary_system, -binary_system "
                             "or +wide_binary.")
    parser.add_argument('--max-stars', type=int, default=None, metavar='N',
                        help=MAX_STARS_HELP)

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
    if args.planets is False and (args.moons or args.max_planets or args.habitable_world or args.asteroid_belt):
        parser.error("-planets cannot be combined with +moons, +max_planets, +habitable_world, or +asteroid_belt.")

    if args.star_type and args.large_star:
        parser.error("--star-type cannot be combined with +large_star.")

    _validate_star_type(args, parser)
    _validate_name(args, parser)

    if args.intelligent_life is not None and args.habitable_world is False:
        parser.error("+intelligent_life/-intelligent_life cannot be combined with -habitable_world.")

    if args.system_file:
        # Read once here so a missing/unreadable/non-JSON/non-object file
        # is a usage error, not a traceback from build_system_config.
        try:
            file_data = run_system.load_system_file(args.system_file)
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
                and 0 <= file_num_orbits <= limits.MAX_NUM_ORBITS):
            parser.error(f"--system-file {args.system_file!r}: num_orbits {file_num_orbits!r} must be a "
                         f"whole number from 0 to {limits.MAX_NUM_ORBITS}.")

    if args.num_orbits is not None and args.num_orbits < 0:
        parser.error("--num-orbits must be zero or a positive integer.")
    if args.num_orbits is not None and args.num_orbits > limits.MAX_NUM_ORBITS:
        parser.error(f"--num-orbits must be at most {limits.MAX_NUM_ORBITS}.")

    if args.num_orbits is not None and args.planets is False:
        parser.error("--num-orbits cannot be combined with -planets.")

    _validate_star_count(args, parser, file_data if args.system_file else {})

    if args.flavor_chance_system is not None and not (0.0 <= args.flavor_chance_system <= 1.0):
        parser.error("--flavor-chance-system must be a float between 0.0 and 1.0.")

    if args.flavor_chance_planet is not None and not (0.0 <= args.flavor_chance_planet <= 1.0):
        parser.error("--flavor-chance-planet must be a float between 0.0 and 1.0.")


def _validate_star_count(args, parser, file_data):
    """
    `--stars` and `--max-stars` (GEN.267) against each other and the pair options, as the run will resolve them: a
    command-line value wins over the system file's. The message is `run_system.stars_problem`'s.
    """
    def resolved(name):
        value = getattr(args, name, None)
        return file_data.get(name) if value is None else value

    problem = run_system.max_stars_problem(resolved("max_stars")) or run_system.stars_problem(
        resolved("stars"), resolved("max_stars"), resolved("binary_system"), resolved("wide_binary"))
    if problem:
        parser.error(problem)


def add_shared_generation_options(parser):
    """
    Adds every per-system generation-tuning option the `sector` and
    `galaxy` subcommands share -- everything `iter_sector_configs()`
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
    # The forcing options are for a single system only (GEN.51). They stay
    # registered, hidden, so one given here gets a clear error, and their
    # attributes stay None ("let the generator decide") for
    # build_system_config.
    for name, _attr, _description in run_system.TRISTATE_OPTIONS:
        parser.add_argument(f'-{name}', f'+{name}', dest=name, action=SingleSystemOnlyAction,
                            nargs=0, default=None, help=argparse.SUPPRESS)

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
    parser.add_argument('--directive', action='append', default=None, metavar='MINIMUM',
                        help="Redraw each sector until it holds at least this much: systems>=N, "
                             "habitable>=N or type:X>=N (X: O B A F G K M wd bh ns). Repeat for more. "
                             "After --directive-attempts draws a habitable>=N minimum still missed is "
                             "met by a forced fill (a habitable world forced into N systems, reported as "
                             "met_forced); any other minimum the closest draw misses is reported. A minimum no "
                             "sector can reach is refused before drawing.")
    parser.add_argument('--directive-attempts', type=int, default=directives.DEFAULT_ATTEMPTS,
                        help="Draws to make at most per sector for --directive (default %(default)s).")
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
    parser.add_argument('--prevalence', type=prevalence_setting, action='append', default=None,
                        metavar='FEATURE=PERCENT',
                        help="How much more or less often every system gets FEATURE than by chance, as a "
                             "percentage: comets=+50 is 1.5 times as often, comets=-100 never. Repeat "
                             "for more features: " + ", ".join(prevalence.FEATURES) + ".")
    parser.add_argument('--stars', type=int, default=None, action=SingleSystemOnlyAction, help=argparse.SUPPRESS)
    parser.add_argument('--max-stars', type=int, default=None, metavar='N', help=MAX_STARS_HELP)
    parser.add_argument('--star-mix', type=finite_float, nargs=3, default=None,
                        metavar=('SINGLE', 'CLOSE', 'WIDE'),
                        help="The share of systems, in percent, with one star, a close binary and a wide "
                             "pair. The three must total 100. Replaces --prevalence for binary_system and "
                             "wide_binary, which it cannot be combined with.")
    parser.add_argument('--workers', type=int, default=None,
                        help="How many sectors to generate at once, each in its own low-priority worker "
                             "process. Default: 80%% of this machine's cores (one fewer when MySQL runs "
                             "here too), or PLANETGEN_WORKERS; 1 generates one sector at a time in this "
                             "process.")
    parser.add_argument('--strict', action='store_true', help=STRICT_HELP)

    add_logging_arguments(parser)


def _impossible_directive(directive, args):
    """
    Why the `--directive` minimums cannot be met by a sector of this run (GEN.180), or `None`: judged before any drawing
    from the sector's count, which is `--num-systems` exactly or, for a `sector` run, about `--density` times the
    sector's expected systems. A `galaxy` or `sectors` run draws a different count in every sector, so the densest sector
    the model has (`program_constants.PEAK_RELATIVE_DENSITY`) is judged: a minimum even that one cannot hold is refused.
    """
    if not directive:
        return None
    if args.num_systems is not None:
        return directives.impossible(directive, systems=args.num_systems)
    command = getattr(args, 'command', None)
    if command in ('sector', 'galaxy', 'sectors'):
        from planetgen.galaxy.sector import SpaceSector
        mean = (args.density if args.density is not None else 1.0) * SpaceSector("directive check").expected_system_count()
        if command != 'sector':
            mean *= program_constants.PEAK_RELATIVE_DENSITY        # a galaxy sector's count follows its place; judge the densest
        elif args.density is None:
            return None
        return directives.impossible(directive, mean=mean)
    return None


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
    problem = run_system.max_stars_problem(getattr(args, "max_stars", None))
    if problem:
        parser.error(problem)
    star_mix = getattr(args, "star_mix", None)
    if star_mix is not None:
        given = {feature for feature, _ in (args.prevalence or ())}
        if given & {"binary_system", "wide_binary"}:
            parser.error("--star-mix sets binary_system and wide_binary; don't also give them to --prevalence")
        try:
            mix = prevalence.star_mix_prevalences(*star_mix)
        except ValueError as exc:
            parser.error(f"--star-mix: {exc}")
        args.prevalence = [*(args.prevalence or ()), *mix.items()]
    if args.density is not None and args.num_systems is not None:
        parser.error("--density cannot be combined with --num-systems.")
    if args.workers is not None and args.workers < 0:
        parser.error("--workers must be 0 (automatic) or more.")
    if args.directive_attempts < 1:
        parser.error("--directive-attempts must be 1 or more.")
    try:
        directive = directives.parse(args.directive)
    except directives.DirectiveError as exc:
        parser.error(f"--directive {exc}")
    reason = _impossible_directive(directive, args)
    if reason:
        parser.error(f"--directive {reason}")

    if args.density is not None and args.density <= 0:
        parser.error("--density must be a positive number.")

    if args.density is None and args.num_systems is None:
        # `galaxy` and `sectors` modes leave both None rather than defaulting to a flat
        # count here -- run_shell_batch/run_local_neighborhood compute each
        # sector's own --density from the galaxy skeleton's real
        # position-based relative_density instead (see their own
        # docstrings and `_resolve_batch_density`), the same mechanism
        # a single-address run uses, applied uniformly across a whole batch run.
        # `getattr` (not `args.command` directly): `_default_generation_args`
        # calls this same validator against a bare, subcommand-less
        # namespace with no `command` attribute at all -- falling through
        # to the flat default there is harmless, since
        # `regenerate_sector_cell` immediately overwrites both fields
        # right after with its own already-computed density anyway.
        if getattr(args, 'command', None) not in ('galaxy', 'sectors'):
            args.num_systems = 10

    if args.num_systems is not None:
        if args.num_systems < 1:
            parser.error("--num-systems must be a positive integer.")

        if args.min_habitable > args.num_systems:
            parser.error("--min-habitable cannot exceed --num-systems.")

    if args.min_habitable < 0:
        parser.error("--min-habitable cannot be negative.")

    _validate_star_type(args, parser)
    _validate_name(args, parser)

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
    parser.add_argument('--yes', action='store_true',
                        help="With --num-sectors > 1: don't ask before generating (the size and time "
                             "question a terminal gets). A run the database disk can't hold is still refused.")
    parser.add_argument('--estimate-only', action='store_true',
                        help="Only show the size and time estimate, then stop without writing anything; "
                             "ends with one 'ESTIMATE {json}' line.")
    store.add_mysql_connection_args(parser)


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


def add_sectors_arguments(parser):
    """
    Adds the options of the commands that generate sectors in one of the modes (`sectors`, and `galaxy`, which adds the
    stages after them) -- on top of whatever `add_shared_generation_options` already added -- to `parser`: the
    `--shell`/`--block`/`--center-sector` mode-selection group, `--band`, `--wedge`, the `--shells`/`--bands`/`--wedges`
    span, the course, `--limit`, `--yes`, `--radius-pc`, `--max-shell`, `--min-start-density`, the random starts and the
    MySQL connection args.

    Args:
        parser (argparse.ArgumentParser): The parser to add options to.
    """
    mode_group = parser.add_mutually_exclusive_group(required=False)
    mode_group.add_argument('--shell', type=int, metavar='I',
                            help="Batch mode: generate every not-yet-generated sector in shell I (0-indexed "
                                 "distance from the galactic center, one sector edge deep; shell 0 is the "
                                 "core ball) that the galaxy's outline holds. Narrow it with --band (one "
                                 "band of the shell) and --wedge (one sector).")
    mode_group.add_argument('--block', metavar='M.I.S.SLAB',
                            help="Block mode: generate every not-yet-generated sector the galaxy's outline "
                                 "allows inside one block of the Galaxy Map's ladder, given by its key "
                                 "(size.shell.band.wedge, as in the map's stage links; the size is 243, "
                                 "81, 27, 9 or 3). Sectors are generated shell by shell outward, with "
                                 "the same progress and cancel as the other modes.")
    mode_group.add_argument('--around', type=finite_float, nargs=3, metavar=('X', 'Y', 'Z'),
                            help="Region mode: generate every not-yet-generated sector a sphere of --radius-pc touches, "
                                 "centred on this point (parsecs from the galactic centre), nearest first, in units of at most "
                                 "--unit-sectors consecutive wedges of one shell and band. The same seed gives the same sectors "
                                 "whatever the point, the radius or the order regions are run in. Rerun it after a stop: charted "
                                 "sectors are skipped.")
    mode_group.add_argument('--center-sector', type=int, metavar='SECTOR_ID',
                            help="Local-neighborhood mode: generate every not-yet-generated sector "
                                 "within --radius-pc of the given, already galaxy-placed sector's own "
                                 "stored center.")

    parser.add_argument('--region', action='store_true',
                        help="With --center-sector or --shell --band --wedge: move the center to the anchor sector of "
                             "the region (9-block) that holds it and fill the standard sphere of that anchor "
                             f"({program_constants.DEFAULT_GENERATE_RADIUS_PC:g} pc, about 109 sectors), so the finished fill charts and names the whole "
                             "region. Replaces --radius-pc.")
    parser.add_argument('--band', type=int, metavar='J',
                        help="With --shell: only this band of the shell (signed; 0 is the equatorial band, "
                             "positive is north, and the two end bands are the polar caps).")
    parser.add_argument('--wedge', type=int, metavar='K',
                        help="Requires --shell I and --band J: single-address mode, generating exactly the "
                             "one not-yet-generated sector at (shell I, band J, wedge K) -- e.g. the address "
                             "behind a designation copied from the interactive Galaxy Map. Skips it "
                             "(exit status 1) if it is outside the galaxy's stored outline, or reports its "
                             "existing sector_id if it was already generated. A sparse address still "
                             "generates (it may hold no systems).")
    course_group = parser.add_argument_group("course mode")
    course_group.add_argument('--course', nargs=2, metavar=('FROM', 'TO'),
                              help="Course mode: generate the uncharted sectors that block the NAV course "
                                   "between two objects (a star's object ID, or <kind>:<ID> for a body or "
                                   "phenomenon): the cells of the route's hops through unknown "
                                   "space, inside the galaxy's outline. Past NAV_CHART_CONFIRM_SECTORS "
                                   "sectors it needs --yes.")
    course_group.add_argument('--course-border', action='store_true',
                              help="With --course: also generate the uncharted cells that share a face "
                                   "with those on the course.")
    span_group = parser.add_argument_group("span mode")
    span_group.add_argument('--shells', metavar='FIRST:LAST',
                            help="Span mode: generate every not-yet-generated sector in these shells (N or "
                                 "FIRST:LAST, inclusive). Default: every shell the outline reaches.")
    span_group.add_argument('--bands', metavar='FIRST:LAST',
                            help="Span mode: only these bands (N or FIRST:LAST, inclusive; negative south of "
                                 "the equator, so write --bands=-2:2). Default: every band in the outline. "
                                 "Alone it fills those bands out to the galaxy's edge.")
    span_group.add_argument('--wedges', metavar='FIRST:LAST',
                            help="Span mode, with exactly one shell and one band: only this arc of wedges, "
                                 "inclusive; FIRST past LAST wraps through wedge 0 (50:5). Default: the "
                                 "whole band.")
    parser.add_argument('--limit', type=int,
                        help="With --shell (or --block, or a span): generate only the first N "
                             "not-yet-generated sectors.")
    parser.add_argument('--yes', action='store_true',
                        help="Don't ask before generating: skips the size and time question asked on a "
                             "terminal before filling more than one sector, and with --shell (or --block, "
                             "or a span) the confirmation normally required before generating "
                             f"more than {_RING_LIMIT} sectors. A run the database disk "
                             "can't hold is still refused.")
    parser.add_argument('--estimate-only', action='store_true',
                        help="Only show the size and time estimate (and whether the disk would refuse the "
                             "run), then stop without writing anything; ends with one 'ESTIMATE {json}' "
                             "line. In random-start mode the estimate is for one random start.")
    parser.add_argument('--radius-pc', type=finite_float,
                        help="With --center-sector: the neighborhood search radius, in parsecs. With "
                             "--shell --band --wedge: after generating that address, also generate its "
                             "neighborhood within this radius. With "
                             "neither --shell nor --center-sector (random-start mode): overrides the "
                             "default 12 pc neighborhood radius around the randomly chosen starting "
                             "sector.")
    parser.add_argument('--grow-shells', type=int, metavar='K',
                        help="With --around: grow the sphere by K whole shells of one sector edge (about the surface times K "
                             "more sectors).")
    parser.add_argument('--unit-sectors', type=int, metavar='N',
                        help="With --around: the most sectors in a unit of work, the retry and resume grain (default "
                             f"{region_run.BLOCK_SECTORS}).")
    parser.add_argument('--max-shell', type=int,
                        help="With neither --shell nor --center-sector (random-start mode): the highest "
                             "shell the randomly chosen starting sector may land in. Default: anywhere "
                             "inside the galaxy's stored outline.")
    parser.add_argument('--min-start-density', type=finite_float,
                        help="With neither --shell nor --center-sector (random-start mode): require the "
                             "randomly chosen starting sector's own real relative_density (the same "
                             "'expected' figure printed alongside each saved sector) to be at least this "
                             "value before accepting it -- e.g. 1.0 for at least as dense as the galaxy's "
                             "own real local density. Retried the same way an already-occupied "
                             "address is (see RANDOM_START_MAX_PLACEMENT_ATTEMPTS). Cannot "
                             "be combined with --density/--num-systems (those override every position's "
                             "density uniformly, leaving no per-position value to compare against).")
    parser.add_argument('--neighborhoods', type=int, default=1, metavar='N',
                        help="Random-start mode: generate N random neighborhoods instead of one. Each "
                             "start is inside the galaxy with its whole neighborhood, and at least twice the "
                             "radius from every other start, so they never overlap.")
    parser.add_argument('--neighborhood-gamma', type=finite_float, default=0.0, metavar='G',
                        help="With --neighborhoods: bias the starts toward dense space. A start is kept with "
                             "probability min(1, relative density) ** G; 0 (the default) keeps every one, "
                             "so starts are uniform by volume.")
    parser.add_argument('--avoid-filled-space', action='store_true',
                        help="Random-start mode: drop any random centre whose neighborhood touches a "
                             "sector that is already generated, so every neighborhood grows in empty space. "
                             "Logs how many candidate centres qualified. Without it, a neighborhood may "
                             "overlap filled space (the sectors already there are kept, not redone).")
    store.add_mysql_connection_args(parser)


def add_galaxy_arguments(parser):
    """
    Adds the `galaxy` subcommand's own options -- on top of whatever `add_shared_generation_options` already added
    -- to `parser`: the sector options (`add_sectors_arguments`) and the options of the scatters after the sectors.

    Args:
        parser (argparse.ArgumentParser): The parser to add options to.
    """
    add_sectors_arguments(parser)
    scatter_group = parser.add_argument_group("the scatters after the sectors")
    scatter_group.add_argument('--bright-star-min-luminosity', type=finite_float,
                               default=program_constants.BRIGHT_STAR_MIN_LUMINOSITY_SOL,
                               help="The scatter's threshold, solar luminosities. "
                                    f"Default: {program_constants.BRIGHT_STAR_MIN_LUMINOSITY_SOL:g}.")
    scatter_group.add_argument('--phenomenon-min-mass', type=finite_float, default=None, metavar='M_SUN',
                               help="The stellar mass limit of the scatters, solar "
                                    "masses; blank keeps the one already stored, else "
                                    f"{program_constants.PHENOMENON_MIN_MASS_SOLAR:g}.")
    scatter_group.add_argument('--compact-min-mass', type=compact_mass_limit, default=None, metavar='M_SUN',
                               help=COMPACT_MIN_MASS_HELP)


def validate_galaxy_args(args, parser):
    """
    Validates `add_galaxy_arguments`'s own options, calling `parser.error`
    (which exits) on the first problem found, then shapes `args` into the
    namespace `generate_sector()`/`iter_sector_configs()` expect.
    Doesn't include `add_shared_generation_options`'s own validation
    (`validate_shared_generation_args`) -- callers run both.

    Args:
        args (argparse.Namespace): Parsed arguments.
        parser (argparse.ArgumentParser): The parser to raise errors
                                          through (so the caller's own
                                          `--help`/usage text is shown).
    """
    check_compact_mass_limit(args, parser)
    try:
        args.bright_star_min_luminosity = luminosity_floor.check(args.bright_star_min_luminosity)
    except ValueError as exc:
        parser.error(f"--bright-star-min-luminosity: {exc}")
    validate_sectors_args(args, parser)


def validate_sectors_args(args, parser):
    """
    Validates `add_sectors_arguments`'s options (the modes and what goes with each), calling `parser.error` on the first
    problem found, then shapes `args` into the namespace `generate_sector()`/`iter_sector_configs()` expect. Doesn't
    include `add_shared_generation_options`'s own validation (`validate_shared_generation_args`).
    """
    if getattr(args, "around", None) is not None:
        _validate_around(args, parser)
        return
    for flag, value in (("--grow-shells", getattr(args, "grow_shells", None)), ("--unit-sectors", getattr(args, "unit_sectors", None))):
        if value is not None:
            parser.error(f"{flag} requires --around.")

    if getattr(args, "block", None) is not None:
        _validate_block(args, parser)
        return

    course = getattr(args, "course", None)
    if course is not None:
        _validate_course(args, parser)
        return
    if getattr(args, "course_border", False):
        parser.error("--course-border requires --course.")

    if any(value is not None for value in (args.shells, args.bands, args.wedges)):
        _validate_span(args, parser)
        return
    random_start = args.shell is None and args.center_sector is None

    if args.shell is not None and args.shell < 0:
        parser.error("--shell must be >= 0.")
    if args.shell is not None and args.shell > limits.MAX_GENERATE_SHELL:
        parser.error(f"--shell must be at most {limits.MAX_GENERATE_SHELL}.")
    if args.band is not None and args.shell is None:
        parser.error("--band requires --shell.")
    if args.wedge is not None and (args.shell is None or args.band is None):
        parser.error("--wedge requires --shell and --band.")
    if args.wedge is not None and args.wedge < 0:
        parser.error("--wedge must be >= 0.")
    if args.band is not None:
        try:
            low, high = spherical.band_range(args.shell)
        except ValueError as exc:
            parser.error(f"--shell: {exc}")
        if not low <= args.band <= high:
            parser.error(f"--band: shell {args.shell} holds bands {low} to {high}.")
        if args.wedge is not None and args.wedge >= spherical.wedge_count(args.shell, args.band):
            parser.error(f"--wedge: shell {args.shell} band {args.band} holds wedges "
                         f"0 to {spherical.wedge_count(args.shell, args.band) - 1}.")
    if args.wedge is not None and args.limit is not None:
        parser.error("--limit only applies to --shell (batch mode), not --shell --band --wedge (single-address mode).")
    if args.wedge is not None and args.yes:
        parser.error("--yes only applies to --shell (batch mode), not --shell --band --wedge (single-address mode).")
    if args.wedge is not None and (args.density is not None or args.num_systems is not None):
        parser.error("--density/--num-systems can't be combined with --shell --band --wedge (single-address "
                     "mode) -- it always uses this address's own real predicted "
                     "density, the same as when the galaxy map's own live view found it.")

    if args.region and args.center_sector is None and args.wedge is None:
        parser.error("--region requires --center-sector or --shell --band --wedge.")
    if args.region and args.radius_pc is not None:
        parser.error(f"--region fills the standard {program_constants.DEFAULT_GENERATE_RADIUS_PC:g} pc sphere of the anchor; it can't be combined with "
                     "--radius-pc.")
    if args.center_sector is not None and args.radius_pc is None and not args.region:
        parser.error("--center-sector requires --radius-pc (or --region).")
    if args.radius_pc is not None and args.shell is not None and args.wedge is None:
        parser.error("--radius-pc only applies to --center-sector, --shell --band --wedge, or random-start mode "
                     "(neither --shell nor --center-sector), not --shell alone.")
    if args.radius_pc is not None and args.radius_pc <= 0:
        parser.error("--radius-pc must be a positive number.")
    if args.radius_pc is not None and args.radius_pc > limits.MAX_GENERATE_RADIUS_PC:
        parser.error(f"--radius-pc must be at most {limits.MAX_GENERATE_RADIUS_PC:g}.")

    if args.limit is not None and args.shell is None:
        parser.error("--limit only applies to --shell.")
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be a positive integer.")
    if args.limit is not None and args.limit > limits.MAX_GENERATE_LIMIT:
        parser.error(f"--limit must be at most {limits.MAX_GENERATE_LIMIT}.")
    if args.yes and args.shell is None:
        parser.error("--yes only applies to --shell.")

    if args.neighborhoods < 1 or args.neighborhoods > limits.MAX_GENERATE_NEIGHBORHOODS:
        parser.error(f"--neighborhoods must be between 1 and {limits.MAX_GENERATE_NEIGHBORHOODS}.")
    if args.neighborhood_gamma < 0:
        parser.error("--neighborhood-gamma must be 0 or more.")
    if (args.neighborhoods > 1 or args.neighborhood_gamma) and not random_start:
        parser.error("--neighborhoods and --neighborhood-gamma only apply to random-start mode (neither "
                     "--shell nor --center-sector).")
    if args.avoid_filled_space and not random_start:
        parser.error("--avoid-filled-space only applies to random-start mode (neither --shell nor --center-sector).")
    if args.max_shell is not None and not random_start:
        parser.error("--max-shell only applies to random-start mode (neither --shell nor --center-sector).")
    if args.max_shell is not None and args.max_shell < 0:
        parser.error("--max-shell must be >= 0.")
    if args.max_shell is not None and args.max_shell > limits.MAX_GENERATE_SHELL:
        parser.error(f"--max-shell must be at most {limits.MAX_GENERATE_SHELL}.")

    if args.min_start_density is not None and not random_start:
        parser.error("--min-start-density only applies to random-start mode (neither --shell nor "
                     "--center-sector).")
    if args.min_start_density is not None and args.min_start_density <= 0:
        parser.error("--min-start-density must be a positive number.")
    if args.min_start_density is not None and (args.density is not None or args.num_systems is not None):
        parser.error("--min-start-density cannot be combined with --density/--num-systems -- those "
                     "override every position's density uniformly, leaving no per-position value for "
                     "--min-start-density to compare against.")

    # generate_sector()/iter_sector_configs() expect a namespace shaped
    # like the `sector` subcommand's own output -- see this function's
    # docstring for why these are always None here.
    args.sector_name = None
    args.system_file = None
    args.num_orbits = None
    args.name = None


def _validate_course(args, parser):
    """Refuses what can't go with `--course` (NAV.48) and shapes `args` like the other modes."""
    others = [flag for flag, value in (
        ("--shell", args.shell), ("--band", args.band), ("--wedge", args.wedge),
        ("--center-sector", args.center_sector), ("--radius-pc", args.radius_pc),
        ("--max-shell", args.max_shell), ("--min-start-density", args.min_start_density),
        ("--region", args.region or None),
        ("--shells", args.shells), ("--bands", args.bands), ("--wedges", args.wedges), ("--limit", args.limit),
    ) if value is not None]
    if others:
        parser.error(f"--course can't be combined with {', '.join(others)}.")
    args.sector_name = args.system_file = args.num_orbits = args.name = None


def _validate_around(args, parser):
    """Checks `--around` (GEN.1128): a start point, a positive `--radius-pc`, whole non-negative `--grow-shells` and a positive `--unit-sectors`,
    nothing that names another region; shapes `args` like the other modes."""
    others = [flag for flag, value in (
        ("--shell", args.shell), ("--band", args.band), ("--wedge", args.wedge), ("--center-sector", args.center_sector),
        ("--max-shell", args.max_shell), ("--min-start-density", args.min_start_density), ("--limit", args.limit),
        ("--shells", args.shells), ("--bands", args.bands), ("--wedges", args.wedges), ("--region", args.region or None),
    ) if value is not None]
    if others:
        parser.error(f"--around can't be combined with {', '.join(others)}.")
    if args.radius_pc is None:
        parser.error("--around requires --radius-pc.")
    if args.radius_pc <= 0:
        parser.error("--radius-pc must be a positive number.")
    if args.radius_pc > limits.MAX_GENERATE_RADIUS_PC:
        parser.error(f"--radius-pc must be at most {limits.MAX_GENERATE_RADIUS_PC:g}.")
    if args.grow_shells is not None and args.grow_shells < 0:
        parser.error("--grow-shells must be >= 0.")
    if args.unit_sectors is not None and args.unit_sectors < 1:
        parser.error("--unit-sectors must be >= 1.")
    args.sector_name = args.system_file = args.num_orbits = args.name = None


def _validate_block(args, parser):
    """Checks `--block` (GEN.310): a key of a block of a stage level (`drill.parse_stage_key`), nothing that names another
    region, a `--limit` in range; shapes `args` like the other modes."""
    others = [flag for flag, value in (
        ("--band", args.band), ("--wedge", args.wedge), ("--radius-pc", args.radius_pc),
        ("--max-shell", args.max_shell), ("--min-start-density", args.min_start_density),
        ("--shells", args.shells), ("--bands", args.bands), ("--wedges", args.wedges),
        ("--region", args.region or None),
    ) if value is not None]
    if others:
        parser.error(f"--block can't be combined with {', '.join(others)}.")
    try:
        block = drill.parse_stage_key(args.block)
    except ValueError as exc:
        parser.error(f"--block: {exc}")
    args.block = drill.format_drill_key(block)
    if args.limit is not None and not 1 <= args.limit <= limits.MAX_GENERATE_LIMIT:
        parser.error(f"--limit must be between 1 and {limits.MAX_GENERATE_LIMIT}.")
    args.sector_name = args.system_file = args.num_orbits = args.name = None


def _validate_span(args, parser):
    """Turns `--shells/--bands/--wedges` into `args.span` (a `Span`) and
    refuses what can't go with a span (ADM.29)."""
    others = [flag for flag, value in (
        ("--shell", args.shell), ("--band", args.band), ("--wedge", args.wedge),
        ("--center-sector", args.center_sector), ("--radius-pc", args.radius_pc),
        ("--max-shell", args.max_shell), ("--min-start-density", args.min_start_density),
        ("--region", args.region or None),
    ) if value is not None]
    if others:
        parser.error(f"--shells/--bands/--wedges can't be combined with {', '.join(others)}.")
    if args.wedges is not None and (args.shells is None or args.bands is None):
        parser.error("--wedges requires --shells and --bands (one shell, one band).")
    try:
        span = Span(
            shells=parse_range(args.shells, "--shells") if args.shells is not None else None,
            bands=parse_range(args.bands, "--bands") if args.bands is not None else None,
            wedges=parse_range(args.wedges, "--wedges") if args.wedges is not None else None)
    except SpanError as exc:
        parser.error(str(exc))
    if span.shells is not None and span.shells[1] > limits.MAX_GENERATE_SHELL:
        parser.error(f"--shells must stay at or below shell {limits.MAX_GENERATE_SHELL}.")
    if args.limit is not None and not 1 <= args.limit <= limits.MAX_GENERATE_LIMIT:
        parser.error(f"--limit must be between 1 and {limits.MAX_GENERATE_LIMIT}.")
    args.span = span
    args.sector_name = args.system_file = args.num_orbits = args.name = None


def _galaxy_seed_arg(text):
    """`--seed`'s type: 32 hex digits, as the galaxy's 16-byte seed."""
    try:
        return galaxySeed.parse_seed(text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from None


def add_plan_arguments(parser):
    """
    Adds every option the `plan` subcommand accepts (besides `--version`)
    to `parser` -- the galaxy shape parameter group plus the scan options
    and the MySQL connection args.

    Args:
        parser (argparse.ArgumentParser): The parser to add options to.
    """
    shape_group = parser.add_argument_group(
        "galaxy shape (the lever table, planetgen.galaxy.dials: every option is Auto, drawn from the galaxy seed, "
        "until given)")
    dials.add_arguments(shape_group)

    parser.add_argument('--seed', type=_galaxy_seed_arg, default=None, metavar="HEX",
                        help="The galaxy's 128-bit seed, as 32 hex digits: the same seed makes the same galaxy "
                             "on the same PlanetGen release. Default: the seed already stored, or one drawn at "
                             "random for a new galaxy. Changing it needs a galaxy with no sectors.")
    parser.add_argument('--reach-shells', type=int, default=DEFAULT_REACH_SHELL, metavar="N",
                        help="The outermost shell the grid has addresses for (the galaxy's reach; a shell is 4 pc "
                             f"deep). Default: {DEFAULT_REACH_SHELL} (262 kpc). The object ID can name shells up to "
                             f"{spherical.MAX_ADDRESSABLE_SHELL}; the galaxy's disk outline must fit inside it.")

    parser.add_argument('--workers', type=int, default=None,
                        help="How many worker processes the run uses at most. Default: 80%% of this machine's cores "
                             "(one fewer when MySQL runs here too), or PLANETGEN_WORKERS.")
    resolution_group = parser.add_argument_group("how deep a dense cell is filled (planetgen.galaxy.resolution)")
    resolution_group.add_argument('--cell-star-ceiling', type=finite_float, default=None, metavar='SYSTEMS',
                                  help="A cell that expects more systems than this is filled only with the stars at "
                                       "or above --fill-luminosity-cap; the rest is unresolved, its count and mass "
                                       "shown on the sector page. Stored with the galaxy by a full plan. Default: "
                                       f"{program_constants.CELL_STAR_CEILING:,.0f}.")
    resolution_group.add_argument('--fill-luminosity-cap', type=finite_float, default=None, metavar='L_SUN',
                                  help="The dimmest star (solar luminosities) a cell above --cell-star-ceiling "
                                       "resolves, at least the brightest white dwarf "
                                       f"({program_constants.WD_LUMINOSITY_RANGE_SOL[1]:g}). Stored with the galaxy "
                                       f"by a full plan. Default: {program_constants.FILL_LUMINOSITY_CAP_SOL:g}.")
    limits_group = parser.add_argument_group("scatter limits by place (planetgen.galaxy.scatter_limits)")
    limits_group.add_argument('--core-min-luminosity', type=finite_float, default=None, metavar='L_SUN',
                              help="The luminosity floor (solar luminosities) of the scatter in the galactic core. "
                                   "Stored with the galaxy by a full plan. Default: the baseline floor times "
                                   f"{program_constants.CORE_TIER_FACTOR:g}, until set.")
    limits_group.add_argument('--core-mass-limit', type=finite_float, default=None, metavar='M_SUN',
                              help="The mass limit (solar masses) of the scatter in the galactic core. Stored with "
                                   "the galaxy by a full plan. Default: the baseline limit times "
                                   f"{program_constants.CORE_TIER_FACTOR:g}, until set.")
    limits_group.add_argument('--halo-min-luminosity', type=finite_float, default=None, metavar='L_SUN',
                              help="The luminosity floor (solar luminosities) of the scatter in the halo, from "
                                   f"{program_constants.HALO_MIN_LUMINOSITY_RANGE_SOL[0]:g} to "
                                   f"{program_constants.HALO_MIN_LUMINOSITY_RANGE_SOL[1]:g}; the halo has no mass "
                                   "limit of its own. Stored with the galaxy by a full plan. Default: "
                                   f"{program_constants.HALO_MIN_LUMINOSITY_SOL:g}.")
    parser.add_argument('--timing', action='store_true',
                        help="Print a table of where the plan's time went, step by step, the costliest first "
                             "(also printed with --debug).")
    store.add_mysql_connection_args(parser)
    add_logging_arguments(parser)


def validate_plan_args(args, parser):
    """
    Validates `add_plan_arguments`'s own options, calling `parser.error`
    (which exits) on the first problem found.

    Args:
        args (argparse.Namespace): Parsed arguments.
        parser (argparse.ArgumentParser): The parser to raise errors through.
    """
    if args.workers is not None and args.workers < 0:
        parser.error("--workers must be 0 (automatic) or more.")
    if args.cell_star_ceiling is not None and args.cell_star_ceiling < 1:
        parser.error("--cell-star-ceiling must be at least 1 system.")
    if args.fill_luminosity_cap is not None:
        try:
            star_population.check_threshold(args.fill_luminosity_cap)
        except ValueError as exc:
            parser.error(f"--fill-luminosity-cap: {exc}")
    for option in ("core_min_luminosity", "core_mass_limit"):
        if getattr(args, option) is not None and getattr(args, option) <= 0:
            parser.error(f"--{option.replace('_', '-')} must be more than 0.")
    if args.halo_min_luminosity is not None:
        try:
            scatter_limits.check_halo_min_luminosity(args.halo_min_luminosity)
        except ValueError as exc:
            parser.error(f"--halo-min-luminosity: {exc}")
    # Anything else that still can't be normalized.
    try:
        shape, _drawn = shapeDraw.plan_shape(args, args.seed)
    except shapeDraw.ScaleError as exc:
        parser.error(f"mass, size and density: {exc}.")
    except (ArithmeticError, ValueError) as exc:
        parser.error(f"these galaxy shape parameters can't be normalized ({exc}).")
    if not (math.isfinite(shape.k_norm) and shape.k_norm > 0):
        parser.error(f"these galaxy shape parameters can't be normalized (k_norm={shape.k_norm!r}).")


SCATTER_KIND_HELP = {
    "mass": "Scatter the massive stars: every star born at or above --phenomenon-min-mass, galaxy-wide. Refused "
            "while the massive stars are stored: run 'planetgen clear-scatter --luminosity' and "
            "'--mass' first.",
    "luminosity": "Scatter the bright stars: every star born lighter than the stored mass limit and at least "
                  "--bright-star-min-luminosity bright (needs --mass run first). Refused while the bright stars are "
                  "stored: run 'planetgen clear-scatter --luminosity' first.",
    "phenomena": "Scatter the phenomena (black holes, neutron stars, planetary nebulae, supernova remnants, "
                 "hypervelocity stars, the nucleus); neutron stars and black holes from --compact-min-mass up. "
                 "Refused while they are stored: run 'planetgen clear-scatter --phenomena' first.",
    "nebula-field": "Draw the nebula field table (the galaxy's molecular cloud field of classes M and N) from the "
                    "galaxy seed. Refused while it is stored: run 'planetgen clear-scatter --nebula-field' first.",
    "band": "Go one level dimmer on the stored bright stars: keep the stars already placed and add only those from "
            "--down-to up to the level already scattered.",
}
"""dict: The help of each kind of `planetgen scatter`."""


def add_scatter_arguments(parser):
    """
    Adds the options of the `scatter` subcommand (PERF.121): exactly one kind of scatter and only that kind's
    parameters, then the run options and the MySQL connection.
    """
    kind = parser.add_argument_group("the scatter to draw (exactly one)")
    chosen = kind.add_mutually_exclusive_group(required=True)
    for name, text in SCATTER_KIND_HELP.items():
        chosen.add_argument(f"--{name}", dest=name.replace("-", "_"), action="store_true", help=text)
    parser.add_argument('--bright-star-min-luminosity', type=finite_float,
                        default=program_constants.BRIGHT_STAR_MIN_LUMINOSITY_SOL,
                        help="For --luminosity (and --mass, which records it): every star at least this bright "
                             "(solar luminosities) is generated and placed galaxy-wide: from "
                             f"{program_constants.BRIGHT_STAR_FLOOR_MIN_SOL:g} to "
                             f"{program_constants.BRIGHT_STAR_FLOOR_MAX_SOL:g} (the Generate page offers presets). "
                             f"Default: {program_constants.BRIGHT_STAR_MIN_LUMINOSITY_SOL:g}.")
    presets = ", ".join(f"{value:g}" for value in program_constants.PHENOMENON_MIN_MASS_PRESETS)
    parser.add_argument('--phenomenon-min-mass', type=finite_float, default=None, metavar='M_SUN',
                        help="For --mass: the stellar mass limit (solar masses) from which every star is placed "
                             "galaxy-wide; a sector draws the lighter ones when it is filled. One of "
                             f"{presets}. Default: the galaxy's stored limit, else "
                             f"{program_constants.PHENOMENON_MIN_MASS_SOLAR:g}.")
    parser.add_argument('--compact-min-mass', type=compact_mass_limit, default=None, metavar='M_SUN',
                        help="For --phenomena: " + COMPACT_MIN_MASS_HELP)
    parser.add_argument('--down-to', type=finite_float, default=None, metavar='L_SUN',
                        help="For --band: the luminosity (solar luminosities) to go down to. Sectors already filled "
                             "are left out. Does nothing when L_SUN is not below the current level.")
    parser.add_argument('--workers', type=int, default=None,
                        help="How many worker processes the run uses at most. Default: 80%% of this machine's cores "
                             "(one fewer when MySQL runs here too), or PLANETGEN_WORKERS.")
    store.add_mysql_connection_args(parser)
    add_logging_arguments(parser)


def validate_scatter_args(args, parser):
    """Validates `add_scatter_arguments`'s options, calling `parser.error` on the first problem found."""
    if args.workers is not None and args.workers < 0:
        parser.error("--workers must be 0 (automatic) or more.")
    kind = "nebula_field" if args.nebula_field else next(name for name in ("mass", "luminosity", "phenomena", "band")
                                                         if getattr(args, name))
    others = {"phenomenon_min_mass": "mass", "compact_min_mass": "phenomena", "down_to": "band"}
    for option, belongs in others.items():
        if getattr(args, option) is not None and kind != belongs:
            parser.error(f"--{option.replace('_', '-')} belongs to --{belongs}, not --{kind.replace('_', '-')}.")
    if kind == "band":
        if args.down_to is None:
            parser.error("--band needs --down-to L_SUN.")
        try:
            bright_star_fraction(args.down_to)
        except ValueError as exc:
            parser.error(f"--down-to: {exc}")
    if args.phenomenon_min_mass is not None and not any(
            math.isclose(args.phenomenon_min_mass, preset) for preset in program_constants.PHENOMENON_MIN_MASS_PRESETS):
        parser.error("--phenomenon-min-mass must be one of "
                     + ", ".join(f"{value:g}" for value in program_constants.PHENOMENON_MIN_MASS_PRESETS) + ".")
    check_compact_mass_limit(args, parser)
    if kind in ("mass", "luminosity"):
        try:
            args.bright_star_min_luminosity = luminosity_floor.check(args.bright_star_min_luminosity)
            bright_star_fraction(args.bright_star_min_luminosity)
        except ValueError as exc:
            parser.error(f"--bright-star-min-luminosity: {exc}")


CLEAR_KIND_HELP = {
    "mass": "Empty the massive stars (refused while the bright stars are stored, which are split from these at the "
            "mass limit: clear them first).",
    "luminosity": "Empty the bright stars, the band with them.",
    "phenomena": "Empty the phenomena scatter.",
    "nebula-field": "Empty the nebula field table.",
}
"""dict: The help of each kind of `planetgen clear-scatter`."""


def add_clear_scatter_arguments(parser):
    """Adds the options of `clear-scatter` (Boss's F7): exactly one kind to clear, then the MySQL connection."""
    kind = parser.add_argument_group("the scatter to clear (exactly one)")
    chosen = kind.add_mutually_exclusive_group(required=True)
    for name, text in CLEAR_KIND_HELP.items():
        chosen.add_argument(f"--{name}", dest=name.replace("-", "_"), action="store_true", help=text)
    store.add_mysql_connection_args(parser)
    add_logging_arguments(parser)


def add_fill_arguments(parser):
    """Adds the options of a command that only runs units on the queue and takes no setting of its own
    (`fill-map-luminosity`, `rescore-habitability`): the worker count, the MySQL connection and the logging."""
    parser.add_argument('--workers', type=int, default=None,
                        help="How many worker processes the run uses at most. Default: 80%% of this machine's cores "
                             "(one fewer when MySQL runs here too), or PLANETGEN_WORKERS.")
    store.add_mysql_connection_args(parser)
    add_logging_arguments(parser)


def _database_time(text):
    """A `--since` value: a time on the database's clock, `YYYY-MM-DD HH:MM:SS` (a `T` between the date and the time
    also works)."""
    import datetime
    try:
        return datetime.datetime.fromisoformat(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{text!r} is not a time like '2026-10-22 17:30:00'") from None


def add_since_arguments(parser):
    """Adds the options of a command that works on the galaxy sectors saved since a time (`link`, `backfill`, `settle`,
    PERF.123): `--since` (the database's clock; a `sectors` or `galaxy` run says its start time), then the common ones."""
    parser.add_argument('--since', type=_database_time, default=None, metavar='TIME',
                        help="Work on the galaxy sectors saved at or after TIME, 'YYYY-MM-DD HH:MM:SS' on the "
                             "database's clock. Default: the start of the newest 'sectors' or 'galaxy' run.")
    add_fill_arguments(parser)


def validate_fill_args(args, parser):
    """Validates `add_fill_arguments`'s options (and that `settle` is given --since or --sector-id, not both)."""
    if args.workers is not None and args.workers < 0:
        parser.error("--workers must be 0 (automatic) or more.")
    if args.command == 'settle' and args.sector_id and args.since is not None:
        parser.error("--since and --sector-id cannot be combined.")


def add_clean_nebula_field_arguments(parser):
    """Adds the options of `clean-nebula-field` (PERF.122): `--apply`, then the common ones."""
    parser.add_argument('--apply', action='store_true',
                        help="Delete the stored field clouds the current keep rate rejects (clearing what pointed at "
                             "them and refreshing containment in the sectors they reached). Without it only counts "
                             "them (a dry run). The nebula field table is drawn again by "
                             "'planetgen scatter --nebula-field'.")
    add_fill_arguments(parser)


def add_phenomenon_arguments(parser):
    """
    Adds every option the `phenomenon` subcommand accepts (besides
    `--version`) to `parser` -- `--type`, `--anchor-system`,
    `--num-orbits`, `--sector-id`, the MySQL connection args,
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
                              "planetgen.db.store.compute_phenomenon_placement). Also links a "
                              "supernova-remnant/rogue-planet/comet to that sector, without a computed "
                              "galaxy position (those types have no placement columns of their own). A "
                              "quasar needs a shell-0 (galactic core) sector and is placed at the galactic "
                              "center; only one per galaxy. Omit to generate it unplaced/unlinked, as before.")

    # Database persistence
    store.add_mysql_connection_args(parser)

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
    if args.num_orbits is not None and args.num_orbits > limits.MAX_NUM_ORBITS:
        parser.error(f"--num-orbits must be at most {limits.MAX_NUM_ORBITS}.")
    if args.anchor_system and args.type not in (None, "black-hole", "neutron-star"):
        parser.error("--anchor-system is only valid with --type black-hole or --type neutron-star.")
    if args.sector_id is not None and args.anchor_system:
        parser.error("--sector-id cannot be combined with --anchor-system -- an anchored compact remnant "
                     "belongs to its own StarSystem, not directly to a sector.")


def add_fingerprint_arguments(parser):
    """
    Adds the `fingerprint` subcommand's options: the region (`--shell`,
    `--sector`; the whole galaxy without either), and the MySQL connection
    and logging args.

    Args:
        parser (argparse.ArgumentParser): The parser to add options to.
    """
    parser.add_argument('--shell', type=int, action='append', metavar='SHELL',
                        help="Fingerprint every sector in this shell (repeatable).")
    parser.add_argument('--sector', type=int, nargs=3, action='append', metavar=('SHELL', 'BAND', 'WEDGE'),
                        help="Fingerprint the sector at this address (repeatable).")
    parser.add_argument('--stored', action='store_true',
                        help="Print the tree of digests stored at save time instead (a line a band and a shell, the "
                             "galaxy root, and the region's digest with its sector count).")
    store.add_mysql_connection_args(parser)
    add_logging_arguments(parser)


def add_population_arguments(parser):
    """
    Adds the `population` subcommand's options: `--rescan`, and the MySQL connection and logging args.

    Args:
        parser (argparse.ArgumentParser): The parser to add options to.
    """
    parser.add_argument('--rescan', action='store_true',
                        help="Forget every species, polity and territory and scan every planet again "
                             "(new names, ages and borders).")
    store.add_mysql_connection_args(parser)
    add_logging_arguments(parser)


def add_territories_arguments(parser):
    """Adds the `territories` subcommand's options (PERF.125): the MySQL connection and logging args and nothing else."""
    store.add_mysql_connection_args(parser)
    add_logging_arguments(parser)


def add_bring_to_epoch_arguments(parser):
    """
    Adds the `bring-to-epoch` subcommand's options: `--dry-run`, and the MySQL connection and logging args.

    Args:
        parser (argparse.ArgumentParser): The parser to add options to.
    """
    parser.add_argument('--dry-run', action='store_true',
                        help="Say what would be generated again, which sectors would be left alone and why, and "
                             "change nothing.")
    store.add_mysql_connection_args(parser)
    add_logging_arguments(parser)


def run_bring_to_epoch(args):
    """`planetgen bring-to-epoch` (GEN.764): brings the stored galaxy to the running generator epoch in part
    (`planetgen.galaxy.bring_to_epoch`); exits 1 when the run refused or a sector failed."""
    config = store.mysql_config_from_args(args)
    result = bringToEpoch.run(config, store.control_mysql_config(config), dry_run=args.dry_run)
    if result["refused"] or result["failed"]:
        sys.exit(1)


def add_refresh_parity_arguments(parser):
    """
    Adds the `refresh-parity` subcommand's options: `--compact`, `--limit`, and the MySQL connection and logging args.

    Args:
        parser (argparse.ArgumentParser): The parser to add options to.
    """
    parser.add_argument('--compact', action='store_true',
                        help="Afterwards rewrite the parity file with only the latest record of each group.")
    parser.add_argument('--limit', type=int, default=None, metavar='N',
                        help="Make at most N sectors current (the rest at the next run).")
    store.add_mysql_connection_args(parser)
    add_logging_arguments(parser)


def run_refresh_parity(args):
    """`planetgen refresh-parity` (DB.46): makes the parity file of the sector storage current for the whole database
    (`planetgen.db.parity.refresh`), and with `--compact` rewrites it without its dead space. Exits 1 when a sector could
    not be protected (a damaged member of its group)."""
    config = store.mysql_config_from_args(args)
    conn = store.get_connection(config)
    try:
        counts = parity.refresh(conn, limit=args.limit)
    finally:
        conn.close()
    log.normal("Parity: " + (", ".join(f"{count:,} {result}" for result, count in sorted(counts.items())) or "nothing to do") + ".")
    if args.compact:
        kept, dropped = parity.compact()
        log.normal(f"Parity: compacted to {kept:,} group record(s)" + (f", {dropped:,} unreadable dropped." if dropped else "."))
    if counts.get(parity.UNPROTECTABLE) or counts.get("failed"):
        sys.exit(1)


def add_repair_db_arguments(parser):
    """
    Adds the `repair-db` subcommand's options: `--shell` and `--sector` (as `check-db`), `--dry-run`, and the MySQL
    connection and logging args.

    Args:
        parser (argparse.ArgumentParser): The parser to add options to.
    """
    parser.add_argument("--shell", type=int, action="append", metavar="SHELL",
                        help="Look only at this shell's sectors (repeatable).")
    parser.add_argument("--sector", type=int, nargs=3, action="append", metavar=("SHELL", "BAND", "WEDGE"),
                        help="Look only at the sector at this address (repeatable).")
    parser.add_argument('--dry-run', action='store_true',
                        help="Say which damaged sectors the parity file could rebuild and write nothing.")
    parser.add_argument('--salvage', metavar='TABLE', default=None,
                        help="Step 1 of rebuilding a table whose clustered index has an unreadable page (run with the server "
                             "started with innodb_force_recovery = 1 and the website and generation stopped): dump its readable "
                             "rows to --dump-dir, then drop the table and make it again empty.")
    parser.add_argument('--reload', metavar='TABLE', default=None,
                        help="Step 2 (after a normal restart): insert the dumped rows of TABLE back, then repair the sectors "
                             "that lost rows from the parity file.")
    parser.add_argument('--restore-sectors', action='store_true',
                        help="With --reload sectors: also put back the sector rows the parity file holds and the table lacks.")
    parser.add_argument('--dump-dir', metavar='DIR', default=None,
                        help="Where --salvage writes and --reload reads the dump (default: a 'salvage' folder in the parity folder).")
    store.add_mysql_connection_args(parser)
    add_logging_arguments(parser)


def run_repair_db(args):
    """`planetgen repair-db` (DB.47, DB.48): rebuilds the sectors that no longer hash to their stored checksum from the parity file
    (`planetgen.db.repair`), and makes the parity file again from the database if it is the damaged part. Exits 1 when a damaged
    sector is left (or, with `--dry-run`, when there is any damage) and 2 when there is no parity file. `--salvage TABLE` and
    `--reload TABLE` are the two steps of rebuilding a table with an unreadable clustered-index page
    (`planetgen.db.table_rebuild`); a reload ends with the repair above."""
    if args.salvage and args.reload:
        _fatal("--salvage and --reload are two steps with a restart between them; give one")
    config = store.mysql_config_from_args(args)
    scope = dbCheck.Scope(tuple(args.shell or ()), tuple(tuple(a) for a in args.sector or ()))
    directory = args.dump_dir or os.path.join(os.path.dirname(parity.parity_path()), "salvage")
    conn = store.get_connection(config)
    try:
        try:
            if args.salvage:
                table_rebuild.salvage(conn, args.salvage, directory)
                return
            if args.reload:
                table_rebuild.reload(conn, args.reload, directory)
                if args.restore_sectors and args.reload == "sectors":
                    table_rebuild.restore_sectors(conn)
            report = dbRepair.repair(conn, scope, dry_run=args.dry_run)
            if any(result.outcome == dbRepair.UNREADABLE for result in report.results):
                tables = table_rebuild.probe_unreadable(conn)
                log.normal("Repair: unreadable table(s): " + (", ".join(tables) or "none found now") + ". To rebuild one, stop the "
                           "website and generation, restart the server with innodb_force_recovery = 1, run `planetgen repair-db "
                           "--salvage TABLE`, restart it normally, then `planetgen repair-db --reload TABLE`.")
        except table_rebuild.RebuildRefused as exc:
            log.error(str(exc))
            sys.exit(1)
        except FileNotFoundError as exc:
            log.error(str(exc))
            sys.exit(2)
    finally:
        conn.close()
    verb = "would be repaired" if args.dry_run else "repaired"
    log.normal(f"Repair: {report.count(dbRepair.REPAIRED) + report.count(dbRepair.RESEALED) + report.count(dbRepair.WOULD_REPAIR):,}"
               f" sector(s) {verb}, {report.count(dbRepair.UNREPAIRABLE) + report.count(dbRepair.UNREADABLE):,} left damaged.")
    if not report.clean:
        sys.exit(1)


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
        prog='planetgen',
        description="Unified Generation CLI",
        epilog="Generates a star system, sector, galaxy, galaxy density skeleton, or exotic phenomenon, "
               "and saves it to the database -- see this script's module docstring for each subcommand's "
               "own section. Run 'planetgen <command> --help' for that command's own full option list.")
    parser.add_argument('--version', action=VersionAction, banner=version_banner('planetgen'))

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

    sectors_parser = subparsers.add_parser(
        'sectors', prefix_chars='-+',
        description="Generate sectors of the planned galaxy, in one of the modes: a shell, a band, one address, a span, "
                    "a course, the neighborhood of a sector or a random start. Only the sectors: 'planetgen link', "
                    "'backfill' and 'settle' (with --since, the start time the run says) do the stages after them.",
        help="Generate sectors of the planned galaxy (no link, scatter, backfill or settle).")
    add_shared_generation_options(sectors_parser)
    add_sectors_arguments(sectors_parser)

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

    scatter_parser = subparsers.add_parser(
        'scatter',
        description="Scatter Pre-placement: one kind of scatter on the stored plan",
        help="Draw one scatter (mass, luminosity, phenomena, nebula-field or band) on the stored plan.")
    add_scatter_arguments(scatter_parser)

    clear_scatter_parser = subparsers.add_parser(
        'clear-scatter',
        description="Empty one kind of scatter's stored rows, so 'planetgen scatter' can draw it again: the "
                    "massive stars, the bright stars (the band with them), the phenomena or the nebula field.",
        help="Empty one kind of scatter (mass, luminosity, phenomena or nebula-field).")
    add_clear_scatter_arguments(clear_scatter_parser)

    fill_map_luminosity_parser = subparsers.add_parser(
        'fill-map-luminosity',
        description="Fill the map luminosity of the scattered phenomena of a galaxy scattered before schema v84, "
                    "from their stored masses, in id-range batches that readers never wait for. Resumable: it only "
                    "fills rows still empty. A new scatter writes the value with each row.",
        help="Fill the map luminosity of the scattered phenomena.")
    add_fill_arguments(fill_map_luminosity_parser)

    rescore_habitability_parser = subparsers.add_parser(
        'rescore-habitability',
        description="Re-derive the stored habitability numbers (the flux swing of a planet in a double-star "
                    "system, the PHI-6 index, the six domains, the equipment tier and the note) of every planet and "
                    "moon from its stored conditions, so a galaxy generated before a change to the rules agrees with "
                    "a new one. Only those columns change; it can be stopped and run again.",
        help="Re-score the stored habitability of every planet and moon.")
    add_fill_arguments(rescore_habitability_parser)

    fill_nebula_blocks_parser = subparsers.add_parser(
        'fill-nebula-blocks',
        description="Make the nebula totals of the Galaxy Map's blocks again: for every block of the ladder, the "
                    "nebulae in it, their volume, centroid, families, dominant class and cover, from the nebula_field "
                    "table and the stored nebulae. A plan draws them with the field; this makes them for a galaxy "
                    "whose field was drawn before. Nothing else changes; it can be run again.",
        help="Make the nebula block totals of the Galaxy Map again.")
    add_fill_arguments(fill_nebula_blocks_parser)

    fill_block_totals_parser = subparsers.add_parser(
        'fill-block-totals',
        description="Make the totals of the Galaxy Map's blocks again, for every block of the ladder: what its "
                    "generated sectors hold (sectors, systems, stars, luminosity) and the light of the pre-placed "
                    "bright stars and scattered phenomena nothing was built from (count, summed luminosity, centre "
                    "of light, mean temperature). A scatter, a clear and the map luminosity fill make the light at "
                    "their end and a galaxy run brings in its sectors; this is for a galaxy whose totals were "
                    "never made, or after a run stopped before it finished. Nothing else changes; it can be run "
                    "again.",
        help="Make the block totals of the Galaxy Map again.")
    add_fill_arguments(fill_block_totals_parser)

    clean_nebula_field_parser = subparsers.add_parser(
        'clean-nebula-field',
        description="Bring an existing galaxy's stored molecular clouds in line with the current nebula keep "
                    "rates: counts the stored field clouds the rule now rejects and, with --apply, deletes them "
                    "(clearing what pointed at them and refreshing containment in the sectors they reached).",
        help="Count (and with --apply delete) the stored field clouds the nebula keep rate rejects.")
    add_clean_nebula_field_arguments(clean_nebula_field_parser)

    link_parser = subparsers.add_parser(
        'link',
        description="Link the galaxy sectors saved since a time to their neighbours: containment, nearest systems and "
                    "quadrants, and the new systems merged into the lists of the sectors around them. Nothing else "
                    "changes; it can be run again.",
        help="Link the sectors saved since a time to their neighbours.")
    add_since_arguments(link_parser)

    backfill_parser = subparsers.add_parser(
        'backfill',
        description="The bright-star backfill out from the galaxy sectors saved since a time: the sectors around them "
                    "get the stars born from 1, 2, 5 and 8 solar masses up, each only where it is not drawn yet.",
        help="Backfill the bright stars out from the sectors saved since a time.")
    add_since_arguments(backfill_parser)

    settle_parser = subparsers.add_parser(
        'settle',
        description="Save the path of every body in the galaxy sectors saved since a time and in the sectors around "
                    "them, so each path sees the final set of neighbours.",
        help="Save the sector paths of the sectors saved since a time.")
    add_since_arguments(settle_parser)
    settle_parser.add_argument('--sector-id', type=int, nargs='+', default=None, metavar='ID',
                               help="Instead of --since: save the paths of these sectors and of the sectors around "
                                    "them, for a sector a black hole or neutron star was added to ('planetgen "
                                    "phenomenon --sector-id' saves no paths itself).")

    run_job_parser = subparsers.add_parser(
        'run-job',
        description="Runs a Generate job as one run: the commands of SPEC, a JSON list of {\"command\", \"argv\"} "
                    "(argv: the command's options, what follows its name on its own command line), put together as the stages "
                    "of one run in the work queue, with one stage list and one estimate. A command's checks are made when "
                    "its own stage opens, so a job can plan the galaxy and scatter it. A job holds "
                    + ", ".join(generationJob.JOB_COMMANDS) + ".",
        help="Run a list of generation commands as one run.")
    run_job_parser.add_argument('spec', metavar='SPEC', help="The job as JSON, or @FILE to read it from a file.")
    run_job_parser.add_argument('--workers', type=int, default=None,
                                help="How many worker processes the run uses at most. Default: 80%% of this machine's cores "
                                     "(one fewer when MySQL runs here too), or PLANETGEN_WORKERS.")
    store.add_mysql_connection_args(run_job_parser)
    add_logging_arguments(run_job_parser)

    finish_run_parser = subparsers.add_parser(
        'finish-run',
        description="Carries on a generation run that stopped partway (a killed process, an out-of-memory database, a "
                    "cancel) from the first stage that did not finish, on the same galaxy and build, instead of "
                    "starting again: the stages it finished are kept, the one that stopped runs again from its start and "
                    "the rest run as they would have.",
        help="Carry on the last generation run that stopped before it was done.")
    finishRun.add_arguments(finish_run_parser)

    phenomenon_parser = subparsers.add_parser(
        'phenomenon',
        description="Exotic Stellar Phenomenon Generation Options",
        help="Generate a single exotic stellar phenomenon.")
    add_phenomenon_arguments(phenomenon_parser)

    check_math_parser = subparsers.add_parser(
        'check-math',
        description="Runs the math check (planetgen/physics/mathcheck.py): reference values from real "
                    "astronomy, identities and sampler distributions. Exits 1 if any check fails.",
        help="Check the generator's math (update.sh runs it).")
    check_math_parser.add_argument('-v', '--verbose', action='store_true',
                                   help="List every check, not only the failures.")
    add_logging_arguments(check_math_parser)

    versions_parser = subparsers.add_parser(
        'versions',
        description="Lists the version key each galaxy database was updated under, newest first.",
        help="List the version keys recorded by updates.")
    store.add_mysql_connection_args(versions_parser)
    add_logging_arguments(versions_parser)

    fingerprint_parser = subparsers.add_parser(
        'fingerprint',
        description="Prints a canonical SHA-256 digest of each sector's generated content, and one for the "
                    "region (the whole galaxy, plan included, without --shell or --sector), skipping row ids, "
                    "timestamps and what is rebuilt from positions (docs/design/generation/reproducible-galaxies.md). "
                    "Two builds of one galaxy are the same when their digests match. Writes nothing.",
        help="Print a digest of the generated content, to compare two builds.")
    add_fingerprint_arguments(fingerprint_parser)

    check_db_parser = subparsers.add_parser(
        'check-db',
        description="Checks the database for damage and changes nothing: schema revision, table health, "
                    "rows without parents, ids, impossible values, sector counts and version keys, and (named) "
                    "the thorough table check and every star system. Each check is a command of its own on the "
                    "work queue (planetgen commands). Ends with one pass or fail line per check; exits 1 when it "
                    "found damage and 2 when a check could not run.",
        help="Check the database for damage (writes nothing).")
    checkdb_cli.add_check_db_arguments(check_db_parser)
    add_logging_arguments(check_db_parser)

    population_parser = subparsers.add_parser(
        'population',
        description="Population and Politics Pass",
        help="Name species, date civilizations and found polities from what is stored.")
    add_population_arguments(population_parser)

    territories_parser = subparsers.add_parser(
        'territories',
        description="Territories Pass",
        help="Work out which polity owns which system from the polities stored.")
    add_territories_arguments(territories_parser)
    bring_to_epoch_parser = subparsers.add_parser(
        'bring-to-epoch',
        description="Brings the stored galaxy to the running generator epoch without a full regeneration: every "
                    "galaxy-placed sector made under an older epoch whose change altered what a sector holds is "
                    "generated again under the current code (its object IDs follow the galaxy's seed), the species "
                    "of an older population epoch are made again, and a sector with admin edits is left alone and "
                    "listed. A change of the galaxy's plan, or a scatter drawn under an older epoch, is refused with "
                    "the reason. An interrupted run is resumed by running it again.",
        help="Bring a stored galaxy to the running generator epoch.")
    add_bring_to_epoch_arguments(bring_to_epoch_parser)
    refresh_parity_parser = subparsers.add_parser(
        'refresh-parity',
        description="Makes the parity file that can rebuild a damaged sector (settings: parity) current for the whole "
                    "database: every sector that was marked or never sealed, every sealed sector the file does not hold, "
                    "and every member of the file whose sector is gone. Saves, edits and deletions keep it current "
                    "themselves once parity.enabled is on; this builds it for a galaxy that was made before, or catches "
                    "up after a crash.",
        help="Make the sector parity file current.")
    add_refresh_parity_arguments(refresh_parity_parser)

    repair_db_parser = subparsers.add_parser(
        'repair-db',
        description="Rebuilds the sectors whose rows no longer hash to their stored checksum from the parity file "
                    "(settings: parity; built by refresh-parity). A sector the parity cannot mend (its group has more "
                    "damage than parity blocks, or the record is behind the database too) is listed, with the reason. "
                    "The columns the orbit update rewrites are left as they are on every row that still exists. A parity "
                    "file that is itself damaged is made again from the database.",
        help="Rebuild damaged sectors from the parity file.")
    add_repair_db_arguments(repair_db_parser)

    benchmark_parser = subparsers.add_parser(
        'benchmark',
        description="Times a small galaxy from the plan to the finished fill in scratch databases and prints "
                    "the time and database work of every stage, at each worker count asked, with an optional "
                    "profile of the fill.",
        help="Time where generation spends its time.")
    benchmark.add_benchmark_arguments(benchmark_parser)
    add_logging_arguments(benchmark_parser)

    warm_parser = subparsers.add_parser(
        'warm',
        description="Builds the site's caches ahead of the first visitor: the Galaxy Map's opening view, the lists, "
                    "the sector, system and phenomenon pages and the search. Each cache family is its own command "
                    "(see 'planetgen commands'); this builds the list of them and the wrapper runs them in order "
                    "on the Redis queue. Run it as the web server's user.",
        help="Build the site's caches.")
    warm_cli.add_warm_arguments(warm_parser)
    add_logging_arguments(warm_parser)

    queue_parser = subparsers.add_parser(
        'queue',
        description="Stops everything and holds the host until told otherwise (stop-all: cancels every run and job, "
                    "empties the queues, ends every planetgen worker and command process; nothing starts until "
                    "resume), lifts the hold (resume) or shows it (status). Run it as root.",
        help="Stop everything, resume, or show the hold.")
    queue_cli.add_queue_arguments(queue_parser)
    add_logging_arguments(queue_parser)

    commands_parser = subparsers.add_parser(
        'commands',
        description="Lists every command the work queue can run: its name, whether it runs as units on the "
                    "queue or inline in the process that started the run, and what it does. Writes nothing.",
        help="List the work queue's commands.")
    add_logging_arguments(commands_parser)

    command_run_parser = subparsers.add_parser(
        'command-run',
        description="Runs one unit of a queue command alone, from its JSON arguments, and prints its JSON result: "
                    "the same function a worker runs, with nothing else around it. The unit's arguments name its "
                    "database and carry no password; the settings supply it.",
        help="Run one unit of a queue command from its JSON arguments.")
    command_run_parser.add_argument('kind', help="The command (see 'planetgen commands').")
    command_run_parser.add_argument('unit_args', metavar='ARGS',
                                    help="The unit's arguments as JSON, or @FILE to read them from a file.")
    add_logging_arguments(command_run_parser)

    return parser, {
        'commands': commands_parser,
        'command-run': command_run_parser,
        'warm': warm_parser,
        'queue': queue_parser,
        'benchmark': benchmark_parser,
        'system': system_parser,
        'sector': sector_parser,
        'sectors': sectors_parser,
        'galaxy': galaxy_parser,
        'plan': plan_parser,
        'scatter': scatter_parser,
        'clear-scatter': clear_scatter_parser,
        'fill-map-luminosity': fill_map_luminosity_parser,
        'rescore-habitability': rescore_habitability_parser,
        'fill-nebula-blocks': fill_nebula_blocks_parser,
        'fill-block-totals': fill_block_totals_parser,
        'clean-nebula-field': clean_nebula_field_parser,
        'link': link_parser,
        'backfill': backfill_parser,
        'settle': settle_parser,
        'run-job': run_job_parser,
        'phenomenon': phenomenon_parser,
        'finish-run': finish_run_parser,
        'population': population_parser,
        'territories': territories_parser,
        'bring-to-epoch': bring_to_epoch_parser,
        'refresh-parity': refresh_parity_parser,
        'repair-db': repair_db_parser,
        'check-math': check_math_parser,
        'versions': versions_parser,
        'fingerprint': fingerprint_parser,
        'check-db': check_db_parser,
    }


def process_args():
    """
    Parses command-line arguments for `planetgen`, then validates
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

    # --mysql-port's range (1 to 65535) is checked by its argparse type,
    # planetgen.db.store._mysql_port (OPS.6).

    validate_command(args, command_parser)
    if args.command == 'run-job':
        try:
            args.job = generationJob.parse_entries(args.spec, args, command_parsers, validate_command)
        except ValueError as exc:
            command_parser.error(str(exc))

    return args


def validate_command(args, command_parser):
    """
    Validates the arguments of the subcommand `args.command` through its own `validate_*_args` function(s), calling
    `command_parser.error` (which exits) on the first problem. A command with none needs none. A job (`run-job`) validates each
    of its commands' with this too.
    """
    if args.command == 'system':
        validate_system_args(args, command_parser)
    elif args.command == 'sector':
        validate_shared_generation_args(args, command_parser)
        validate_sector_args(args, command_parser)
    elif args.command == 'galaxy':
        validate_shared_generation_args(args, command_parser)
        validate_galaxy_args(args, command_parser)
    elif args.command == 'sectors':
        validate_shared_generation_args(args, command_parser)
        validate_sectors_args(args, command_parser)
    elif args.command == 'plan':
        validate_plan_args(args, command_parser)
    elif args.command == 'scatter':
        validate_scatter_args(args, command_parser)
    elif args.command in ('fill-map-luminosity', 'rescore-habitability', 'fill-nebula-blocks', 'fill-block-totals',
                          'clean-nebula-field', 'finish-run',
                          'link', 'backfill', 'settle'):
        validate_fill_args(args, command_parser)
    elif args.command == 'phenomenon':
        validate_phenomenon_args(args, command_parser)
    elif args.command == 'run-job' and args.workers is not None and args.workers < 0:
        command_parser.error("--workers must be 0 (automatic) or more.")


def run_fingerprint(args):
    """`planetgen fingerprint`: prints one line a sector (its address and
    digest), the plan's digest for the whole galaxy, and the region's
    (`planetgen.db.fingerprint`, GEN.58); with `--stored` the tree of digests
    stored at save time (`planetgen.db.fingerprint_tree`, GEN.305)."""
    conn = store.get_connection(store.mysql_config_from_args(args))
    try:
        if args.stored:
            result = fingerprintTree.lines(conn, shells=args.shell, addresses=args.sector)
        else:
            result = fingerprint.region_fingerprint(conn, shells=args.shell, addresses=args.sector).lines()
        warning = version_check.galaxy_warning(conn, compute_battery=True)
    finally:
        conn.close()
    if warning:
        log.normal(f"Warning: {warning}")
    for line in result:
        print(line)


def run_check_math(args):
    """`planetgen check-math`: runs the `check-math` command as a check run through the wrapper, prints its report and
    exits 1 if any check failed (TEST.68)."""
    from planetgen.queue import commands, redisqueue, runrecord, units, wrapper
    step = commands.Step("math", "check-math", "Check the math", {"verbose": bool(args.verbose)})
    found = []
    previous_term = _exit_on_sigterm()
    try:
        wrapper.run("Math check", [step], lambda index, plan: None, lambda index, unit, result: found.append(result),
                    kind=runrecord.CHECK, argv=sys.argv[1:])
    except (redisqueue.Unavailable, runrecord.Busy, units.RunStopped) as exc:
        log.error(f"error: the math check could not run: {exc}")
        raise SystemExit(1)
    finally:
        if previous_term is not None:
            signal.signal(signal.SIGTERM, previous_term)
    result = found[0]
    if result["failed"]:
        log.error(result["report"])
        raise SystemExit(1)
    log.normal(result["report"])


def run_versions(args):
    """`planetgen versions`: prints the version-key history (OPS.13)."""
    from planetgen.cli import version_history
    log.normal(version_history.format_history(store.mysql_config_from_args(args)))


def run_commands(args):
    """`planetgen commands`: lists the work queue's commands (PERF.118)."""
    for kind, command in sorted(queueCommands.REGISTRY.items()):
        log.normal(f"{kind:24} {'inline' if command.inline else 'units':7} {command.summary}")


def run_command_unit(args):
    """`planetgen command-run KIND ARGS`: runs one unit of a queue command alone and prints its JSON result."""
    text = args.unit_args
    if text.startswith("@"):
        with open(text[1:], encoding="utf-8") as handle:
            text = handle.read()
    try:
        unit_args = json.loads(text)
    except ValueError as exc:
        _fatal(f"the unit's arguments are not JSON: {exc}")
    try:
        result = queueCommands.run_alone(args.kind, unit_args)
    except KeyError:
        _fatal(f"no command {args.kind!r}; 'planetgen commands' lists them")
    except ValueError as exc:
        _fatal(str(exc))
    log.normal(json.dumps(result, sort_keys=True))


READ_ONLY_COMMANDS = ("check-db", "check-math", "fingerprint", "versions", "benchmark", "commands", "command-run", "warm", "queue",
                      "refresh-parity")
"""tuple: Subcommands that write nothing generated (`warm` only builds caches), so they aren't logged as runs
(no activity log lines or run history)."""


BULK_COMMANDS = ("galaxy", "plan", "scatter", "population", "territories", "finish-run", "run-job")
"""tuple: Subcommands that always generate in bulk, so the data check runs
first (GEN.590); `sector` joins them when it makes more than one sector
(`is_bulk_run`)."""


def is_bulk_run(args):
    """Whether this run is a bulk generation the data check must gate."""
    if args.command in BULK_COMMANDS:
        return True
    return args.command == "sector" and (getattr(args, "num_sectors", 1) or 1) > 1


_COMMAND_HANDLERS = {
    'benchmark': benchmark.run_benchmark,
    'commands': run_commands,
    'command-run': run_command_unit,
    'warm': warm_cli.run_warm,
    'queue': lambda args: sys.exit(queue_cli.run_queue(args) or 0),
    'check-math': run_check_math,
    'versions': run_versions,
    'check-db': checkdb_cli.run_check_db,
    'fingerprint': run_fingerprint,
    'system': run_system.run_system,
    'sector': run_sector.run_sector,
    'sectors': run_galaxy.run_sectors,
    'galaxy': run_galaxy.run_galaxy,
    'plan': run_plan.run_plan,
    'scatter': run_plan.run_scatter,
    'clear-scatter': run_plan.run_clear_scatter,
    'fill-map-luminosity': run_plan.run_fill_map_luminosity,
    'rescore-habitability': run_plan.run_rescore_habitability,
    'fill-nebula-blocks': run_plan.run_fill_nebula_blocks,
    'fill-block-totals': run_plan.run_fill_block_totals,
    'clean-nebula-field': run_plan.run_clean_nebula_field,
    'link': run_galaxy.run_link,
    'backfill': run_galaxy.run_backfill,
    'settle': run_galaxy.run_settle,
    'run-job': generationJob.run_job,
    'phenomenon': run_phenomenon.run_phenomenon,
    'finish-run': finishRun.run_finish,
    'population': run_population.run_population,
    'territories': run_population.run_territories,
    'bring-to-epoch': run_bring_to_epoch,
    'refresh-parity': run_refresh_parity,
    'repair-db': run_repair_db,
}


DATA_CHANGES = {"plan": ("everything",), "scatter": ("everything",), "clear-scatter": ("everything",),
                "fill-map-luminosity": ("everything",), "finish-run": ("everything",),
                "rescore-habitability": ("everything",), "fill-nebula-blocks": ("everything",), "fill-block-totals": ("everything",),
                "clean-nebula-field": ("everything",),
                "population": ("population",), "territories": ("population",),
                "bring-to-epoch": ("generated", "population"), "refresh-parity": (), "repair-db": ()}
"""dict: The kinds of data change (`planetgen.db.dataversion.CHANGES`) a command makes when it is not just
generation; any other writing command is `generated`."""


def _note_data_changed(args, database):
    """Tells the page cache what this run changed (PERF.104), so the pages read it at once instead of at the next
    check of the galaxy's content stamp. A run that stopped part-way changed what it had written."""
    commands = [command for command, _args in args.job] if args.command == "run-job" else [args.command]
    kinds = list(dict.fromkeys(kind for command in commands for kind in DATA_CHANGES.get(command, ("generated",))))
    for kind in kinds:
        dataversion.record(database, kind)


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
    elif args.debug is not None or log.config_debug():
        level = log.DEBUG       # --debug, or debug mode in config.json: every command writes its decisions to stdout
    else:
        level = log.NORMAL
    log.apply_record_arguments(args)
    try:
        # A command whose own result goes to stdout (`--output -`) keeps stdout for it; the debug lines go to stderr.
        log.configure(level, debug_file=(args.debug or None), to_stderr=getattr(args, "output", None) == "-",
                      areas=args.debug_areas)
    except OSError as exc:
        _fatal(f"cannot open --debug file {args.debug!r}: {exc.strerror or exc}", logger_ready=False)
    logged = (not getattr(args, "output", None) and args.command not in READ_ONLY_COMMANDS
              and not (args.command == "bring-to-epoch" and args.dry_run))
    with log.untraced():    # the seed read must not put SQL lines in front of the run line
        run_line = version_key.run_line(_run_line_seed(args) if logged else None, " ".join(_run_argv(sys.argv[1:])))
    log.normal(run_line)
    log.debug("Command: %s, options: %s", args.command,
              {key: ("<withheld>" if "password" in key else value) for key, value in sorted(vars(args).items())})

    seed = secrets.randbits(128)
    draw.set_run_seed(seed)
    log.debug(f"Seeded the run's random number generator with {seed}; galaxy sectors and bright stars draw "
              f"from the galaxy's own seed instead.")

    # GEN.590: a bulk run writes nothing at all (not even the activity log's
    # start line) when a downloaded data file is missing. The math check is
    # not run here (PERF.170): `update.sh` runs it at every update and the
    # admin page runs it on request.
    fetch.set_disregard_missing_data(getattr(args, "disregard_missing_data", False))
    if is_bulk_run(args):
        try:
            for line in run_galaxy.disregard_notice(run_galaxy.require_datasets()):
                log.normal(line)
        except run_galaxy.DatasetsMissing as exc:
            _fatal(str(exc))

    # One start and one finish line per run in the activity log (SEC.28);
    # a `system --output` run writes no database, so it isn't logged, nor
    # are `check-math` and `fingerprint`, which write nothing (`logged`, above).
    try:
        database = store.mysql_config_from_args(args).database
    except AttributeError:  # a subcommand without the --mysql-* options
        database = store.DEFAULT_MYSQL_CONFIG.database
    started = time.monotonic()
    if logged:
        activity_log.event("GEN", "generate.start", user=_run_user(), command=args.command, db=database)
    status = "failed"
    failure = None      # why the run's open stage failed (UX.122), when it did
    if logged:
        updateguard.refuse_while_updating("A run")      # OPS.53: nothing starts while the site is being updated
    previous_term = _exit_on_sigterm() if logged else None
    history = _start_history(args, seed) if logged and not getattr(args, "estimate_only", False) else None
    args.generation_run_id = history   # the run the per-object timings (PERF.92) are saved under
    try:
        _COMMAND_HANDLERS[args.command](args)
        status = "ok"
    except pymysql.err.MySQLError as exc:
        failure = f"database error: {exc}"
        _fatal(failure, traceback=True)
    except OSError as exc:
        # e.g. an unwritable --output path.
        where = f" ({exc.filename})" if getattr(exc, "filename", None) else ""
        failure = f"{exc.strerror or exc}{where}"
        _fatal(failure, traceback=True)
    except KeyboardInterrupt:
        status, failure = "interrupted", "interrupted"
        raise
    except SystemExit as exc:
        # SIGTERM inside a work queue (Cancel) exits with 128 + signal.
        if isinstance(exc.code, int) and exc.code >= 128:
            status, failure = "interrupted", "cancelled"
        elif exc.code in (None, 0):
            status = "ok"
        else:
            failure = failure or f"exited with status {exc.code}"
        raise
    except Exception as exc:  # noqa: BLE001 -- Python prints the traceback; at a terminal, hold it on screen first
        failure = f"{type(exc).__name__}: {exc}"
        _wait_for_enter()
        raise
    finally:
        if status != "ok":
            stages.fail(failure or "failed")
        if previous_term is not None:
            signal.signal(signal.SIGTERM, previous_term)
        if history is not None:
            _finish_history(args, history, status)
        if logged:
            _note_data_changed(args, database)
            activity_log.event("GEN", "generate.finish", user=_run_user(), command=args.command, db=database,
                              status=status, seconds=round(time.monotonic() - started, 1),
                              **{key: run_common.RUN_COUNTS[key] for key in ("sectors", "systems", "phenomena")})


def _exit_on_sigterm():
    """
    Makes a SIGTERM (Cancel on the Generate page, `kill`) end the run as `SystemExit(128 + 15)`, so the run's own
    cleanup runs (its run record ends as cancelled at once, not after the stale time) as Ctrl+C does.

    Returns:
        The handler it replaced, to put back when the run ends, or `None` when it is not the main thread.
    """
    if threading.current_thread() is not threading.main_thread():
        return None

    def on_term(signum, _frame):
        raise SystemExit(128 + signum)

    return signal.signal(signal.SIGTERM, on_term)


RUN_ARGV_WITHHELD = ("--mysql-host", "--mysql-port", "--mysql-user", "--mysql-password", "--mysql-database",
                     "--debug", "--debug-areas")
"""tuple: Options left out of the command line a run records (`_run_argv`): the database is the site's own when the admin
page runs it again, and a password never goes in the control database."""


def _run_argv(argv):
    """`argv` (after the script name) without `RUN_ARGV_WITHHELD`'s
    options and their values."""
    kept, skip = [], None
    for arg in argv:
        if skip == "value" or (skip == "optional" and not arg.startswith("-")):
            skip = None
            continue
        skip = None
        name = arg.split("=", 1)[0]
        if name in RUN_ARGV_WITHHELD:
            if "=" not in arg:
                skip = "optional" if name == "--debug" else "value"
            continue
        kept.append(arg)
    return kept


def _run_line_seed(args):
    """
    The galaxy seed the run's first line names (OPS.10): a `plan --seed`'s
    own, else the stored one. Read without touching the schema; `None`
    when there is none yet or the database can't be read (the run then
    reports that itself).
    """
    if getattr(args, "seed", None) is not None and args.command == "plan":
        return args.seed
    if args.command == "run-job":       # a job that plans the galaxy names the seed it plans with
        seeds = [entry.seed for command, entry in args.job if command == "plan" and getattr(entry, "seed", None) is not None]
        if seeds:
            return seeds[0]
    try:
        conn = store.get_connection(store.mysql_config_from_args(args), ensure_schema=False)
    except Exception as exc:  # noqa: BLE001 -- the run reports its own database errors
        log.debug(f"Run line: can't open the database ({exc}).")
        return None
    try:
        return store.get_galaxy_seed(conn)
    except Exception as exc:  # noqa: BLE001 -- no galaxy_shape table yet
        log.debug(f"Run line: no galaxy seed to show ({exc}).")
        return None
    finally:
        conn.close()


def _start_history(args, run_seed):
    """
    Records this run in the galaxy's run history (`generation_runs`, DB.6):
    its subcommand and command line (`_run_argv`), its own seed, and the
    code's version key. A galaxy is built by a series of runs, not by its
    seed alone, so this is what a rebuild replays. Never stops the run.

    Returns:
        int or None: The row's id, or `None` when it couldn't be written.
    """
    try:
        conn = store.get_connection(store.mysql_config_from_args(args))
    except Exception as exc:  # noqa: BLE001 -- the run reports its own database errors
        log.debug(f"Run history: can't open the database ({exc}).")
        return None
    try:
        return store.start_generation_run(conn, args.command, _run_argv(sys.argv[1:]), run_seed)
    except Exception as exc:  # noqa: BLE001 -- the history never fails a run
        log.debug(f"Run history: can't record this run ({exc}).")
        return None
    finally:
        conn.close()


def _finish_history(args, run_id, status):
    """Records how the run ended (`_start_history`); never raises."""
    try:
        conn = store.get_connection(store.mysql_config_from_args(args))
        try:
            store.finish_generation_run(conn, run_id, status)
        finally:
            conn.close()
    except Exception as exc:  # noqa: BLE001 -- the history never fails a run
        log.debug(f"Run history: can't record how this run ended ({exc}).")


def _run_user():
    """Who ran this: the login name (`getpass.getuser`), or `None`."""
    try:
        return getpass.getuser()
    except Exception:  # noqa: BLE001 -- no user name in this environment
        return None


def _fatal(message, logger_ready=True, traceback=False):
    """A one-line error and exit status 1, instead of a traceback unless
    `traceback` is set (ADM.25: the error being handled is shown in full,
    on the console, in a job's web log and in the debug log, so it can be
    copied into a report). Goes through `log.error` (shown even under
    --quiet/--silent) once the logger is configured, else straight to
    stderr. A run at a terminal then waits for Enter (ADM.24), so a
    console window that closes with the run doesn't take the error with
    it."""
    if logger_ready:
        if traceback:
            log.exception(f"Error: {message}")
        else:
            log.error(f"Error: {message}")
    else:
        print(f"planetgen: error: {message}", file=sys.stderr)
        if traceback:
            import traceback as traceback_module

            traceback_module.print_exc()
    _wait_for_enter()
    raise SystemExit(1)


def _wait_for_enter():
    """After a failed run at an interactive terminal, holds the output on
    screen until Enter is pressed (ADM.24); a run with its input or output
    redirected (a web job, a script, a test) exits at once."""
    if not (sys.stdin.isatty() and sys.stdout.isatty()):
        return
    try:
        input("The run failed. Press Enter to continue...")
    except (EOFError, KeyboardInterrupt):
        pass


if __name__ == "__main__":
    main()
