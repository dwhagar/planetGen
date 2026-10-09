# planetgen/cli/generate.py

"""
Generate
========

The single entry point for every generator, dispatched by subcommand
(run from the checkout's `src/` as `python3 -m planetgen.cli.generate`,
or as the installed `planetgen` command):

    planetgen system [options]      -- one star system
    planetgen sector [options]      -- one or more independent sectors
    planetgen galaxy [options]      -- many sectors placed as one galaxy
    planetgen plan [options]        -- the galaxy's density skeleton
    planetgen phenomenon [options]  -- one exotic stellar phenomenon
    planetgen population [options]  -- species, civilizations and territories
    planetgen check-math            -- the math check bulk runs start with

Run `planetgen <command> --help` for that command's own full option
list. This module is the command line: every command's options and
their validation, the parser, the run history, and `main`. Each
command's work is in `planetgen.generation.run_<command>`.
"""

import argparse
import getpass
import math
import secrets
import sys
import time

import pymysql

from planetgen.queue import work as workQueue
from planetgen.db import store
from planetgen.admin import activity_log
from planetgen.generation import limits, prevalence
from planetgen.galaxy import seed as galaxySeed, version_key
from planetgen.physics import mathcheck
from planetgen import tuning as program_constants
from planetgen.util import log
from planetgen._version import VersionAction, version_banner
from planetgen.galaxy.density import build_galaxy_shape
from planetgen.galaxy.drill import parse_drill_key
from planetgen.galaxy.skeleton import DEFAULT_MAX_RING
from planetgen.generation.star import STAR_TYPE_PATTERN
from planetgen.generation.star_population import bright_star_fraction
from planetgen.generation import run_common
from planetgen.generation import run_galaxy
from planetgen.generation import run_phenomenon
from planetgen.generation import run_plan
from planetgen.generation import run_population
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


STRICT_HELP = ("Stop instead of warning: a run the console would otherwise go ahead with after a "
               "warning (the database disk too small, a ring or block past "
               "LARGE_RING_WARNING_THRESHOLD sectors without --limit or --yes, an address outside the "
               "galaxy's outline, --min-habitable above a sector's drawn count, a forced body no "
               "system had room for) ends with an error and nothing generated, as before GEN.81.")


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

    if args.flavor_chance_system is not None and not (0.0 <= args.flavor_chance_system <= 1.0):
        parser.error("--flavor-chance-system must be a float between 0.0 and 1.0.")

    if args.flavor_chance_planet is not None and not (0.0 <= args.flavor_chance_planet <= 1.0):
        parser.error("--flavor-chance-planet must be a float between 0.0 and 1.0.")


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
    parser.add_argument('--prevalence', type=prevalence_setting, action='append', default=None,
                        metavar='FEATURE=PERCENT',
                        help="How much more or less often every system gets FEATURE than by chance, as a "
                             "percentage: comets=+50 is 1.5 times as often, comets=-100 never. Repeat "
                             "for more features: " + ", ".join(prevalence.FEATURES) + ".")
    parser.add_argument('--workers', type=int, default=None,
                        help="How many sectors to generate at once, each in its own low-priority worker "
                             "process. Default: 80%% of this machine's cores (one fewer when MySQL runs "
                             "here too), or PLANETGEN_WORKERS; 1 generates one sector at a time in this "
                             "process.")
    parser.add_argument('--population', action='store_true',
                        help="Also run the population pass (species, civilizations, territories) "
                             "after the sectors are saved. Off by default; 'planetgen population' "
                             "runs it any time.")
    parser.add_argument('--strict', action='store_true', help=STRICT_HELP)

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
                             "(exit status 1) if it is outside the galaxy's stored outline, or reports its "
                             "existing sector_id if it was already generated. A sparse address still "
                             "generates (it may hold no systems).")
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
                        help="Don't ask before generating: skips the size and time question asked on a "
                             "terminal before filling more than one sector, and with --ring (or --ring "
                             "--shell, or --block) the confirmation normally required before generating "
                             "more than LARGE_RING_WARNING_THRESHOLD sectors. A run the database disk "
                             "can't hold is still refused.")
    parser.add_argument('--estimate-only', action='store_true',
                        help="Only show the size and time estimate (and whether the disk would refuse the "
                             "run), then stop without writing anything; ends with one 'ESTIMATE {json}' "
                             "line. In random-start mode the estimate is for one random start.")
    parser.add_argument('--radius-pc', type=finite_float,
                        help="With --center-sector: the neighborhood search radius, in parsecs. With "
                             "--ring --slot: after generating that address, also generate its "
                             "neighborhood within this radius. With "
                             "neither --ring nor --center-sector (random-start mode): overrides the "
                             "default 12 pc neighborhood radius around the randomly chosen starting "
                             "sector. Once every sector is generated, the bright stars within 100 ly "
                             "of the requested sector are backfilled (see --backfill-from).")
    parser.add_argument('--max-ring', type=int,
                        help="With neither --ring nor --center-sector (random-start mode): the highest "
                             "ring the randomly chosen starting sector may land in. Default: anywhere "
                             "inside the galaxy's stored outline.")
    parser.add_argument('--min-start-density', type=finite_float,
                        help="With neither --ring nor --center-sector (random-start mode): require the "
                             "randomly chosen starting sector's own real relative_density (the same "
                             "'expected' figure printed alongside each saved sector) to be at least this "
                             "value before accepting it -- e.g. 1.0 for at least as dense as the galaxy's "
                             "own real local density. Retried the same way an already-occupied "
                             "address is (see RANDOM_START_MAX_PLACEMENT_ATTEMPTS). Cannot "
                             "be combined with --density/--num-systems (those override every position's "
                             "density uniformly, leaving no per-position value to compare against).")
    parser.add_argument('--no-settle', action='store_true',
                        help="Skip the last step of the run, which saves the path every star system, "
                             "rogue planet and comet takes through its sector (for the sectors the run "
                             "created and the ones around them). 'planetgen orbits' saves them later.")
    backfill_group = parser.add_argument_group("bright stars after the run (GEN.30)")
    backfill_group.add_argument('--backfill-from', choices=run_galaxy.BACKFILL_FROM_CHOICES, default="requested",
                                help="Once every sector of the run is generated, backfill the bright stars "
                                     "around the requested sector only ('requested', the default: the "
                                     "random start, --center-sector or --slot address, else the generated "
                                     "sector nearest the middle of the run), around every sector the run "
                                     "generated ('all', reaching 100 ly past the farthest one), or not at "
                                     "all ('none'). The backfill goes down to 100 L_sun under 10 ly, 250 "
                                     "under 25 ly, 500 under 50 ly and 750 out to 100 ly, and never adds "
                                     "stars to a sector already generated.")
    backfill_group.add_argument('--then-scatter', action='store_true',
                                help="After the sectors and before the backfill, scatter the bright stars "
                                     "galaxy-wide (as 'planetgen plan --bright-stars-only'), leaving out "
                                     "every sector already generated. A new galaxy uses this, so the "
                                     "scatter never draws stars for sectors it would fill anyway.")
    backfill_group.add_argument('--bright-star-min-luminosity', type=finite_float,
                                default=program_constants.BRIGHT_STAR_MIN_LUMINOSITY_SOL,
                                help="With --then-scatter: the scatter's threshold, solar luminosities. "
                                     f"Default: {program_constants.BRIGHT_STAR_MIN_LUMINOSITY_SOL:g}.")
    store.add_mysql_connection_args(parser)


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
        if block_layer is not None and block_layer not in run_galaxy.block_layers(args.block):
            layers = run_galaxy.block_layers(args.block)
            parser.error(f"--block-layer must be one of the block's layers, {layers[0]} to {layers[-1]}.")
        if args.limit is not None and not 1 <= args.limit <= limits.MAX_GENERATE_LIMIT:
            parser.error(f"--limit must be between 1 and {limits.MAX_GENERATE_LIMIT}.")
        args.sector_name = args.system_file = args.num_orbits = args.name = None
        return
    if block_layer is not None:
        parser.error("--block-layer requires --block.")

    random_start = args.ring is None and args.center_sector is None

    if args.ring is not None and args.ring < 0:
        parser.error("--ring must be >= 0.")
    if args.ring is not None and args.ring > limits.MAX_GENERATE_RING:
        parser.error(f"--ring must be at most {limits.MAX_GENERATE_RING}.")
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
    if args.radius_pc is not None and args.radius_pc > limits.MAX_GENERATE_RADIUS_PC:
        parser.error(f"--radius-pc must be at most {limits.MAX_GENERATE_RADIUS_PC:g}.")

    if args.limit is not None and args.ring is None:
        parser.error("--limit only applies to --ring.")
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be a positive integer.")
    if args.limit is not None and args.limit > limits.MAX_GENERATE_LIMIT:
        parser.error(f"--limit must be at most {limits.MAX_GENERATE_LIMIT}.")
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
    if args.max_ring is not None and args.max_ring > limits.MAX_GENERATE_RING:
        parser.error(f"--max-ring must be at most {limits.MAX_GENERATE_RING}.")

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
    shape_group = parser.add_argument_group("galaxy shape (galaxyDensity.GalaxyShape)")
    shape_group.add_argument('--disk-scale-length-pc', type=finite_float, default=2600.0,
                             help="Thin disk exponential scale length, parsecs (the thick disk's is 0.77x). "
                                  "Default: 2600 "
                                  "(real Milky Way scale).")
    shape_group.add_argument('--disk-scale-height-pc', type=finite_float, default=300.0,
                             help="Thin disk exponential scale height, parsecs (the thick disk's is 3x). "
                                  "Default: 300 (real Milky Way).")
    shape_group.add_argument('--bulge-scale-radius-pc', type=finite_float, default=1580.0,
                             help="Bar bulge scale length along the bar, parsecs (0.39x across it, 0.27x "
                                  "vertically). Default: 1580 (COBE/DIRBE fit, Dwek et al. 1995).")
    shape_group.add_argument('--bulge-amplitude', type=finite_float, default=3.11,
                             help="Bulge central density, relative to the thin disk's at the center. "
                                  "Default: 3.11 (bulge 31%% of the stars, as in the Milky Way).")
    shape_group.add_argument('--arm-count', type=int, default=2,
                             help="Number of spiral arms. Default: 2 (grand-design).")
    shape_group.add_argument('--pitch-angle-deg', type=finite_float, default=15.0,
                             help="Spiral arm pitch angle, degrees. Default: 15.")
    shape_group.add_argument('--arm-amplitude', type=finite_float, default=0.4,
                             help="Arm/inter-arm density contrast amplitude, in [0, 1). Default: 0.4.")
    shape_group.add_argument('--calibration-radius-pc', type=finite_float, default=None,
                             help="In-plane radius the relative_density=1.0 calibration point sits at. "
                                  "Defaults to build_galaxy_shape's own default (3.15x disk scale length, the Sun's radius).")

    parser.add_argument('--seed', type=_galaxy_seed_arg, default=None, metavar="HEX",
                        help="The galaxy's 128-bit seed, as 32 hex digits: the same seed makes the same galaxy "
                             "on the same PlanetGen release. Default: the seed already stored, or one drawn at "
                             "random for a new galaxy. Changing it needs a galaxy with no sectors.")
    parser.add_argument('--max-ring', type=int, default=DEFAULT_MAX_RING,
                        help=f"Hard cap on how far out a layer is scanned. Default: {DEFAULT_MAX_RING}.")

    bright_group = parser.add_argument_group("bright-star pre-placement (planetgen.generation.bright_stars)")
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
                              help="No longer needed: the scatter always leaves filled sectors out (GEN.30).")
    bright_group.add_argument('--workers', type=int, default=None,
                              help="How many layers of bright stars to draw at once, each in its own "
                                   "low-priority worker process. Default: 80%% of this machine's cores "
                                   "(one fewer when MySQL runs here too), or PLANETGEN_WORKERS.")
    bright_group.add_argument('--bright-stars-down-to', type=finite_float, default=None, metavar='L_SUN',
                              help="Go one layer dimmer on the stored plan: keep the bright stars already "
                                   "placed and add only those from L_SUN up to the level already scattered. "
                                   "Sectors already filled are left out (their own systems already reach "
                                   "that bright). Does nothing when L_SUN is not below the current level.")
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
    if args.arm_amplitude < 0 or args.arm_amplitude >= 1:
        parser.error("--arm-amplitude must be in [0, 1).")
    if args.workers is not None and args.workers < 0:
        parser.error("--workers must be 0 (automatic) or more.")
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
    if args.bright_stars_down_to is not None:
        if args.no_bright_stars or args.bright_stars_only:
            parser.error("--bright-stars-down-to can't be combined with --no-bright-stars or --bright-stars-only.")
        try:
            bright_star_fraction(args.bright_stars_down_to)
        except ValueError as exc:
            parser.error(f"--bright-stars-down-to: {exc}")
    if not args.no_bright_stars:
        try:
            bright_star_fraction(args.bright_star_min_luminosity)
        except ValueError as exc:
            parser.error(f"--bright-star-min-luminosity: {exc}")


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
                              "planetgen.db.store.compute_phenomenon_placement). Also links a "
                              "supernova-remnant/rogue-planet/comet to that sector, without a computed "
                              "galaxy position (those types have no placement columns of their own). A "
                              "quasar needs a ring-0, layer-0 (galactic core) sector and is placed at the galactic "
                              "center; only one per galaxy. Omit to generate it unplaced/unlinked, as before.")

    parser.add_argument('--no-settle', action='store_true',
                        help="With --sector-id: skip saving the sector paths of that sector and the sectors "
                             "around it, which a black hole or neutron star can bend ('planetgen orbits' "
                             "saves them later).")

    # Database persistence
    store.add_mysql_connection_args(parser)

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
    if args.num_orbits is not None and args.num_orbits > limits.MAX_NUM_ORBITS:
        parser.error(f"--num-orbits must be at most {limits.MAX_NUM_ORBITS}.")
    if args.anchor_system and args.type not in (None, "black-hole", "neutron-star"):
        parser.error("--anchor-system is only valid with --type black-hole or --type neutron-star.")
    if args.sector_id is not None and args.anchor_system:
        parser.error("--sector-id cannot be combined with --anchor-system -- an anchored compact remnant "
                     "belongs to its own StarSystem, not directly to a sector.")


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
    store.add_mysql_connection_args(parser)
    add_logging_arguments(parser)


def validate_population_args(args, parser):
    """`--rescan` and `--territories-only` contradict each other."""
    if args.rescan and args.territories_only:
        parser.error("--rescan and --territories-only cannot be combined.")


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

    check_math_parser = subparsers.add_parser(
        'check-math',
        description="Runs the math check (planetgen/physics/mathcheck.py): reference values from real "
                    "astronomy, identities and sampler distributions. Exits 1 if any check fails.",
        help="Check the generator's math before generating.")
    check_math_parser.add_argument('-v', '--verbose', action='store_true',
                                   help="List every check, not only the failures.")
    add_logging_arguments(check_math_parser)

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
        'check-math': check_math_parser,
    }


def reject_lost_option_values(args, parser):
    """
    Python 3.9's argparse drops a lone `--` given as an option's value
    (`--workers=--`) and stores an empty list instead of an error
    (PERF.27; later versions fix it). Any single-value option holding a
    list it didn't default to gets the usage error newer versions give.

    Args:
        args (argparse.Namespace): The parsed arguments.
        parser (argparse.ArgumentParser): The command's own parser.

    Raises:
        SystemExit: Through `parser.error`, exit status 2.
    """
    for action in parser._actions:
        if type(action) is not argparse._StoreAction or action.nargs is not None or not action.option_strings:
            continue
        value = getattr(args, action.dest, None)
        if isinstance(value, list) and not isinstance(action.default, list):
            parser.error(f"argument {'/'.join(action.option_strings)}: expected one argument")


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

    reject_lost_option_values(args, command_parser)
    validate_logging_args(args, command_parser)

    # --mysql-port's range (1 to 65535) is checked by its argparse type,
    # planetgen.db.store._mysql_port (OPS.6).

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


def run_check_math(args):
    """`planetgen check-math`: prints the math check's report and exits
    1 if any check failed (TEST.68)."""
    results = mathcheck.run_all()
    report = mathcheck.format_report(results, verbose=args.verbose)
    if mathcheck.failures(results):
        log.error(report)
        raise SystemExit(1)
    log.normal(report)


BULK_COMMANDS = ("galaxy", "plan", "population")
"""tuple: Subcommands that always generate in bulk, so the math check runs
first (TEST.68); `sector` joins them when it makes more than one sector
(`is_bulk_run`)."""


def is_bulk_run(args):
    """Whether this run is a bulk generation the math check must gate."""
    if args.command in BULK_COMMANDS:
        return True
    return args.command == "sector" and (getattr(args, "num_sectors", 1) or 1) > 1


_COMMAND_HANDLERS = {
    'check-math': run_check_math,
    'system': run_system.run_system,
    'sector': run_sector.run_sector,
    'galaxy': run_galaxy.run_galaxy,
    'plan': run_plan.run_plan,
    'phenomenon': run_phenomenon.run_phenomenon,
    'population': run_population.run_population,
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
    logged = not getattr(args, "output", None) and args.command != "check-math"
    log.normal(version_key.run_line(_run_line_seed(args) if logged else None, " ".join(_run_argv(sys.argv[1:]))))
    log.debug("Command: %s, options: %s", args.command,
              {key: ("<withheld>" if "password" in key else value) for key, value in sorted(vars(args).items())})

    seed = secrets.randbits(128)
    draw.set_run_seed(seed)
    log.debug(f"Seeded the run's random number generator with {seed}; galaxy sectors and bright stars draw "
              f"from the galaxy's own seed instead (GEN.39).")

    # TEST.68: a bulk run checks the math first and writes nothing at all
    # (not even the activity log's start line) when it fails.
    if is_bulk_run(args):
        try:
            run_galaxy.require_math_check()
        except run_galaxy.MathCheckFailed as exc:
            _fatal(str(exc))

    # One start and one finish line per run in the activity log (SEC.28);
    # a `system --output` run writes no database, so it isn't logged, nor
    # is `check-math`, which writes nothing (`logged`, above).
    try:
        database = store.mysql_config_from_args(args).database
    except AttributeError:  # a subcommand without the --mysql-* options
        database = store.DEFAULT_MYSQL_CONFIG.database
    started = time.monotonic()
    if logged:
        activity_log.event("GEN", "generate.start", user=_run_user(), command=args.command, db=database)
    status = "failed"
    root = _open_run_node(args, database) if logged else None
    history = _start_history(args, seed) if logged and not getattr(args, "estimate_only", False) else None
    try:
        _COMMAND_HANDLERS[args.command](args)
        status = "ok"
    except pymysql.err.MySQLError as exc:
        _fatal(f"database error: {exc}", traceback=True)
    except OSError as exc:
        # e.g. an unwritable --output path.
        where = f" ({exc.filename})" if getattr(exc, "filename", None) else ""
        _fatal(f"{exc.strerror or exc}{where}", traceback=True)
    except KeyboardInterrupt:
        status = "interrupted"
        raise
    except SystemExit as exc:
        # SIGTERM inside a work queue (Cancel) exits with 128 + signal.
        if isinstance(exc.code, int) and exc.code >= 128:
            status = "interrupted"
        elif exc.code in (None, 0):
            status = "ok"
        raise
    except Exception:  # noqa: BLE001 -- Python prints the traceback; at a terminal, hold it on screen first
        _wait_for_enter()
        raise
    finally:
        if root is not None:
            workQueue.close_node(root, {"ok": "done", "interrupted": "cancelled"}.get(status, "failed"))
        if history is not None:
            _finish_history(args, history, status)
        if logged:
            activity_log.event("GEN", "generate.finish", user=_run_user(), command=args.command, db=database,
                              status=status, seconds=round(time.monotonic() - started, 1),
                              **{key: run_common.RUN_COUNTS[key] for key in ("sectors", "systems", "phenomena")})


RUN_ARGV_WITHHELD = ("--mysql-host", "--mysql-port", "--mysql-user", "--mysql-password", "--mysql-database",
                     "--debug")
"""tuple: Options left out of the command line a run records in its job
tree root (`_run_argv`): the database is the site's own when the admin
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


def _open_run_node(args, database):
    """The job tree root of this run (ADM.12): its command, what it was
    asked to do, and the command line that would run it again. A run the
    Generate page started goes under that page job's step
    (`workQueue.PARENT_ENV_VAR`). Recorded in the control database when
    there is one; never stops the run."""
    argv = _run_argv(sys.argv[1:])
    title = " ".join(["planetgen", *argv])[:255]
    try:
        control = store.control_mysql_config(store.mysql_config_from_args(args))
    except AttributeError:  # a subcommand without the --mysql-* options
        control = store.control_mysql_config()
    return workQueue.open_node(args.command, title, control, argv=argv, database=database)


def _run_line_seed(args):
    """
    The galaxy seed the run's first line names (OPS.10): a `plan --seed`'s
    own, else the stored one. Read without touching the schema; `None`
    when there is none yet or the database can't be read (the run then
    reports that itself).
    """
    if getattr(args, "seed", None) is not None and args.command == "plan":
        return args.seed
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
