import argparse
import copy
import logging
import os
import random
import secrets
import sys

# stellarObjects lives at src/stellarObjects (src layout) -- add src/ to the
# import path so this keeps working without requiring `pip install .`
# first, matching how html/'s CGI scripts fall back to a no-install layout.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))

import phenomenonGen
import systemGen
from stellarObjects import _db, program_constants
from stellarObjects._version import VersionAction, version_banner
from stellarObjects.asteroidFieldData import AsteroidField
from stellarObjects.compactRemnant import BlackHole, NeutronStar
from stellarObjects.config import SystemConfig
from stellarObjects.names import SECTOR_NAMES, SECTOR_PREFIXES, SECTOR_SUFFIXES
from stellarObjects.nebulaData import Nebula
from stellarObjects.roguePlanetData import InterstellarComet, RoguePlanet
from stellarObjects.spaceSector import SpaceSector, _sample_poisson_count
from stellarObjects.supernovaRemnantData import SupernovaRemnant
from stellarObjects.systemData import StarSystem
from stellarObjects.utils import generate_phoneme_salad_name

# Suppress transformers warnings
logging.getLogger("transformers").setLevel(logging.ERROR)

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
    Adds every per-system generation-tuning option this script and
    `galaxyGen.py` share -- everything `build_sector_configs()` (and, via
    it, `systemGen.build_system_config`) needs from parsed args -- to
    `parser`. Factored out so the two scripts' CLI surfaces can never
    silently drift apart on what a given flag means, per this track's
    "call into sectorGen.py rather than duplicating its argparse logic"
    brief.

    Deliberately excludes `--name`/`-n` (sector naming is this script's
    own per-invocation concept; `galaxyGen.py` names each generated sector
    itself, once per shell slot / neighborhood address) and `--output`/
    `--db-path` (each script's own I/O conventions differ enough that
    sharing them would obscure more than it saves).

    Args:
        parser (argparse.ArgumentParser): The parser to add options to.
                                          Must accept `+`/`-` prefix chars
                                          (`prefix_chars='-+'`) for the
                                          tri-state options below.
    """
    for name, _attr, description in systemGen.TRISTATE_OPTIONS:
        parser.add_argument(f'-{name}', f'+{name}', dest=name, action=systemGen.TristateAction,
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


def validate_shared_generation_args(args, parser):
    """
    Validates every option `add_shared_generation_options` added, calling
    `parser.error` (which exits) on the first problem found. Shared with
    `galaxyGen.py` for the same reason `add_shared_generation_options` is.

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


def process_args():
    """
    Parses command-line arguments for generating a whole sector of star
    systems at once.

    Every option `systemGen.py` accepts for tuning a single system's
    generation -- the `+name`/`-name` tri-state flags (see
    `systemGen.TRISTATE_OPTIONS`), `--star-type`, `--age`, `--markdown`, and
    the flavor-text overrides -- applies here too, uniformly, to every
    system the sector contains (see `add_shared_generation_options`).
    `--system-file`, `--num-orbits`, and systemGen.py's own per-system
    `--name` are deliberately not offered here: those describe one
    specific, hand-crafted system (exact orbital slots, an exact object
    count, a single fixed name), which contradicts generating a whole
    sector of varied, independently-random systems. Use
    `systemGen.py --system-file` directly for that, and stitch its output
    into a sector by hand if needed. (This script's own `--name`/`-n`,
    below, is a different setting entirely -- it names the *sector*, not
    any one system in it.)

    Sector-specific options, on top of everything reused from `systemGen.py`:
    - `--num-systems`: How many star systems the sector contains. Defaults
      to 10.
    - `--name` / `-n`: Hard-sets the sector's own name, overriding the
      default random two-word name (see `generate_sector_name`). Stored
      under `args.sector_name`, not `args.name` -- `args.name` is reserved
      for `systemGen.build_system_config`'s per-*system* forced-name option,
      deliberately left `None` below (see that assignment's comment); the
      sector's name and an individual system's name are unrelated settings
      that happen to share an obvious flag spelling.
    - `--min-habitable`: Guarantees at least this many of the sector's
      systems have a habitable world, chosen randomly among them, without
      forcing *every* system to have one the way a uniform
      `+habitable_world` would. Cannot exceed `--num-systems`, and cannot be
      combined with a uniform `-habitable_world` (which forbids habitable
      worlds sector-wide).
    - `--density`: An alternative to `--num-systems` -- a multiplier on
      real local stellar density (see `SpaceSector.expected_system_count`)
      that gets randomly (Poisson) sampled into a concrete count per
      sector, so different invocations -- or different sectors within one
      `--num-sectors` run -- can be meaningfully denser or sparser than
      each other rather than always generating the same flat count.
      Cannot be combined with `--num-systems`.
    - `--num-sectors`: Generates this many independent sectors in one run,
      each with no galactic positioning (unlike `galaxyGen.py`), all saved
      into the same database. Defaults to 1. Cannot be combined with
      `--name`/`-n`, since every generated sector would otherwise share the
      same forced name.
    - `--console`: Prints each sector's rendered Markdown/wikitext to the
      console. Off by default -- without it, `main()` only prints a short
      status line per saved sector, regardless of `--num-sectors`.

    Returns:
        argparse.Namespace: An object containing the parsed command-line arguments.
    """
    additional_info = [
        "Additional Information:",
        "This tool generates a whole sector of independently-random star systems in one pass, reusing systemGen.py's",
        "own generation logic and options for each one. Options like +habitable_world/-habitable_world, --star-type,",
        "and --age apply uniformly to every system in the sector; use --min-habitable instead if you just want a",
        "guaranteed number of habitable systems among an otherwise varied sector. Use --density instead of",
        "--num-systems for a physically-grounded, meaningfully denser-or-sparser-than-another sector. For one",
        "specific, hand-crafted system (exact orbital slots, an exact name), use systemGen.py --system-file",
        "directly instead."
    ]
    additional_info = " ".join(additional_info)

    parser = argparse.ArgumentParser(
        description="Sector Generation Options",
        epilog=additional_info,
        prefix_chars='-+')

    parser.add_argument('--version', action=VersionAction, banner=version_banner('sectorGen.py'))

    add_shared_generation_options(parser)

    parser.add_argument('--name', '-n', dest='sector_name', type=str,
                        help="Force the name of the sector, overriding the default random two-word name. "
                             "Cannot be combined with --num-sectors > 1.")
    parser.add_argument('--num-sectors', type=int, default=1,
                        help="Generate this many independent sectors, each with no galactic positioning, "
                             "saving all of them into the same database. Defaults to 1.")
    parser.add_argument('--output', '-o', type=str, help="Output to a file.")
    parser.add_argument('--console', action='store_true',
                        help="Also print each sector's rendered Markdown/wikitext to the console. By "
                             "default, only status messages are printed.")
    _db.add_mysql_connection_args(parser)

    args = parser.parse_args()

    validate_shared_generation_args(args, parser)

    if args.num_sectors < 1:
        parser.error("--num-sectors must be a positive integer.")

    if args.num_sectors > 1 and args.sector_name:
        parser.error("--name/-n cannot be combined with --num-sectors > 1 (every generated sector would "
                     "share the same forced name).")

    # systemGen.build_system_config() expects a namespace shaped like its own
    # process_args() output, including these three -- deliberately not
    # exposed as sector-level flags (see this function's docstring), so they
    # get the same "not given" default systemGen.py's own parser would.
    args.system_file = None
    args.num_orbits = None
    args.name = None

    return args


def generate_sector_name():
    """
    Generates a random two-word sector name, each word independently drawn
    from the same phoneme-salad name generator used for star/planet/moon
    names -- using the sector-flavored `SECTOR_NAMES`/`SECTOR_PREFIXES`/
    `SECTOR_SUFFIXES` base lists instead, so generated sectors draw on real
    astronomical regions (galactic arms, superclusters, nebulae) and
    science-fiction sector names rather than reusing star names verbatim.
    No literal "Sector" suffix. Overridden entirely by `--name`/`-n` (see
    `process_args`), which hard-sets the whole name instead.

    Returns:
        str: A newly generated sector name, e.g. "Voranthis Kelmoor" --
        always exactly two words.
    """
    # allow_split=False: generate_phoneme_salad_name can itself split a
    # long result into two words (e.g. "Xyleth Anore"). Since this
    # function already joins two independent calls into one name, leaving
    # splitting on could silently produce 3-4 words instead of 2.
    # syllable_fraction=0.5 trims each word's base syllables by about
    # half before the prefix/suffix are attached -- many SECTOR_NAMES
    # entries (e.g. "Sagittarius", "Metropolis") are long real place
    # names, and two of them joined together made for unwieldy sector
    # names. max_length=7 backstops that: prefixes, suffixes, and the
    # occasional spliced-in universal phoneme are fixed-ish overhead that
    # doesn't shrink with syllable_fraction, so a long base name could
    # still slip through longer than intended without a hard cap too.
    first_word = generate_phoneme_salad_name(SECTOR_NAMES, SECTOR_PREFIXES, SECTOR_SUFFIXES, allow_split=False, syllable_fraction=0.5, max_length=7)
    second_word = generate_phoneme_salad_name(SECTOR_NAMES, SECTOR_PREFIXES, SECTOR_SUFFIXES, allow_split=False, syllable_fraction=0.5, max_length=7)
    return f"{first_word} {second_word}"


def build_sector_configs(args):
    """
    Builds one `SystemConfig` per system in the sector, sharing the same
    tri-state/value options across all of them (via
    `systemGen.build_system_config`, so this can never silently drift from
    what those options mean for a single system), then -- if
    `--min-habitable` was given and not already guaranteed by a uniform
    `+habitable_world` -- forces `HABITABLE_WORLD = True` on that many
    randomly-chosen configs among the rest.

    Args:
        args (argparse.Namespace): Parsed arguments from `process_args()`,
            with `args.num_systems` already resolved to a concrete count --
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
    configs = [systemGen.build_system_config(args)[0] for _ in range(args.num_systems)]

    if args.min_habitable > len(configs):
        raise SystemExit(
            f"Error: --min-habitable ({args.min_habitable}) exceeds this sector's generated system "
            f"count ({len(configs)}); with --density, the count is randomly sampled per sector and can "
            f"land below --min-habitable. Try a smaller --min-habitable, a higher --density, or "
            f"--num-systems for an exact count instead."
        )

    if args.min_habitable > 0:
        already_habitable = [i for i, cfg in enumerate(configs) if cfg.HABITABLE_WORLD is True]
        still_needed = args.min_habitable - len(already_habitable)
        if still_needed > 0:
            candidates = [i for i in range(len(configs)) if i not in already_habitable]
            for i in random.sample(candidates, k=still_needed):
                configs[i].HABITABLE_WORLD = True
                # Mirrors systemGen.build_system_config's own habitable-world +
                # asteroid-belt normalization, reapplied here since it ran before
                # this override existed.
                if configs[i].ASTEROID_BELT is True:
                    if configs[i].LARGE_STAR is False:
                        raise SystemExit(
                            "Error: --min-habitable requires forcing a habitable world onto a system that also "
                            "has +asteroid_belt forced sector-wide; that combination needs a large star, but "
                            "-large_star was also forced sector-wide."
                        )
                    configs[i].LARGE_STAR = True

    return configs


def generate_sector_phenomena(sector, args, galactic_center_dist_ly=None):
    """
    Populates an already-built `sector` with a realistically sparse
    population of exotic stellar phenomena (`phenomenonGen.py`'s own seven
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
    rendering, printing, or saving anything -- the shared core `main()`
    and `galaxyGen.py` both build on, per this track's "expose a clean
    function-level entry point for building and saving one sector at a
    given galaxy position" requirement.

    Args:
        args (argparse.Namespace): Parsed arguments from `process_args()`
            (or an equivalently-shaped namespace a caller like
            `galaxyGen.py` builds itself -- see
            `add_shared_generation_options`/`validate_shared_generation_args`
            for the option surface this function actually reads, via
            `build_sector_configs`).
        galactic_center_dist_ly (float, optional): This sector's distance
            from the galactic center, in light-years, threaded down into
            every generated system's Hill-sphere calculation (see
            `Star.calculate_system_perimeter`'s docstring). `None` (the
            default) falls back to the fixed `GALACTIC_CENTER_DISTANCE_LY`
            constant -- this script's own standalone CLI (`main()`, no
            galaxy context) always calls this with the default, so it
            keeps producing "unplaced" sectors exactly as before.

    Returns:
        tuple: `(sector_name, SpaceSector)` -- `sector_name` is
              `args.sector_name` or a freshly generated one (see
              `generate_sector_name`); the `SpaceSector` already has every
              system added (`SpaceSector.add_system`'s Hill-sphere-based
              random placement), plus a realistically sparse population of
              exotic phenomena (see `generate_sector_phenomena`).
    """
    sector_name = args.sector_name or generate_sector_name()
    sector = SpaceSector(name=sector_name)

    # `--density` resolves to a concrete system count per sector (this
    # sector's own volume, sampled fresh each call) rather than once at
    # parse time -- this is what lets each sector under `--num-sectors`/
    # `galaxyGen.py` vary independently instead of sharing one fixed count.
    # A copy avoids mutating the caller's shared `args` namespace, since
    # this function runs once per sector.
    if args.density is not None:
        args = copy.copy(args)
        args.num_systems = _sample_poisson_count(sector.expected_system_count() * args.density)

    configs = build_sector_configs(args)
    systems = [
        StarSystem(system_config=cfg, galactic_center_dist_ly=galactic_center_dist_ly)
        for cfg in configs
    ]

    for system, cfg in zip(systems, configs):
        sector.add_system(system, system_config=cfg)

    generate_sector_phenomena(sector, args, galactic_center_dist_ly=galactic_center_dist_ly)

    return sector_name, sector


def render_sector_text(sector_name, systems, phenomena, markdown):
    """
    Renders one generated sector's systems (and any generated exotic
    phenomena) as a Markdown or wikitext blob: a sector header, a short
    summary line, an index of every system's name and star type (plus,
    when present, every phenomenon's name and type), then each system's
    and phenomenon's own full rendering, all divided per `markdown`.
    Factored out of `main()` so it can be skipped entirely when neither
    `--console` nor `--output` was given -- rendering a large sector isn't
    free, and by default this script only reports status.

    Args:
        sector_name (str): The sector's name.
        systems (list): The sector's `StarSystem` instances.
        phenomena (list): The sector's `SectorPhenomenonEntry` instances
            (see `generate_sector_phenomena`) -- realistically empty for
            most sectors.
        markdown (bool): `True` for Markdown, `False` for MediaWiki wikitext.

    Returns:
        str: The fully rendered sector text.
    """
    habitable_count = sum(1 for s in systems if s.hab_count > 0)

    output_parts = []
    if markdown:
        output_parts.append(f"# {sector_name}\n\n")
    else:
        output_parts.append(f"= {sector_name} =\n\n")

    system_word = "system" if len(systems) == 1 else "systems"
    habitable_verb = "harbors" if habitable_count == 1 else "harbor"
    output_parts.append(
        f"This sector contains {len(systems)} star {system_word}, "
        f"{habitable_count} of which {habitable_verb} a potentially habitable world.\n\n"
    )

    bullet = "-" if markdown else "*"
    for system in systems:
        output_parts.append(f"{bullet} {system.star.name} ({system.star.type})\n")
    output_parts.append("\n")

    if phenomena:
        phenomenon_word = "phenomenon" if len(phenomena) == 1 else "phenomena"
        output_parts.append(f"It also holds {len(phenomena)} notable stellar {phenomenon_word}:\n\n")
        for entry in phenomena:
            label = phenomenonGen.TYPE_LABELS[entry.phenomenon_type]
            output_parts.append(f"{bullet} {entry.phenomenon.name} ({label})\n")
        output_parts.append("\n")

    divider = "\n\n---\n\n" if markdown else "\n\n----\n\n"
    blobs = [str(system) for system in systems] + [str(entry.phenomenon) for entry in phenomena]
    output_parts.append(divider.join(blobs))

    return "".join(output_parts)


def main():
    """
    The main entry point for the sector generation script.

    Parses command-line arguments, then generates `--num-sectors` sectors
    (1 by default), each one a fully populated `SpaceSector` (see
    `generate_sector` -- one `SystemConfig` per system via
    `build_sector_configs`, a full `StarSystem` from each, all placed in
    the sector via `SpaceSector.add_system`'s Hill-sphere-based placement)
    saved into the same database. No `galactic_center_dist_ly` is passed
    to `generate_sector` here -- this script's own CLI has no galaxy
    context (see `galaxyGen.py` for that), so every sector it produces is
    "unplaced", regardless of `--num-sectors`.

    Each sector is rendered under a single sector header with a short
    summary and an index of every system's name and star type. That
    rendering is written to `--output` (if given, all sectors appended to
    the same file, divided the same way systems within a sector are) and/or
    printed to the console (only if `--console` was given) -- by default,
    with neither flag, only a short status line per saved sector is
    printed, so a large `--num-sectors` run doesn't flood the console with
    rendered text nobody asked to see.
    """
    random.seed(secrets.randbits(128))

    args = process_args()
    divider = "\n\n---\n\n" if args.markdown else "\n\n----\n\n"

    for i in range(args.num_sectors):
        sector_name, sector = generate_sector(args)
        systems = [entry.star_system for entry in sector.entries]

        if args.output or args.console:
            output_text = render_sector_text(sector_name, systems, sector.phenomena, args.markdown)

            if args.output:
                with open(args.output, 'w' if i == 0 else 'a') as f:
                    if i > 0:
                        f.write(divider)
                    f.write(output_text)

            if args.console:
                print(output_text)

        mysql_config = _db.mysql_config_from_args(args)
        sector_id = _db.save_sector(sector, config=mysql_config)
        phenomena_note = f", {len(sector.phenomena)} phenomena" if sector.phenomena else ""
        print(f"Saved sector '{sector_name}' to the database (sector_id={sector_id}, "
              f"{len(systems)} systems{phenomena_note}, "
              f"{mysql_config.database}@{mysql_config.host}:{mysql_config.port}).")

    if args.num_sectors > 1:
        print(f"Generated {args.num_sectors} sectors.")


if __name__ == "__main__":
    main()
