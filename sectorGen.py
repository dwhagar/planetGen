import argparse
import logging
import random
import secrets

import systemGen
from stellarObjects import _db, program_constants
from stellarObjects._version import VersionAction, version_banner
from stellarObjects.names import SECTOR_NAMES, SECTOR_PREFIXES, SECTOR_SUFFIXES
from stellarObjects.spaceSector import SpaceSector
from stellarObjects.systemData import StarSystem
from stellarObjects.utils import generate_phoneme_salad_name

# Suppress transformers warnings
logging.getLogger("transformers").setLevel(logging.ERROR)


def process_args():
    """
    Parses command-line arguments for generating a whole sector of star
    systems at once.

    Every option `systemGen.py` accepts for tuning a single system's
    generation -- the `+name`/`-name` tri-state flags (see
    `systemGen.TRISTATE_OPTIONS`), `--star-type`, `--age`, `--markdown`, and
    the flavor-text overrides -- applies here too, uniformly, to every
    system the sector contains. `--system-file`, `--num-orbits`, and
    systemGen.py's own per-system `--name` are deliberately not offered
    here: those describe one specific, hand-crafted system (exact orbital
    slots, an exact object count, a single fixed name), which contradicts
    generating a whole sector of varied, independently-random systems. Use
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

    Returns:
        argparse.Namespace: An object containing the parsed command-line arguments.
    """
    additional_info = [
        "Additional Information:",
        "This tool generates a whole sector of independently-random star systems in one pass, reusing systemGen.py's",
        "own generation logic and options for each one. Options like +habitable_world/-habitable_world, --star-type,",
        "and --age apply uniformly to every system in the sector; use --min-habitable instead if you just want a",
        "guaranteed number of habitable systems among an otherwise varied sector. For one specific, hand-crafted",
        "system (exact orbital slots, an exact name), use systemGen.py --system-file directly instead."
    ]
    additional_info = " ".join(additional_info)

    parser = argparse.ArgumentParser(
        description="Sector Generation Options",
        epilog=additional_info,
        prefix_chars='-+')

    parser.add_argument('--version', action=VersionAction, banner=version_banner('sectorGen.py'))

    for name, _attr, description in systemGen.TRISTATE_OPTIONS:
        parser.add_argument(f'-{name}', f'+{name}', dest=name, action=systemGen.TristateAction,
                            nargs=0, default=None,
                            help=f"+{name} forces every system in the sector to have {description}; "
                                 f"-{name} forces every system in the sector to not have {description}.")

    parser.add_argument('--num-systems', type=int, default=10,
                        help="The number of star systems to generate in the sector. Defaults to 10.")
    parser.add_argument('--name', '-n', dest='sector_name', type=str,
                        help="Force the name of the sector, overriding the default random two-word name.")
    parser.add_argument('--min-habitable', type=int, default=0,
                        help="Guarantee at least this many systems in the sector have a habitable world, "
                             "chosen randomly among them, without requiring every system to have one.")

    parser.add_argument('--output', '-o', type=str, help="Output to a file.")
    parser.add_argument('--db-path', type=str,
                        help="Path to the SQLite database file the generated sector is saved to. "
                             "Defaults to stellarObjects._db.DEFAULT_DB_PATH (db/planetgen.db).")
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

    args = parser.parse_args()

    if args.num_systems < 1:
        parser.error("--num-systems must be a positive integer.")

    if args.min_habitable < 0:
        parser.error("--min-habitable cannot be negative.")

    if args.min_habitable > args.num_systems:
        parser.error("--min-habitable cannot exceed --num-systems.")

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
    first_word = generate_phoneme_salad_name(SECTOR_NAMES, SECTOR_PREFIXES, SECTOR_SUFFIXES, allow_split=False)
    second_word = generate_phoneme_salad_name(SECTOR_NAMES, SECTOR_PREFIXES, SECTOR_SUFFIXES, allow_split=False)
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
        args (argparse.Namespace): Parsed arguments from `process_args()`.

    Returns:
        list: A list of `num_systems` `SystemConfig` instances, one per
              system the sector will contain.
    """
    configs = [systemGen.build_system_config(args)[0] for _ in range(args.num_systems)]

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


def main():
    """
    The main entry point for the sector generation script.

    Parses command-line arguments, builds one `SystemConfig` per system in
    the sector (see `build_sector_configs`), generates a full `StarSystem`
    from each, places them all in a `SpaceSector` (positions auto-assigned
    by `SpaceSector.add_system`'s Hill-sphere-based placement), and renders
    them together under a single sector header with a short summary and an
    index of every system's name and star type.

    If an output file was specified, the whole rendered sector is written
    there; otherwise it's printed directly to the console. Either way, the
    generated sector is then saved to the database (`stellarObjects._db`,
    `--db-path` to override the default location) -- this is the "fill a
    sector" step: every run of this script populates the database, not
    just stdout/a file.
    """
    random.seed(secrets.randbits(128))

    args = process_args()
    sector_name = args.sector_name or generate_sector_name()

    configs = build_sector_configs(args)
    systems = [StarSystem(system_config=cfg) for cfg in configs]

    sector = SpaceSector(name=sector_name)
    for system, cfg in zip(systems, configs):
        sector.add_system(system, system_config=cfg)

    habitable_count = sum(1 for s in systems if s.hab_count > 0)

    output_parts = []
    if args.markdown:
        output_parts.append(f"# {sector_name}\n\n")
    else:
        output_parts.append(f"= {sector_name} =\n\n")

    system_word = "system" if len(systems) == 1 else "systems"
    habitable_verb = "harbors" if habitable_count == 1 else "harbor"
    output_parts.append(
        f"This sector contains {len(systems)} star {system_word}, "
        f"{habitable_count} of which {habitable_verb} a potentially habitable world.\n\n"
    )

    bullet = "-" if args.markdown else "*"
    for system in systems:
        output_parts.append(f"{bullet} {system.star.name} ({system.star.type})\n")
    output_parts.append("\n")

    divider = "\n\n---\n\n" if args.markdown else "\n\n----\n\n"
    output_parts.append(divider.join(str(system) for system in systems))

    output_text = "".join(output_parts)

    if args.output:
        with open(args.output, 'w') as f:
            f.write(output_text)
    else:
        print(output_text)

    sector_id = _db.save_sector(sector, db_path=args.db_path)
    db_path = args.db_path or _db.DEFAULT_DB_PATH
    print(f"Saved sector '{sector_name}' to the database (sector_id={sector_id}, {db_path}).")


if __name__ == "__main__":
    main()
