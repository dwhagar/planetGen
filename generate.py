#!/usr/bin/env python
# generate.py

"""
Unified Generation CLI
=========================

One command-line entry point for every generator in this project.
`generate.py <command> [options]` replaces running `systemGen.py`,
`sectorGen.py`, `galaxyGen.py`, `galaxyPlan.py`, or `phenomenonGen.py`
directly -- every command still saves to the same MySQL database
(`stellarObjects._db`), and every command's own option surface is
exactly what its standalone script already offered:

    generate.py system [options]      -- one star system (systemGen.py)
    generate.py sector [options]      -- one or more sectors (sectorGen.py)
    generate.py galaxy [options]      -- many sectors as one galaxy (galaxyGen.py)
    generate.py plan [options]        -- the galaxy density skeleton (galaxyPlan.py)
    generate.py phenomenon [options]  -- a single exotic phenomenon (phenomenonGen.py)

Run `generate.py <command> --help` for that command's own full option
list (identical to running e.g. `sectorGen.py --help` directly).

This file adds no generation behavior of its own -- each subcommand
reuses that module's own argument-adding/validating function pair
(`add_*_arguments`/`validate_*_args`, factored out of each module's own
`process_args()` for exactly this purpose) and business logic
(`systemGen.build_system_config`, `sectorGen.generate_sector`,
`galaxyGen.run_shell_batch`/`run_local_neighborhood`/`run_random_start`,
`galaxyPlan.build_skeleton`, `phenomenonGen.generate_phenomenon`), so a
system/sector/galaxy/skeleton/phenomenon generated through `generate.py`
is built by the exact same code path as one generated through that
module's own standalone CLI. Those standalone scripts remain in place
(and stay directly importable -- `sectorGen.py` itself calls into
`systemGen.py`, `galaxyGen.py` calls into `sectorGen.py`, and the test
suite exercises several of them directly) for backward compatibility,
but `generate.py` is the single recommended entry point for actually
running generation from the command line.
"""

import argparse
import logging
import os
import random
import secrets
import sys

# stellarObjects lives at src/stellarObjects (src layout) -- add src/ to the
# import path so this keeps working without requiring `pip install .`
# first, matching every other root-level entry script.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))

import galaxyGen
import galaxyPlan
import phenomenonGen
import sectorGen
import systemGen
from stellarObjects import _db
from stellarObjects._version import VersionAction, version_banner

# Suppress transformers warnings
logging.getLogger("transformers").setLevel(logging.ERROR)


def build_parser():
    """
    Builds the top-level parser and every subcommand's own subparser, by
    calling straight into each module's `add_*_arguments` function -- so
    this can never silently drift from what each standalone script's own
    parser actually accepts.

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
               "own script of origin. Run 'generate.py <command> --help' for that command's own full "
               "option list.")
    parser.add_argument('--version', action=VersionAction, banner=version_banner('generate.py'))

    subparsers = parser.add_subparsers(dest='command', required=True)

    system_parser = subparsers.add_parser(
        'system', prefix_chars='-+',
        description="System Generation Options",
        help="Generate a single star system (see systemGen.py).")
    systemGen.add_system_arguments(system_parser)

    sector_parser = subparsers.add_parser(
        'sector', prefix_chars='-+',
        description="Sector Generation Options",
        help="Generate one or more independent sectors (see sectorGen.py).")
    sectorGen.add_shared_generation_options(sector_parser)
    sectorGen.add_sector_arguments(sector_parser)

    galaxy_parser = subparsers.add_parser(
        'galaxy', prefix_chars='-+',
        description="Galaxy Generation Options",
        help="Generate many sectors as one galaxy (see galaxyGen.py).")
    sectorGen.add_shared_generation_options(galaxy_parser)
    galaxyGen.add_galaxy_arguments(galaxy_parser)

    plan_parser = subparsers.add_parser(
        'plan',
        description="Galaxy Density Skeleton Builder",
        help="Build/replace the galaxy's density skeleton (see galaxyPlan.py).")
    galaxyPlan.add_plan_arguments(plan_parser)

    phenomenon_parser = subparsers.add_parser(
        'phenomenon',
        description="Exotic Stellar Phenomenon Generation Options",
        help="Generate a single exotic stellar phenomenon (see phenomenonGen.py).")
    phenomenonGen.add_phenomenon_arguments(phenomenon_parser)

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
    whichever subcommand was chosen through that module's own
    `validate_*_args` function(s) -- exactly the validation its
    standalone script's `process_args()` would run.

    Returns:
        argparse.Namespace: Parsed (and validated) arguments, with
            `args.command` set to the chosen subcommand name.
    """
    parser, command_parsers = build_parser()
    args = parser.parse_args()
    command_parser = command_parsers[args.command]

    if args.command == 'system':
        systemGen.validate_system_args(args, command_parser)
    elif args.command == 'sector':
        sectorGen.validate_shared_generation_args(args, command_parser)
        sectorGen.validate_sector_args(args, command_parser)
    elif args.command == 'galaxy':
        sectorGen.validate_shared_generation_args(args, command_parser)
        galaxyGen.validate_galaxy_args(args, command_parser)
    elif args.command == 'plan':
        galaxyPlan.validate_plan_args(args, command_parser)
    elif args.command == 'phenomenon':
        phenomenonGen.validate_phenomenon_args(args, command_parser)

    return args


def run_system(args):
    """
    Generates one star system and saves it to the database -- mirrors
    `systemGen.main()`'s own body exactly.

    Args:
        args (argparse.Namespace): Validated arguments (`command ==
            "system"`).
    """
    system_config, output_path = systemGen.build_system_config(args)
    system = systemGen.StarSystem(system_config=system_config)

    if output_path:
        with open(output_path, 'w') as f:
            f.write(str(system))
    else:
        print(system)

    mysql_config = _db.mysql_config_from_args(args)
    star_system_id = _db.save_system(system, system_config, config=mysql_config)
    print(f"Saved system '{system.star.name}' to the database (star_system_id={star_system_id}, "
          f"{mysql_config.database}@{mysql_config.host}:{mysql_config.port}).")


def run_sector(args):
    """
    Generates `args.num_sectors` sectors and saves each to the database --
    mirrors `sectorGen.main()`'s own body exactly.

    Args:
        args (argparse.Namespace): Validated arguments (`command ==
            "sector"`).
    """
    divider = "\n\n---\n\n" if args.markdown else "\n\n----\n\n"

    for i in range(args.num_sectors):
        sector_name, sector = sectorGen.generate_sector(args)
        systems = [entry.star_system for entry in sector.entries]

        if args.output or args.console:
            output_text = sectorGen.render_sector_text(sector_name, systems, sector.phenomena, args.markdown)

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


def run_galaxy(args):
    """
    Dispatches to shell-batch, local-neighborhood, or random-start mode --
    mirrors `galaxyGen.main()`'s own body exactly.

    Args:
        args (argparse.Namespace): Validated arguments (`command ==
            "galaxy"`).
    """
    edge_pc = galaxyGen._edge_pc()

    if args.shell is not None:
        galaxyGen.run_shell_batch(args, edge_pc)
    elif args.center_sector is not None:
        galaxyGen.run_local_neighborhood(args, edge_pc)
    else:
        galaxyGen.run_random_start(args, edge_pc)


def run_plan(args):
    """
    Builds and persists the galaxy's density skeleton -- mirrors
    `galaxyPlan.main()`'s own body exactly.

    Args:
        args (argparse.Namespace): Validated arguments (`command ==
            "plan"`).
    """
    summary = galaxyPlan.build_skeleton(args)
    print(
        f"Skeleton built in {summary['elapsed_s']:.2f}s: scanned {summary['shells_scanned']} shells, "
        f"outer edge = shell {summary['outer_shell_index']}, {summary['total_bands']} band(s) stored, "
        f"~{summary['total_candidate_slots']:,} candidate sector slots."
    )
    if not summary["edge_confirmed"]:
        print(
            f"WARNING: reached --max-shell ({args.max_shell}) without a run of "
            f"{args.empty_streak_to_stop} consecutive empty shells -- the galaxy's true edge was not "
            f"confirmed. Re-run with a larger --max-shell if these shape parameters really do produce "
            f"a galaxy this large."
        )


def run_phenomenon(args):
    """
    Generates one exotic phenomenon and saves it to the database --
    mirrors `phenomenonGen.main()`'s own body exactly.

    Args:
        args (argparse.Namespace): Validated arguments (`command ==
            "phenomenon"`).
    """
    phenomenon_type = args.type or random.choice(phenomenonGen.program_constants.PHENOMENON_TYPE_CHOICES)

    system_config = phenomenonGen.SystemConfig()
    system_config.MARKDOWN = args.markdown
    if args.num_orbits is not None:
        system_config.NUM_ORBITS = args.num_orbits

    phenomenon = phenomenonGen.generate_phenomenon(phenomenon_type, system_config, args.anchor_system, name=args.name)

    if args.output:
        with open(args.output, 'w') as f:
            f.write(str(phenomenon))
    else:
        print(phenomenon)

    mysql_config = _db.mysql_config_from_args(args)
    phenomenon_id = _db.save_phenomenon(phenomenon, system_config, phenomenon_type, config=mysql_config,
                                         sector_id=args.sector_id)
    print(f"Saved {phenomenonGen.TYPE_LABELS[phenomenon_type]} to the database (id={phenomenon_id}, "
          f"{mysql_config.database}@{mysql_config.host}:{mysql_config.port}).")


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
    validates command-line arguments, seeds the random number generator
    cryptographically (as every other generation script's own `main()`
    does), then dispatches to the chosen subcommand's own `run_*`
    function.
    """
    random.seed(secrets.randbits(128))

    args = process_args()
    _COMMAND_HANDLERS[args.command](args)


if __name__ == "__main__":
    main()
