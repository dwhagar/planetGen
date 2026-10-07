# planetgen/generation/run_system.py

"""
Run System
==========

The `system` command: one star system, built from `SystemConfig` (the
tri-state `+name`/`-name` options and `--system-file`). The sector
command builds each of its systems through `build_system_config`.
"""

import sys

from planetgen.db import store
from planetgen import tuning as program_constants
from planetgen.util import log
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.db.render import render_star_system
from planetgen.generation import run_common


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


def single_system_conflict(system_config):
    """
    The contradiction in a single system's resolved options, or `None`:
    `-planets` (no planets or belts at all) with anything that needs a
    planet or a belt (GEN.50). `validate_system_args` already rejects these
    on the command line; this also catches them from a `--system-file`.

    Args:
        system_config (SystemConfig): The resolved configuration.

    Returns:
        str or None: The error message.
    """
    if system_config.PLANETS is not False:
        return None
    forced = [flag for attr, flag in (("MOONS", "+moons"), ("MAX_PLANETS", "+max_planets"),
                                      ("HABITABLE_WORLD", "+habitable_world"),
                                      ("ASTEROID_BELT", "+asteroid_belt"))
              if getattr(system_config, attr) is True]
    if system_config.NUM_ORBITS:
        forced.append("num_orbits")
    if system_config.SLOTS:
        forced.append("slots")
    if not forced:
        return None
    return f"-planets cannot be combined with {', '.join(forced)}."


def run_system(args):
    """
    Generates one star system and saves it to the database, or, with
    `--output`, writes its page to a file (or stdout) and touches no
    database at all -- the admin site's one-off system page runs it that
    way (`src/planetgen/web/system_page.py`).

    A forced option the system can't meet is never saved silently
    (GEN.49): contradictory options are refused up front, and a system
    still missing a forced body after `SINGLE_SYSTEM_GENERATION_ATTEMPTS`
    whole systems ends the run with an error, saving nothing.

    Args:
        args (argparse.Namespace): Validated arguments (`command ==
            "system"`).

    Raises:
        SystemExit: On contradictory options, or when no attempt met every
            forced option.
    """
    system_config = build_system_config(args)
    conflict = single_system_conflict(system_config)
    if conflict:
        log.error(f"Error: {conflict}")
        raise SystemExit(1)

    attempts = program_constants.SINGLE_SYSTEM_GENERATION_ATTEMPTS
    for attempt in range(1, attempts + 1):
        system = StarSystem(system_config=system_config)
        if not system.unmet_requirements:
            break
        log.debug(f"System attempt {attempt}/{attempts}: no {' or '.join(system.unmet_requirements)} "
                  f"around {system.primary_star.type}; trying a new system")
    else:
        star = f"--star-type {system_config.STAR_TYPE}" if system_config.STAR_TYPE else "the drawn star"
        # GEN.81: the last try is kept, without the body it had no room for.
        run_common._refuse_or_warn(args, f"No system with {' and '.join(system.unmet_requirements)} came out of "
                              f"{attempts} tries for {star}; the last one is kept without it. Hot, "
                              f"short-lived and very large stars often have no room for one.")

    if args.output:
        text = render_star_system(system, "markdown" if system_config.MARKDOWN else "wikitext")
        if args.output == "-":
            sys.stdout.write(text if text.endswith("\n") else text + "\n")
        else:
            with open(args.output, "w", encoding="utf-8") as f:
                f.write(text)
            log.normal(f"Wrote system '{system.name}' to {args.output} (not saved to the database).")
        return

    mysql_config = store.mysql_config_from_args(args)
    star_system_id = store.save_system(system, system_config, config=mysql_config)
    run_common.RUN_COUNTS["systems"] += 1
    log.normal(f"Saved system '{system.name}' to the database (star_system_id={star_system_id}, "
               f"{mysql_config.database}@{mysql_config.host}:{mysql_config.port}).")
